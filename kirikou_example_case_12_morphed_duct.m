%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Kirikou -   A simple specialized 2d Vorticity Equation Solver for     %
%               Actuator Disk Flows (Kirikou-Dogoro Suite)                %
%                                                                         %
%   Date    :   June 2014 to March 2017                                   %
%   Author  :   Gael de Oliveira                                          %
%                                                                         %
%   License :   Case by case written agreement limited to specific        %
%               applications. Distribution to any individual or           %
%               organization requires explicit written agreement from     %
%               original author.                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Environment (a bit brutally, for sure!)
clear all; close all; clc; %#ok<CLALL>


%% Actuator Code Inputs:
%       u_inf = 1     m/s       Freestream Speed
%       rho   = 1.225 kg/m3     Fluid Density (incomp.)
%       d     = 1     m         Diameter
Ct           =  6/9;     % [-- ] Actuator Force Coefficient (C_F_a)

%% Inputs
% Geometry
% Airfoil definition
airfoil_file = 'airfoils/donqio.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
flip_airfoil = true;                   % If true airfoil is flipped upside down (boolean!) (pre-morphing)
a12_morphing = -1  ;                   % Controls airfoil morphing (1= no morphing, >1= bulge outer side, <1= carve outer side, post-flipping)


% Duct: Definition
px_c         =  1/4;                   % Duct airfoil rotation point (x coordinate, in airfoil units)
py_c         =  0;                     % Duct airfoil rotation point (y coordinate, in airfoil units)
alpha_duct   = -5*pi()/180;            % Duct airfoil geometric angle of attack (in radians)
x_duct       = -0.20;                  % Translation of airfoil rotation point from origin (x direction)
y_duct       =  0.68;                  % Translation of airfoil rotation point from origin (y direction) (clearance at 2
c_duct       =  1.00;                  % Chord of of duct airfoil, as ratio of input file  (assumed to have chord 1)

% Inflow: Straight Free-Stream Definition
u_inf        = 1;                      % Free-Stream Magnitude
alpha        = 0*pi/180;               % Set angle of attack in radians
% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise)
c_over_R     = 0.1;                    % Now specify rotational effects
x_cp         = 1/4;                    % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;                  % Is effect of rotation accounted for ? ()

%%  Preprocessing
% Load Airfoil Geometry
coord        = load(airfoil_file);
% Extract single airfoil coordinates
px_raw           = coord(:,1);
py_raw           = coord(:,2);
% Flipt and morph airfoil coordinates
[px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = ...
            morph_airfoil_coordinates(px_raw, py_raw,  a12_morphing, flip_airfoil);

[px_duct_scaled  , py_duct_scaled , i_start_scaled , i_end_scaled ] = make_duct_from_airfoil(px_air_scaled , py_air_scaled ,alpha_duct, c_duct, x_duct, y_duct, px_c, py_c);
[px_duct_morphed , py_duct_morphed, i_start_morphed, i_end_morphed] = make_duct_from_airfoil(px_air_morphed, py_air_morphed,alpha_duct, c_duct, x_duct, y_duct, px_c, py_c);

% Plot current setup
figure(7);
plot(px_duct_scaled , py_duct_scaled ); hold on;
plot(px_duct_morphed, py_duct_morphed); grid on;
plot([0 0], 0.5*[-1 1]);title('Setup Overview');
legend('Original Donqi', 'Morphed Donqi', 'Actuator Disk', 'Location', 'East');
print -dpdf donqi_full_duct_morphed


%% Create Body panel object (Base scaled case)
% Create Case Object
ipc_scaled   = inviscid_panel_case_multielement(px_duct_scaled , py_duct_scaled, i_start_scaled , i_end_scaled);
% Describe Straight Inflow
ipc_scaled.u_inf    = u_inf;
ipc_scaled.alpha    = alpha;
% Describe Still-Air Rotation
ipc_scaled.c_over_R = c_over_R;
ipc_scaled.x_cp     = x_cp;
ipc_scaled.rotation = rotation;
% Force, Solve and Postprocess
ipc_scaled.generate_solution;

%% Create Body panel object (Morphed case)
% Create Case Object
ipc_morphed  = inviscid_panel_case_multielement(px_duct_morphed, py_duct_morphed, i_start_morphed, i_end_morphed);
% Describe Straight Inflow
ipc_morphed.u_inf    = u_inf;
ipc_morphed.alpha    = alpha;
% Describe Still-Air Rotation
ipc_morphed.c_over_R = c_over_R;
ipc_morphed.x_cp     = x_cp;
ipc_morphed.rotation = rotation;
% Force, Solve and Postprocess
ipc_morphed.generate_solution;

%% Create Solver Objects
% % Create Solver Object (with body)
DSD_scaled  = kirikou_single_actuator_coupled_solver(Ct, ipc_scaled);
% % Set less stringent convergence requirements for faster solution (with body)
DSD_scaled.max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD_scaled.max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)

% % Create Solver Object (with body)
DSD_morphed = kirikou_single_actuator_coupled_solver(Ct, ipc_morphed);
% % Set less stringent convergence requirements for faster solution (with body)
DSD_morphed.max_RES_stretch    = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD_morphed.max_RES_shape      = 1e-10;                        % Maximum Shape residual (volume flow accross segments)

%% Now solve
DSD_scaled.preprocess_run_postprocess;
DSD_morphed.preprocess_run_postprocess;

%% And display some comparison
disp('--- Comparison between cases ---');
% Actuator Force
disp(['Scaled  > C_F_a = ' , num2str(DSD_scaled.RES.C_F_a  )]);
disp(['Morphed > C_F_a = ' , num2str(DSD_morphed.RES.C_F_a )]);
% Body Force
disp(['Scaled  > C_F_b = ' , num2str(DSD_scaled.RES.C_F_b  )]);
disp(['Morphed > C_F_b = ' , num2str(DSD_morphed.RES.C_F_b )]);
% Power Coefficient
disp(['Scaled  > C_P   = ' , num2str(DSD_scaled.RES.Cp_num) ]);
disp(['Morphed > C_P   = ' , num2str(DSD_morphed.RES.Cp_num)]);


%%% Create Case Object for actuator less comparison
ipc_ref_scaled  = inviscid_panel_case_multielement(px_air_scaled, py_air_scaled, 1, length(py_air_scaled));
% Describe Straight Inflow
ipc_ref_scaled.u_inf    = u_inf;
ipc_ref_scaled.alpha    = 0;
% Describe Still-Air Rotation
ipc_ref_scaled.c_over_R = c_over_R;
ipc_ref_scaled.x_cp     = x_cp;
ipc_ref_scaled.rotation = rotation;
% Force, Solve and Postprocess
ipc_ref_scaled.generate_solution;

%%% Create Case Object for actuator less comparison
ipc_ref_morphed  = inviscid_panel_case_multielement(px_air_morphed, py_air_morphed, 1, length(py_air_morphed));
% Describe Straight Inflow
ipc_ref_morphed.u_inf    = u_inf;
ipc_ref_morphed.alpha    = 0;
% Describe Still-Air Rotation
ipc_ref_morphed.c_over_R = c_over_R;
ipc_ref_morphed.x_cp     = x_cp;
ipc_ref_morphed.rotation = rotation;
% Force, Solve and Postprocess
ipc_ref_morphed.generate_solution;


% % Explore the fields of the DSD object, with a particular focus on the
% DSD.RES structure, which summarizes case results!
% Typical DSD.RES values (for example inputs):
%             r_a: 0.5000               % Actuator Radius
%           phi_a: -0.4444              % Actuator Force (Loading) Density
%             x_v: -0.2000              % x-stance of lifting vortices
%             r_v: 0.4000               % y-stance of lifting vortices
%         gamma_v: 0.4000               % circulation strenght of lifting vortices
%           C_F_a: 0.8889               % Actuator Force Coefficient
%           C_F_b: 0.2088               % Vortices Streamwise Force Coefficient (computed with generalized Kutta-Joukowski-Lagally theorem using numerical velocity)
%          Cp_num: -0.7348              % Numerical Cp, computed by integrating velocity field over actuator
%         Cp_theo: -0.7318              % Theoretical Cp, computed from C_F_b using the de Vries power coefficient law (as in my Torque 2016)
%            e_Cp: 0.0041               % Absolute Difference between Numerical Cp values
%          ua_num: 0.8266               % Average velocity on actuator plane (numerical)
%         ua_theo: 0.8233               % Average velocity on actuator plane (theorethical, given C_F_a and C_F_b)
%          a1_num: 0.1734               % Induction Factor on actuator plane (numerical)
%         a1_theo: 0.1767               % Induction Factor on actuator plane (theoretical)
%     RES_stretch: 9.9787e-13           % Solution Residual Filament Vorticity Evolution (RMS)
%       RES_shape: 8.8902e-13           % Solution Residual Wake Flow Tangency Condition (RMS)
%              VS: [1x1 constant_strenght_vortex_segment_2d]