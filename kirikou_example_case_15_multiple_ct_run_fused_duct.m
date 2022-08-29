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
Ct           =  0.1;     % [-- ] Actuator Force Coefficient (C_F_a) (in this script: only used for initialization)

%% Inputs
% Geometry
% Airfoil:  definition
airfoil_file = 'airfoils/donqio.air';   % Filename of reference airfoil
BPO_upper    = 16  ;                    % Number of parameters (order=degree+1) for description of airfoil upper (suction ) side shape
BPO_lower    = 16  ;                    % Number of parameters (order=degree+1) for description of airfoil lower (pressure) side shape
eta_airfoil  = -.2 ;                    % Fusing parameter 
flip_airfoil = true;                    % % If true airfoil is flipped upside down (boolean!) (post-morphing)

% Duct: layout definition
px_c         =  1/4;            % Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         =  0;              % Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct       = -0.20;           % Translation of airfoil rotation point from origin (x direction)
y_duct       =  0.68;           % Translation of airfoil rotation point from origin (y direction) (clearance at 2
c_duct       =  1.00;           % Chord of of duct airfoil, as ratio of input file  (assumed to have chord 1)
alpha_duct   = -5*pi()/180;     % Duct Airfoils Geometric Angle of Attack (in radians)

% Inflow: Straight Free-Stream Definition
u_inf        = 1;               % Free-Stream Magnitude
alpha        = 0*pi/180;        % Set angle of attack in radians

% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise)
c_over_R     = 0.1;             % Now specify rotational effects
x_cp         = 1/4;             % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;           % Is effect of rotation accounted for ? ()

% Set of cases: Ct_ranges
Ct_min       = 0.0001  ;        % Minimum Ct of actuator (positive and negative both work, but solution is singular for exactly 0)
Ct_max       = 8/9     ;        % Maximum Ct of actuator (cannot exceed one, convergence degrades when Ct>8/9, still doable (slowly) around 0.91-0.93 with a lot of under-relaxation)
n_Ct_cases   = 8       ;        % Number of different actuator Cts between Ct_min and Ct_max for which solution is generated (about 300s-900s per pair of cases, depends on Ct, longer near and above Ct=8/9)


%% Add sources to matlab path
fs = filesep();                   % Folder separator is OS dependent
addpath([cd fs 'src_kirikou' ]);  % Add optimizer source code folder
addpath([cd fs 'src_optiflow']);  % Add optimizer source code folder

%%  Preprocessing
% Scale, morph and flip airfoil coordinates
[px_air_morphed, py_air_morphed] = ...
            fuse_airfoil_coordinates(BPO_upper, BPO_lower, airfoil_file, eta_airfoil, flip_airfoil);
% Construct duct from airfoil coordinates
[px_duct_morphed , py_duct_morphed, i_start_morphed, i_end_morphed] = make_duct_from_airfoil(px_air_morphed, py_air_morphed,alpha_duct, c_duct, x_duct, y_duct, px_c, py_c);

% Plot current setup
figure(301); 
plot(px_duct_morphed, py_duct_morphed); grid on; hold on;
plot([0 0], 0.5*[-1 1]);title('Setup Overview');
legend(['Fused Donqi duct (IntQi, eta=', num2str(eta_airfoil), ')'], 'Actuator Disk', 'Location', 'East');
print -dpdf fig/301


%% Create body panel object
% Create Case Object
ipc          = inviscid_panel_case_multielement(px_duct_morphed , py_duct_morphed, i_start_morphed , i_end_morphed);
% Describe Straight Inflow
ipc.u_inf    = u_inf;
ipc.alpha    = alpha;
% Describe Still-Air Rotation
ipc.c_over_R = c_over_R;
ipc.x_cp     = x_cp;
ipc.rotation = rotation;
% Force, Solve and Postprocess
ipc.generate_solution;


%% Create Solver Objects
% % Create Solver Object (with body)
DSD_body = kirikou_single_actuator_coupled_solver(Ct, ipc);
% % Create Solver Object (no body, no lifting vortices - hence the zeros!)
DSD_free = kirikou_single_actuator_solver(Ct, 0,0,0);
% % Set less stringent convergence requirements for faster solution (with body)
DSD_body.max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD_body.max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
% % Set less stringent convergence requirements for faster solution (no body, no lifting vortices - hence the zeros!)
DSD_free.max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD_free.max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)


%% Solve Cases (about 300s-900s per pair of cases, depends on Ct)
Ct_array = linspace(Ct_min , Ct_max , n_Ct_cases);
RES_body_cell = cell(size(Ct_array));
RES_free_cell = cell(size(Ct_array));
tic
for n_Ct = 1:length(Ct_array)
    % % Set Ct
    DSD_body.Ct = Ct_array(n_Ct);
    DSD_free.Ct = Ct_array(n_Ct);
    % % Solve!
    DSD_body.preprocess_run_postprocess();
    DSD_free.preprocess_run_postprocess();
    % % Store results!
    RES_body_cell{n_Ct} = DSD_body.RES;
    RES_free_cell{n_Ct} = DSD_free.RES;
end
toc


%% Get cell fields into arrays
% Preallocate
C_F_b_body_array  = zeros(size(Ct_array));

Cp_num_body_array = zeros(size(Ct_array));
Cp_num_free_array = zeros(size(Ct_array));

ua_num_body_array = zeros(size(Ct_array));
ua_num_free_array = zeros(size(Ct_array));

% Fill
for n_Ct = 1:length(Ct_array)
    C_F_b_body_array( n_Ct) = RES_body_cell{n_Ct}.C_F_b;
    Cp_num_body_array(n_Ct) = RES_body_cell{n_Ct}.Cp_num;
    Cp_num_free_array(n_Ct) = RES_free_cell{n_Ct}.Cp_num;
    
    ua_num_body_array(n_Ct) = RES_body_cell{n_Ct}.ua_num;
    ua_num_free_array(n_Ct) = RES_free_cell{n_Ct}.ua_num;
end

%% Plot
figure(301)
subplot(211)
plot(Ct_array, C_F_b_body_array, '.-'); grid on;
xlabel('C_T - Thrust Coefficient of Actuator');
ylabel('C_X');
title('Body Force Coefficient')
legend('With Duct', 'Location', 'SouthEast')
axis([0 8/9 0 0.25])
subplot(212)
plot(Ct_array, -Cp_num_body_array, '.-'); hold on;
plot(Ct_array, -Cp_num_free_array, '.-'); grid on;
xlabel('C_T - Thrust Coefficient of Actuator');
ylabel('C_P ');
title('Power Coefficient');
legend('With Duct', 'Without Duct', 'Location', 'SouthEast')
axis([0 8/9 0 0.8])

savefig fig/301
set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'fig/301.pdf')


figure(303)
subplot(211)
plot(Ct_array,  ua_num_body_array, '.-'); hold on;
plot(Ct_array,  ua_num_free_array, '.-'); grid on;
xlabel('C_T - Thrust Coefficient of Actuator');
ylabel('U_a / U_0');
title('Average Normal Speed over Actuator')
axis([0 8/9 0.6 1.6])
legend('With Duct', 'Without Duct', 'Location', 'NorthEast')

savefig fig/303
set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'fig/303.pdf')

save(['full_data_case_15_alpha_duct_' , num2str(alpha_duct*180/pi(),2) , 'deg__eta_', num2str(eta_airfoil)])


% %% Plotting
% % % Extract Results By Element (for case with wake influence)
% % First Element (Top, with wake influence)
% n_element = 1;
% element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
% element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
% element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
% element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% % Second Element (Bottom, with wake influence)
% n_element = 2;
% element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
% element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
% element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
% element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% 
% % Plot Wake Shape
% figure(8)
% % With Body
% plot(DSD_body.x_start(DSD_body.id_start(1):DSD_body.id_end(1)), DSD_body.y_start(DSD_body.id_start(1):DSD_body.id_end(1)), 'x-' , 'Color', [0 0.4470 0.7410]); hold on
% plot(DSD_body.x_start(DSD_body.id_start(2):DSD_body.id_end(2)), DSD_body.y_start(DSD_body.id_start(2):DSD_body.id_end(2)), 'x-' , 'Color', [0 0.4470 0.7410]); 
% % Without Body
% plot(DSD_free.x_start(DSD_free.id_start(1):DSD_free.id_end(1)), DSD_free.y_start(DSD_free.id_start(1):DSD_free.id_end(1)), 'x-' , 'Color', [0.8500 0.3250 0.0980]); hold on
% plot(DSD_free.x_start(DSD_free.id_start(2):DSD_free.id_end(2)), DSD_free.y_start(DSD_free.id_start(2):DSD_free.id_end(2)), 'x-' , 'Color', [0.8500 0.3250 0.0980]); 
% % And plot body
% plot(element(1).px_middle , element(1).py_middle, 'x-' , 'Color', [0 0.4470 0.7410]);
% plot(element(2).px_middle , element(2).py_middle, 'x-' , 'Color', [0 0.4470 0.7410]);
% % And plot actuator
% plot([0 0], 0.5*[-1 1], 'k')
% 
% % Bells and Whistles
% title('Wake Shape')
% legend('Upper Side - with body', 'Lower Side - with body', 'Upper Side - without body', 'Lower Side - without body', 'Location', 'East')
% %axis([-1 max(DSD_free.x_start) -1 1])
% axis([-1 6 -1 1])
% grid on;

%print -dpng -r300 kirikou_example_case_10_wake_shape.png





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