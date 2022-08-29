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
Ct           =  0.1;     % [-- ] Actuator Force Coefficient (C_F_a)

%% Inputs
% Geometry
%airfoil_file = 'airfoils/naca0012.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
airfoil_file = 'airfoils/donqio.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c         =  1/4;            % Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         =  0;              % Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct       = -0.20;           % Translation of airfoil rotation point from origin (x direction)
y_duct       =  0.68;           % Translation of airfoil rotation point from origin (y direction) (clearance at 2
c_duct       =  1.00;           % Chord of of duct airfoil, as ratio of input file  (assumed to have chord 1)
alpha_duct   = -5*pi()/180;     % Duct Airfoils Geometric Angle of Attack (in radians)
flip_airfoil = true;           % If true airrfoil is flipped upside down (boolean!)
% Inflow: Straight Free-Stream Definition
u_inf        = 1;               % Free-Stream Magnitude
alpha        = 0*pi/180;        % Set angle of attack in radians
% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise)
c_over_R     = 0.1;             % Now specify rotational effects
x_cp         = 1/4;             % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;           % Is effect of rotation accounted for ? ()

%%  Preprocessing
% Load Airfoil Geometry
coord        = load(airfoil_file);
% Extract single airfoil coordinates
px_raw           = coord(:,1);
py_raw           = coord(:,2);
% Make actual airfoil from raw coordinates
if not(flip_airfoil)
    % Keep airfoil as is
    px_air           = px_raw;
    py_air           = py_raw;
else
    % Flip airfoil and maintain CC point ordering
    px_air           = flip(  px_raw);
    py_air           = flip( -py_raw);
end

% Make Duct Top Airfoil
%   Rotate    (+alpha around px_c,py_c)
%   Translate (to x_duct,+y_duct)
px_top       = (  cos( alpha_duct) * (px_air-px_c) + sin( alpha_duct) *   (py_air - py_c)  +   px_c ) * c_duct + x_duct;
py_top       = (- sin( alpha_duct) * (px_air-px_c) + cos( alpha_duct) *   (py_air - py_c)  +   py_c ) * c_duct + y_duct;
% Make Duct Bottom Airfoil
%   Reflect   (-py over chordline)
%   Rotate    (-alpha around px_c,py_c)
%   Translate (to x_duct,-y_duct)
px_bot       = (  cos(-alpha_duct) * (px_air-px_c) + sin(-alpha_duct) * (-(py_air - py_c)) +   px_c ) * c_duct + x_duct;
py_bot       = (- sin(-alpha_duct) * (px_air-px_c) + cos(-alpha_duct) * (-(py_air - py_c)) + (-py_c)) * c_duct - y_duct;
%   Reorder   (to maintain counterclockwise ordering after horizontal reflection)
px_bot       = px_bot(end:-1:1);
py_bot       = py_bot(end:-1:1);

% Make Coordinates for Duct Geometry from airfoil geometry
px_duct      = [px_top ; px_bot];
py_duct      = [py_top ; py_bot];

% Make Index of Trailing Edge Start and Ending points
i_start      = [1 ;length(px_top)+1];
i_end        = [length(px_top) ;length(px_top)+length(px_bot)];

% Plot current setup
figure(7); plot(px_duct, py_duct, [0 0], 0.5*[-1 1]); grid on;
title('Setup Overview')


%% Create Body panel object
% Create Case Object
ipc          = inviscid_panel_case_multielement(px_duct , py_duct, i_start , i_end);
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


%% Solve Cases (about 300s per pair of cases)
% 80 cases, means 24000s, 400min, 6hours and 40min!
Ct_array = linspace(0.00001,1,8)*8/9;
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
figure(17)
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

savefig torque2018_vinit_fig17
set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'torque2018_vinit_fig17.pdf')


figure(85)
subplot(211)
plot(Ct_array,  ua_num_body_array, '.-'); hold on;
plot(Ct_array,  ua_num_free_array, '.-'); grid on;
xlabel('C_T - Thrust Coefficient of Actuator');
ylabel('U_a / U_0');
title('Average Normal Speed over Actuator')
axis([0 8/9 0.6 1.6])
legend('With Duct', 'Without Duct', 'Location', 'NorthEast')

savefig torque2018_vinit_fig85
set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'torque2018_vinit_fig85.pdf')

save(['full_data_Torque2018_alpha_duct_' , num2str(alpha_duct*180/pi(),2) , 'deg'])


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