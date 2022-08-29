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
% Add needed paths
addpath src_kirikou/
addpath src_optiflow/

%% Actuator Disk Inputs
u_inf        =  30   ;                       % [m/s  ] Freestream Speed
rho          =  1.225;                       % [kg/m3] Fluid Density (incomp.)
d1           =  0.086;                       % [m    ] Actuator Disk Diameter
Ct           =  0.65 ;                       % [--   ] Actuator Disk Force Coefficient (C_F_a)
%% Solver Parameters
dt                 = 1e-6;              % [--   ] Artificial timestep for shape deformation (transverse)
dtx                = 0;                       % [--   ] Artificial timestep for shape deformation (longitudinal)
relax              = 1e-4;              % [--   ] Artificial timestep for vorticity adaptation
n_stances          = 60    ;                       % [--   ] Number of stances
N_iterations       = 40000 ; % round(40000 *u_inf/d1);
information_period = 400  ; % round(400   *u_inf/d1);
%% Duct Inputs
% Geometry of first duct element
airfoil_file = 'naca4412p280.air';        % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c         =   1/4        ;             % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         =   0          ;             % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct       = - 0.360*0.200+0.009;             % [m   ] Translation of airfoil rotation point from origin (x direction)
y_duct       =   0.044+0.12*0.200-0.004 ;             % [m   ] Translation of airfoil rotation point from origin (y direction)
alpha_duct   = - 2*pi()/180 ;             % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
chord_ratio  =   0.200      ;             % [--  ] Ratio for scaling airfoil chord
flip_airfoil = true         ;             % [bool] Flip airfoil upside down
% Geometry of second duct segment
airfoil_file2  = 'naca4412p280.air'; % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c2          = 1/4;                     % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
py_c2          = 0;                       % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct2        = x_duct + chord_ratio + 0.0015;                     % [m   ] Translation of airfoil rotation point from origin (x direction)
y_duct2        = 0.082;                 % [m   ] Translation of airfoil rotation point from origin (y direction)
alpha_duct2    = -30*pi()/180;            % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
chord_ratio2   = 0.035;                     % [--  ] Ratio for scaling airfoil chord
flip_airfoil2  = true;                    % [bool] Flip airfoil upside down
% Geometry of actuator:
% Inflow: Straight Free-Stream Definition
% u_inf        = 1;                       % Free-Stream Magnitude
alpha        = 0*pi/180;                  % [rad ] Set inflow angle of attack in radians (always to 0 here)
c_over_R     = 0.1;                       % [--  ] Now specify rotational effects
x_cp         = 1/4;                       % [--  ] Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;                     % [bool] Is effect of rotation accounted for ? ()
% Postprocessing Options
graph_viz    = false;                     % [bool] Make a beautiful picture of the velocity field? (setting to true takes and awful lot of time!)


%% Make duct geometry
% Geometry of first duct segment       % Geometry of second duct segment
segment1.airfoil_file  = airfoil_file  ; segment2.airfoil_file2 = airfoil_file2 ; 
segment1.px_c          = px_c          ; segment2.px_c2         = px_c2         ; 
segment1.py_c          = py_c          ; segment2.py_c2         = py_c2         ; 
segment1.x_duct        = x_duct        ; segment2.x_duct2       = x_duct2       ; 
segment1.y_duct        = y_duct        ; segment2.y_duct2       = y_duct2       ; 
segment1.alpha_duct    = alpha_duct    ; segment2.alpha_duct2   = alpha_duct2   ; 
segment1.chord_ratio   = chord_ratio   ; segment2.chord_ratio2  = chord_ratio2  ; 
segment1.flip_airfoil  = flip_airfoil  ; segment2.flip_airfoil2 = flip_airfoil2 ; 

% Now make it
[px_duct , py_duct, i_start, i_end] = make_multielement_duct(segment1, segment2);

%% Plot problem geometry
figure(1985)
for n_element = 1:length(i_start)
    plot(px_duct(i_start(n_element):i_end(n_element)), py_duct(i_start(n_element):i_end(n_element)));
    hold on; grid on;
end
plot([0 0], [-d1, d1]/2, 'k-x')

%% Solve wake without bodies
% % Create Solver Object (no lifting vortices, hence the zeros!)
DSD2 = kirikou_single_actuator_solver(Ct, 0,0,0);
% Environment settings
DSD2.d1        = d1;     % Diameter
DSD2.u_inf     = u_inf;  % Freestream
DSD2.rho       = rho;    % Density
% Solver settings
DSD2.dt        = dt       ; % [--   ] Artificial timestep for shape deformation (transverse)
DSD2.dtx       = dtx      ; % [--   ] Artificial timestep for shape deformation (longitudinal)
DSD2.relax     = relax    ; % [--   ] Artificial timestep for vorticity adaptation
DSD2.n_stances = n_stances; % [--   ] Number of stances
DSD2.N_iterations = N_iterations;
% Diagnostics
DSD2.information_period = information_period;
DSD2.n_stances          = n_stances;
% % Solve!
DSD2.preprocess_run_postprocess();


%% Create Panel Case
% Create Case Object
ipc          = inviscid_panel_case_multielement(px_duct , py_duct, i_start , i_end);
% Describe Straight Inflow
ipc.u_inf    = u_inf;
ipc.rho      = rho;
ipc.alpha    = alpha;
% Describe Still-Air Rotation
ipc.c_over_R = c_over_R;
ipc.x_cp     = x_cp;
ipc.rotation = rotation;
% Force, Solve and Postprocess
ipc.generate_solution;


%% Create Solver Object
DSD       = kirikou_single_actuator_coupled_solver(Ct, ipc);
% Environment settings
DSD.d1        = d1;     % Diameter
DSD.u_inf     = u_inf;  % Freestream
DSD.rho       = rho;    % Density
% Solver settings
DSD.dt        = dt       ; % [--   ] Artificial timestep for shape deformation (transverse)
DSD.dtx       = dtx      ; % [--   ] Artificial timestep for shape deformation (longitudinal)
DSD.relax     = relax    ; % [--   ] Artificial timestep for vorticity adaptation
DSD.n_stances = n_stances; % [--   ] Number of stances
DSD.N_iterations = N_iterations;
% Put this in monitoring mode
DSD.information_period = information_period;
DSD.n_stances          = n_stances;
% %Debug quickly
% DSD.max_RES_shape   = 1e-6;
% DSD.max_RES_stretch = 1e-6;

%% Prepare Case
% Prepocess Inputs
DSD.preprocess_inputs();
% Discretize with Straight Wake Tube as Initial Guess
DSD.define_discretization_and_initial_geometry();
% Define Reference Panel Strenghts
DSD.define_reference_panel_strenghts();
% Now formulate Initial Guesses
DSD.formulate_initial_guesses();
% Create Vortex Element Object (with state provided by initial guesses)
DSD.create_vortex_element_object();

%% Plot vortex systems

%% Solve Case
% Now Start Solver
DSD.run_solver();
% Now compute Machine Parameters
DSD.compute_machine_parameters();

%% Plot Velocity Fields
% (this takes a lot of time)
if graph_viz == true
    % Generate Velocity Field Meshes
    DSD.generate_velocity_fields();
    % Plot them with streamlines, bells and whistles!
    DSD.plot_velocity_field_with_streamlines();
    % Print
    print -dpng -r300 kirikou_example_case_17_velocity_field.png
end

%% Plot streamwise velocity along centerline
% Define Centerline
x_centerline = linspace(-2,2);
y_centerline = zeros(size(x_centerline));
% Compute Induced Speed Contributions on Centerline
[                     u_ipc_centerline , v_ipc_centerline] = DSD.ipc.induced_speed_on_many_points(x_centerline,y_centerline);
[pot_VS_centerline  , u_VS_centerline  , v_VS_centerline ] = DSD.VS.induced_speed_on_many_points(x_centerline,y_centerline);
% Compute Induced Speeds on Centerline
u_induced_centerline = u_ipc_centerline + u_VS_centerline;
v_induced_centerline = v_VS_centerline  + v_VS_centerline;
% Add Free-Stream to Induced Speeds
u_centerline = u_induced_centerline + u_inf;
v_centerline = v_induced_centerline;

% Plot Speed along centerline
figure(5)
subplot(211)
plot(x_centerline, u_centerline);
xlabel('x/c'); ylabel('u/u_inf'); grid on;
title('Streamwise Velocity along Centerline')
axis([min(x_centerline), max(x_centerline), 1/3, 4/3]);
h = subplot(212);
plot(x_centerline, v_centerline);
xlabel('x/c'); ylabel('v/u_inf'); grid on;
title('Transverse Velocity along Centerline')
h.Position = [0.1300 0.2600 0.7750 0.1912];
print -dpng -r300 kirikou_example_case_17_centerline_speed.png

% Plot Induction breakdown along centerline
figure(6)
subplot(211)
plot(x_centerline, u_induced_centerline , x_centerline, u_ipc_centerline, x_centerline, u_VS_centerline); hold on;
plot([min(x_centerline), max(x_centerline)], [0 0], 'k-.');
xlabel('x/c'); ylabel('u/u_inf'); grid on;
title('Streamwise Induction along Centerline')
axis([min(x_centerline), max(x_centerline), -2/3, 1/3 ]);
subplot(212)
plot(x_centerline, v_centerline, x_centerline, v_ipc_centerline, x_centerline, v_VS_centerline);
xlabel('x/c'); ylabel('v/u_inf'); grid on;
title('Transverse Induction along Centerline')
legend('Total Induction', 'Diffuser Induction', 'Wake Induction', 'Location', 'SouthOutside');
print -dpng -r300 kirikou_example_case_17_induction_breakdown.png


%% Plot Cp Distribution!
% % Make second reference for body
% Create Case Object
ipc2          = inviscid_panel_case_multielement(px_duct , py_duct, i_start , i_end);
% Describe Straight Inflow
ipc2.u_inf    = u_inf;
ipc2.alpha    = alpha;
ipc2.rho      = rho;
% Describe Still-Air Rotation
ipc2.c_over_R = c_over_R;
ipc2.x_cp     = x_cp;
ipc2.rotation = rotation;
% Force, Solve and Postprocess
ipc2.generate_solution;

% % Extract Results By Element (for case with wake influence)
% First Element (Top, with wake influence)
n_element = 1;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% Second Element (Bottom, with wake influence)
n_element = 2;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% Third Element (Top, second segment, with wake influence)
n_element = 3;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% Fourth Element (Bottom, second segment, with wake influence)
n_element = 4;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% % Extract Results By Element (for case with wake influence)
% First Element (Top, with wake influence)
n_element = 1;
element2(n_element).indices   = ipc2.te_first_indices(n_element):ipc.te_last_indices(n_element);
element2(n_element).px_middle = ipc2.px_middle(element(n_element).indices);
element2(n_element).py_middle = ipc2.py_middle(element(n_element).indices);
element2(n_element).cp_plot   = ipc2.cp_plot(  element(n_element).indices);
% Second Element (Bottom, with wake influence)
n_element = 2;
element2(n_element).indices   = ipc2.te_first_indices(n_element):ipc.te_last_indices(n_element);
element2(n_element).px_middle = ipc2.px_middle(element2(n_element).indices);
element2(n_element).py_middle = ipc2.py_middle(element2(n_element).indices);
element2(n_element).cp_plot   = ipc2.cp_plot(  element2(n_element).indices);
% Third Element (Top, with wake influence)
n_element = 3;
element2(n_element).indices   = ipc2.te_first_indices(n_element):ipc.te_last_indices(n_element);
element2(n_element).px_middle = ipc2.px_middle(element(n_element).indices);
element2(n_element).py_middle = ipc2.py_middle(element(n_element).indices);
element2(n_element).cp_plot   = ipc2.cp_plot(  element(n_element).indices);
% Fourth Element (Bottom, with wake influence)
n_element = 4;
element2(n_element).indices   = ipc2.te_first_indices(n_element):ipc.te_last_indices(n_element);
element2(n_element).px_middle = ipc2.px_middle(element2(n_element).indices);
element2(n_element).py_middle = ipc2.py_middle(element2(n_element).indices);
element2(n_element).cp_plot   = ipc2.cp_plot(  element2(n_element).indices);

% Plot!
figure(7)
plot(element(1).px_middle  , element(1).cp_plot, 'o-', 'Color', [0 0.4470 0.7410]); hold on;
plot(element(2).px_middle  , element(2).cp_plot, '^-', 'Color', [0 0.4470 0.7410]); 
plot(element2(1).px_middle , element2(1).cp_plot, 'o-', 'Color', [0.8500 0.3250 0.0980]); 
plot(element2(2).px_middle , element2(2).cp_plot, '^-', 'Color', [0.8500 0.3250 0.0980]); 
grid on;
xlabel('x/c'); ylabel('Cp');
title('Inviscid Cp Plot - Upstream duct segment')
legend('Element 1 - with Actuator', 'Element 2 - with Actuator', 'Element 1 - without Actuator', 'Element 2 - without Actuator');
print -dpng -r300 kirikou_example_case_17_cp_distribution.png

% Plot!
figure(17)
plot(element(3).px_middle  , element(3).cp_plot, 'o-', 'Color', [0 0.4470 0.7410]); hold on;
plot(element(4).px_middle  , element(4).cp_plot, '^-', 'Color', [0 0.4470 0.7410]); 
plot(element2(3).px_middle , element2(3).cp_plot, 'o-', 'Color', [0.8500 0.3250 0.0980]); 
plot(element2(4).px_middle , element2(4).cp_plot, '^-', 'Color', [0.8500 0.3250 0.0980]); 
grid on;
xlabel('x/c'); ylabel('Cp');
title('Inviscid Cp Plot - Downstream duct segment')
legend('Element 1 - with Actuator', 'Element 2 - with Actuator', 'Element 1 - without Actuator', 'Element 2 - without Actuator');
print -dpng -r300 kirikou_example_case_17_cp_distribution.png


%% Plot Comparison between Wake Shapes


% Plot Wake Shape
figure(8)
% With Body
plot(DSD.x_start(DSD.id_start(1):DSD.id_end(1)), DSD.y_start(DSD.id_start(1):DSD.id_end(1)), 'x-' , 'Color', [0 0.4470 0.7410]); hold on
plot(DSD.x_start(DSD.id_start(2):DSD.id_end(2)), DSD.y_start(DSD.id_start(2):DSD.id_end(2)), 'x-' , 'Color', [0 0.4470 0.7410]); 
% Without Body
plot(DSD2.x_start(DSD2.id_start(1):DSD2.id_end(1)), DSD2.y_start(DSD2.id_start(1):DSD2.id_end(1)), 'x-' , 'Color', [0.8500 0.3250 0.0980]); hold on
plot(DSD2.x_start(DSD2.id_start(2):DSD2.id_end(2)), DSD2.y_start(DSD2.id_start(2):DSD2.id_end(2)), 'x-' , 'Color', [0.8500 0.3250 0.0980]); 
% And plot body
plot(element(1).px_middle , element(1).py_middle, 'x-' , 'Color', [0      0.4470 0.7410]);
plot(element(2).px_middle , element(2).py_middle, 'x-' , 'Color', [0      0.4470 0.7410]);
plot(element(3).px_middle , element(3).py_middle, 'x-' , 'Color', [0.8500 0.3250 0.0980]);
plot(element(4).px_middle , element(4).py_middle, 'x-' , 'Color', [0.8500 0.3250 0.0980]);
% Bells and Whistles
title('Wake Shape')
legend('Upper Side - with body', 'Lower Side - with body', 'Upper Side - without body', 'Lower Side - without body', 'Location', 'East')
axis([-1 max(DSD2.x_start) -1 1])
grid on;

print -dpng -r300 kirikou_example_case_17_wake_shape.png





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