%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   File    : Compute DAWT flow with for Donqi duct and fused (2nd        %
%             generation morphing) derivatives                            %
%   Purpose : Study effect of camber on duct performance                  %
%   Authors : Gael de Oliveira and Vinit Dighe                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Kirikou -   A simple specialized 2d Vorticity Equation Solver for     %
%               Actuator Disk Flows (Kirikou-Dogoro Suite)                %
%                                                                         %
%   Date    :   June 2014 to March 2017                                   %
%   Author  :   Gael de Oliveira                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   OptiFLOW - Parallel MultiObjective Airfoil Optimization System        %
%                                                                         %
%   Date    :   June 2014 to March 2017                                   %
%   Author  :   Gael de Oliveira                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all

% Clean Environment (a bit brutally, for sure!)
clear all; close all; clc; %#ok<CLALL>


%% Actuator Code Inputs:
%       u_inf = 1     m/s       Freestream Speed
%       rho   = 1.225 kg/m3     Fluid Density (incomp.)
%       d     = 1     m         Diameter
Ct           =   -8/9; % -0.005;     % [-- ] Actuator Force Coefficient (C_F_a)
x_v          =   0  ;
r_v          =   1.0;
gamma_v      =   0  ;


%% Add sources to matlab path
fs = filesep();                   % Folder separator is OS dependent
addpath([cd fs 'src_kirikou' ]);  % Add optimizer source code folder
addpath([cd fs 'src_optiflow']);  % Add optimizer source code folder

%% Create Planar Solver Object
% % Create Solver Object (no body)
DSD  = kirikou_single_actuator_solver_hover(Ct, x_v, r_v, gamma_v);
% % Set less stringent convergence requirements for faster solution (with body)
DSD.max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD.max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
DSD.n_stances = 120;

%% Solve Planar Flow Model (for u_inf=1)
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
% Now Start Solver
DSD.run_solver();
% Now compute Machine Parameters
DSD.compute_machine_parameters();
% Generate Velocity Fields
DSD.generate_velocity_fields();

%% Create Planar Solver Object (for Hover)
% % Create Solver Object (no body)
DSD_hover  = kirikou_single_actuator_solver_hover(Ct, x_v, r_v, gamma_v);
% % Set less stringent convergence requirements for faster solution (with body)
DSD_hover.max_RES_stretch     = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD_hover.max_RES_shape       = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
DSD_hover.n_stances = 120;

%% Solve Planar Flow Model (for Hover)
% Prepocess Inputs
DSD_hover.preprocess_inputs();
% Discretize with Straight Wake Tube as Initial Guess
DSD_hover.define_discretization_and_initial_geometry();
% Define Reference Panel Strenghts
DSD_hover.define_reference_panel_strenghts();
% Now formulate Initial Guesses
DSD_hover.formulate_initial_guesses();
% Create Vortex Element Object (with state provided by initial guesses)
DSD_hover.create_vortex_element_object();
% Now set speed to hovering
DSD_hover.u_inf = 0;
% Now Start Solver
DSD_hover.run_solver();
% And restore speed reference
DSD_hover.u_inf = 1;
% Now compute Machine Parameters
DSD_hover.compute_machine_parameters();
% Set speed back to hovering
DSD_hover.u_inf = 0;
% Generate Velocity Fields
DSD_hover.generate_velocity_fields();

% Plot Wake
plot(DSD.x_end , DSD.y_end, '.-'); axis([-1 5 -1 1]); grid on;

% Plot Streamlines
close all
n_lines = 10; streamline(DSD.x_mesh,DSD.y_mesh, DSD.u_mesh, DSD.v_mesh, 0 * ones(n_lines,1), 0.6 * linspace(-1,1, n_lines))
hold on
plot([0 0], 0.5*[-1 1], 'r.-')
n_lines = 10; streamline(DSD.x_mesh,DSD.y_mesh, -DSD.u_mesh, -DSD.v_mesh, 0 * ones(n_lines,1), 0.6 * linspace(-1,1, n_lines))

%% Solve Planar Flow Model (for Hover, inhomogeneous actuator loading)
DSD  = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v);
% % Set less stringent convergence requirements for faster solution (with body)
DSD.max_RES_stretch      = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD.max_RES_shape        = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
DSD.n_filaments_per_side = 10;
DSD.n_stances = 20;
% Prepocess Inputs
DSD.preprocess_inputs();
% Make loading inhomogeneous
DSD.DL = disk_loading(DSD.f1_average);
DSD.DL.ab2 = -1; % (tend to zero at tip)

% Discretize with Straight Wake Tube as Initial Guess
DSD.define_discretization_and_initial_geometry();
% Define Reference Panel Strenghts
DSD.define_reference_panel_strenghts();
% Now formulate Initial Guesses
DSD.formulate_initial_guesses();
% Create Vortex Element Object (with state provided by initial guesses)
DSD.create_vortex_element_object();
% Now set speed to hovering
DSD.u_inf = 0.1;
% Now Start Solver
DSD.run_solver();
% And restore speed reference
DSD.u_inf = 1;
% Now compute Machine Parameters
DSD.compute_machine_parameters();
% Set speed back to hovering
DSD.u_inf = 0.1;
% Generate Velocity Fields
DSD.generate_velocity_fields();


plot(DSD.x_end , DSD.y_end, '.-'); axis([-1 5 -1 1]); grid on;
y_over_r = linspace(0,1); plot(y_over_r, DSD.DL.f_w_function(y_over_r))

%% Solve Planar Flow Model (for Hover, inhomogeneous actuator loading)
DSD  = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v);
% % Set less stringent convergence requirements for faster solution (with body)
DSD.max_RES_stretch      = 1e-10;                        % Maximum Stretching (extension) residual (vorticity conservation)
DSD.max_RES_shape        = 1e-10;                        % Maximum Shape residual (volume flow accross segments)
DSD.n_filaments_per_side = 10;
DSD.n_stances            = 20;
% Prepocess Inputs
DSD.preprocess_inputs();
% Make loading inhomogeneous
DSD.DL                   = disk_loading(DSD.f1_average);
DSD.DL.ab2               = -0.5; % (tend to zero at tip)

% Discretize with Straight Wake Tube as Initial Guess
DSD.define_discretization_and_initial_geometry();
% Define Reference Panel Strenghts
DSD.define_reference_panel_strenghts();
% Now formulate Initial Guesses
DSD.formulate_initial_guesses();
% Create Vortex Element Object (with state provided by initial guesses)
DSD.create_vortex_element_object();
% Now set speed to hovering
DSD.u_inf = 0;
% Make relaxation slower
DSD.dt    = 0.0001;
DSD.relax = 0.0001;
% Now Start Solver
DSD.run_solver();
% And go on another time
DSD.run_solver();
% And restore speed reference
DSD.u_inf = 1;
% Now compute Machine Parameters
DSD.compute_machine_parameters();
% Set speed back to hovering
DSD.u_inf = 0;
% Generate Velocity Fields
DSD.generate_velocity_fields();


plot(DSD.x_end , DSD.y_end, '.-'); axis([-1 5 -1 1]); grid on;
y_over_r = linspace(0,1); plot(y_over_r, DSD.DL.f_w_function(y_over_r))
