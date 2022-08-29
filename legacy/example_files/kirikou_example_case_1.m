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


% % Inputs are adimensionalized to references:
%       u_inf = 1     m/s       Freestream Speed
%       rho   = 1.225 kg/m3     Fluid Density (incomp.)
%       d     = 1     m         Diameter
Ct                  =  8/9;     % [-- ] Actuator Force Coefficient (C_F_a)
x_v                 =  0;       % [-- ] x-stance of lifting vortex pair (symmetry on x axis)
r_v                 =  0.5 + 0.05;    % [-- ] distance of lifting vortex pair to x axis (symmetry on x axis)
gamma_v             =  0.25;    % [-- ] strenght of each vortex in pair

% % Create Solver Object
DSD = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v);

% % Refine Discretization
DSD.n_stances = 60;

% % Solve!
DSD.preprocess_run_postprocess();

% % Enjoy results!
DSD.plot_velocity_field_with_streamlines();


% A more refined power coefficient
%[a, m, p] = DSD.compute_induction_on_cross_section(0, DSD.r_v-eps(), DSD.Ct/2)

%% Postprocessing (Compute Induction and Power Coefficient of Actuator) (New Method)
% Compute force on patches from force steps!
            f1_actual_patches = - cumsum(DSD.delta_f_w .* sign(DSD.y_release));
            % Reinterpolate through it with nearest neighbour interpolar
            f1_actual = @(y) interp1(DSD.y_release, f1_actual_patches , y, 'next', 0);
DSD.x1       = 0;
y_range1 = linspace(-(DSD.r1-eps), (DSD.r1-eps), 1000);
x_range1 = DSD.x1 * ones(size(y_range1));
% Compute wake effects on main extractor (all wakes!)
[~ , u_actuator1   , ~]   = induced_speed_on_many_points(DSD.VS, x_range1 , y_range1);
% Compute normal speed on actuator (only x component is needed,
% as we would take dot product of complete vector with actuator normal
% which is aligned with x axis!)
u_a1 = u_actuator1 + DSD.u_inf;
% Now compute induction factor
% a_actuator   = mean(u_inf-u_a_actuator(not(isnan(u_a_actuator))));
a1   = trapz(y_range1 , DSD.u_inf - u_a1) ./ (DSD.d1 *DSD.u_inf);
% Compute mass flow on surface
%m1   = trapz(y_range1, u_a1);
% And perform integral for power!
p1   = trapz(y_range1, f1_actual(y_range1) .* u_a1);
% Reach power coefficient!
Cp1  = DSD.p1 / (0.5*DSD.rho*DSD.S1*(DSD.u_inf^3));  


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