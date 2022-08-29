% IPC convergence study

%% Actuator Code Inputs:
u_inf =  1    ;                         % [m/s  ] Freestream Speed
rho   =  1    ;                         % [kg/m3] Fluid Density (incomp.)
d1    =  0.5  ;                         % [m    ] Actuator Disk Diameter
Ct    =  0.65 ;                         % [--   ] Actuator Disk Force Coefficient (C_F_a)

%% Inputs
% Geometry of first duct element
airfoil_file = 'donqio_reconstruct.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c         =  1/4;                      % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         =  0;                        % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct       = -0.20;                     % [m   ] Translation of airfoil rotation point from origin (x direction)
y_duct       =  0.68;                     % [m   ] Translation of airfoil rotation point from origin (y direction)
c_duct       =  1.00;                     % [m   ] Chord of of duct airfoil, as ratio of input file  (assumed to have chord 1)
alpha_duct   = -5*pi()/180;               % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
chord_ratio  = 1;                         % [--  ] Ratio for scaling airfoil chord
flip_airfoil = true;                      % [bool] Flip airfoil upside down
% Geometry of second duct segment
airfoil_file2  = 'airfoils/naca0012.air'; % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c2          = 1/4;                     % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
py_c2          = 0;                       % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct2        = 1.2;                     % [m   ] Translation of airfoil rotation point from origin (x direction)
y_duct2        = 3/4+0.1;                 % [m   ] Translation of airfoil rotation point from origin (y direction)
alpha_duct2    = -20*pi()/180;            % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
chord_ratio2   = 0.5;                     % [--  ] Ratio for scaling airfoil chord
flip_airfoil2  = true;                    % [bool] Flip airfoil upside down
% Geometry of actuator:
% Inflow: Straight Free-Stream Definition
% u_inf        = 1;                       % Free-Stream Magnitude
alpha        = 0*pi/180;                  % [rad ] Set angle of attack in radians
% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise) 
% (stuff that never needs to be touched, no effect)
c_over_R     = 0.1;                       % [--  ] Now specify rotational effects
x_cp         = 1/4;                       % [--  ] Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;                     % [bool] Is effect of rotation accounted for ? ()
% Postprocessing Options
graph_viz    = false;                     % [bool] Make a beautiful picture of the velocity field? (setting to true takes and awful lot of time!)


%%  Preprocessing: Construct first duct segment points
% % Load Airfoil Geometry
% coord        = load(airfoil_file);
% % Extract single airfoil coordinates
% px           = coord(:,1);
% py           = coord(:,2);
% Load Airfoil Geometry and flip it if needed
coord        = load(airfoil_file);
% Extract single airfoil coordinates
px_raw           = coord(:,1);
py_raw           = coord(:,2);
% Make actual airfoil from raw coordinates
if not(flip_airfoil)
    % Keep airfoil as is
    px           = px_raw;
    py           = py_raw;
else
    % Flip airfoil and maintain CC point ordering
    px           = flip(  px_raw);
    py           = flip( -py_raw);
end

% Make Duct Top Airfoil
%   Rotate    (+alpha around px_c,py_c)
%   Scale     (new_chord = original_chord * chord_ratio)
%   Translate (to x_duct,+y_duct)
px_top       = (  cos(alpha_duct) * (px-px_c) + sin(alpha_duct) * (py - py_c) + px_c ) * chord_ratio + x_duct;
py_top       = (- sin(alpha_duct) * (px-px_c) + cos(alpha_duct) * (py - py_c) + py_c ) * chord_ratio + y_duct;
% Make Duct Bottom Airfoil
%   Reflect   (-py over chordline)
%   Rotate    (-alpha around px_c,py_c)
%   Scale     (new_chord = original_chord * chord_ratio)
%   Translate (to x_duct,-y_duct)
px_bot       = (  cos(-alpha_duct) * (px-px_c) + sin(-alpha_duct) * (-(py - py_c)) +   px_c ) * chord_ratio + x_duct;
py_bot       = (- sin(-alpha_duct) * (px-px_c) + cos(-alpha_duct) * (-(py - py_c)) + (-py_c)) * chord_ratio - y_duct;
%   Reorder   (to maintain counterclockwise ordering after horizontal reflection)
px_bot       = px_bot(end:-1:1);
py_bot       = py_bot(end:-1:1);

%%  Preprocessing: Construct second duct segment points
% % Load Airfoil Geometry
% coord2        = load(airfoil_file2);
% % Extract single airfoil coordinates
% px2           = coord2(:,1);
% py2           = coord2(:,2);
% Load Airfoil Geometry and flip it if needed
coord2        = load(airfoil_file2);
% Extract single airfoil coordinates
px_raw2           = coord2(:,1);
py_raw2           = coord2(:,2);
% Make actual airfoil from raw coordinates
if not(flip_airfoil2)
    % Keep airfoil as is
    px2           = px_raw2;
    py2           = py_raw2;
else
    % Flip airfoil and maintain CC point ordering
    px2           = flip(  px_raw2);
    py2           = flip( -py_raw2);
end

% Make Duct Top Airfoil
%   Rotate    (+alpha around px_c,py_c)
%   Scale     (new_chord = original_chord * chord_ratio)
%   Translate (to x_duct,+y_duct)
px_top2       = (   cos(alpha_duct2) * (px2-px_c2) + sin(alpha_duct2) * (py2 - py_c2) + px_c2) * chord_ratio2 + x_duct2;
py_top2       = ( - sin(alpha_duct2) * (px2-px_c2) + cos(alpha_duct2) * (py2 - py_c2) + py_c2) * chord_ratio2 + y_duct2;
% Make Duct Bottom Airfoil
%   Reflect   (-py over chordline)
%   Rotate    (-alpha around px_c,py_c)
%   Scale     (new_chord = original_chord * chord_ratio)
%   Translate (to x_duct,-y_duct)
px_bot2       = (  cos(-alpha_duct2) * (px2-px_c2) + sin(-alpha_duct2) * (-(py2 - py_c2)) +   px_c2 ) * chord_ratio2 + x_duct2;
py_bot2       = (- sin(-alpha_duct2) * (px2-px_c2) + cos(-alpha_duct2) * (-(py2 - py_c2)) + (-py_c2)) * chord_ratio2 - y_duct2;
%   Reorder   (to maintain counterclockwise ordering after horizontal reflection)
px_bot2       = px_bot2(end:-1:1);
py_bot2       = py_bot2(end:-1:1);

%%  Preprocessing: Join everything together in a multi-element duct geometry

% Make Coordinates for Duct Geometry from airfoil geometry
px_duct      = [px_top ; px_bot; px_top2; px_bot2];
py_duct      = [py_top ; py_bot; py_top2 ; py_bot2;];

% Make Index of Trailing Edge Start and Ending points
i_start1 = 1;
i_start2 = i_start1 + length(px_top );
i_start3 = i_start2 + length(px_bot );
i_start4 = i_start3 + length(px_top2);

i_end1   =          length(px_top );
i_end2   = i_end1 + length(px_bot );
i_end3   = i_end2 + length(px_top2);
i_end4   = i_end3 + length(px_bot2);

i_start      = [i_start1 ; i_start2 ; i_start3 ; i_start4];
i_end        = [i_end1   ; i_end2   ; i_end3   ; i_end4];

% i_start      = [1              ;length(px_top)+1             ];
% i_end        = [length(px_top) ;length(px_top)+length(px_bot)];


%% Dimensional consistency tests: Single element
% Naca0012, 140 elements
ipc0_condition0          = inviscid_panel_case_multielement(px2, py2, 1, length(px2));
% Describe Straight Inflow
ipc0_condition0.u_inf    = 1    ;
ipc0_condition0.rho      = 1    ;
ipc0_condition0.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_condition0.c_over_R = c_over_R;
ipc0_condition0.x_cp     = x_cp;
ipc0_condition0.rotation = rotation;
% Force, Solve and Postprocess
ipc0_condition0.generate_solution;



ipc0_condition1          = inviscid_panel_case_multielement(px2, py2, 1, length(px2));
% Describe Straight Inflow
ipc0_condition1.u_inf    = 2       ;
ipc0_condition1.rho      = 1       ;
ipc0_condition1.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_condition1.c_over_R = c_over_R;
ipc0_condition1.x_cp     = x_cp;
ipc0_condition1.rotation = rotation;
% Force, Solve and Postprocess
ipc0_condition1.generate_solution;


ipc0_condition2          = inviscid_panel_case_multielement(px2, py2, 1, length(px2));
% Describe Straight Inflow
ipc0_condition2.u_inf    = 1       ;
ipc0_condition2.rho      = 2       ;
ipc0_condition2.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_condition2.c_over_R = c_over_R;
ipc0_condition2.x_cp     = x_cp;
ipc0_condition2.rotation = rotation;
% Force, Solve and Postprocess
ipc0_condition2.generate_solution;


num2str([ipc0_condition0.CF_y, ipc0_condition1.CF_y, ipc0_condition2.CF_y])
num2str([ipc0_condition0.F_y , ipc0_condition1.F_y , ipc0_condition2.F_y ])

%% Now refinement convergence tests: Single element

% Coarsen refined mesh twice
px2_coarsened = coarsen_even_vector(px2);
py2_coarsened = coarsen_even_vector(py2);

% Naca0012, 71 elements
ipc0_coarsened          = inviscid_panel_case_multielement(px2_coarsened , py2_coarsened, 1, length(py2_coarsened));
% Describe Straight Inflow
ipc0_coarsened.u_inf    = 1    ;
ipc0_coarsened.rho      = 1    ;
ipc0_coarsened.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_coarsened.c_over_R = c_over_R;
ipc0_coarsened.x_cp     = x_cp;
ipc0_coarsened.rotation = rotation;
% Force, Solve and Postprocess
ipc0_coarsened.generate_solution;

% Refine mesh
% % Now, make a finer mesh
px2_refined0 = px2; %refine_vector(px2_coarsened);
py2_refined0 = py2; %refine_vector(py2_coarsened);

% Naca0012, 140 elements
ipc0_refined0          = inviscid_panel_case_multielement(px2_refined0 , py2_refined0, 1, length(px2_refined0));
% Describe Straight Inflow
ipc0_refined0.u_inf    = 1    ;
ipc0_refined0.rho      = 1    ;
ipc0_refined0.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined0.c_over_R = c_over_R;
ipc0_refined0.x_cp     = x_cp;
ipc0_refined0.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined0.generate_solution;

% Refine mesh
% % Now, make a finer mesh
px2_refined1 = refine_vector(px2_refined0);
py2_refined1 = refine_vector(py2_refined0);

% Naca0012, 280 elements
ipc0_refined1          = inviscid_panel_case_multielement(px2_refined1 , py2_refined1 , 1, length(px2_refined1));
% Describe Straight Inflow
ipc0_refined1.u_inf    = 1    ;
ipc0_refined1.rho      = 1    ;
ipc0_refined1.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined1.c_over_R = c_over_R;
ipc0_refined1.x_cp     = x_cp;
ipc0_refined1.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined1.generate_solution;


% Refine mesh
% % Now, make an even finer mesh
px2_refined2 = refine_vector(px2_refined1);
py2_refined2 = refine_vector(py2_refined1);

% Naca0012, 560 elements
ipc0_refined2          = inviscid_panel_case_multielement(px2_refined2 , py2_refined2 , 1, length(px2_refined2));
% Describe Straight Inflow
ipc0_refined2.u_inf    = 1    ;
ipc0_refined2.rho      = 1    ;
ipc0_refined2.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined2.c_over_R = c_over_R;
ipc0_refined2.x_cp     = x_cp;
ipc0_refined2.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined2.generate_solution;

% Refine mesh
% % Now, make an even finer mesh
px2_refined3 = refine_vector(px2_refined2);
py2_refined3 = refine_vector(py2_refined2);

% Naca0012, 1020 elements
ipc0_refined3          = inviscid_panel_case_multielement(px2_refined3 , py2_refined3 , 1, length(px2_refined3));
% Describe Straight Inflow
ipc0_refined3.u_inf    = 1    ;
ipc0_refined3.rho      = 1    ;
ipc0_refined3.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined3.c_over_R = c_over_R;
ipc0_refined3.x_cp     = x_cp;
ipc0_refined3.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined3.generate_solution;

% Refine mesh
% % Now, make an even finer mesh
px2_refined4 = refine_vector(px2_refined3);
py2_refined4 = refine_vector(py2_refined3);

% Naca0012, 2040 elements
ipc0_refined4          = inviscid_panel_case_multielement(px2_refined4 , py2_refined4 , 1, length(px2_refined4));
% Describe Straight Inflow
ipc0_refined4.u_inf    = 1    ;
ipc0_refined4.rho      = 1    ;
ipc0_refined4.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined4.c_over_R = c_over_R;
ipc0_refined4.x_cp     = x_cp;
ipc0_refined4.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined4.generate_solution;


% Refine mesh
% % Now, make an even finer mesh
px2_refined5 = refine_vector(px2_refined4);
py2_refined5 = refine_vector(py2_refined4);

% Naca0012, 4080 elements
ipc0_refined5          = inviscid_panel_case_multielement(px2_refined5 , py2_refined5 , 1, length(px2_refined5));
% Describe Straight Inflow
ipc0_refined5.u_inf    = 1    ;
ipc0_refined5.rho      = 1    ;
ipc0_refined5.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined5.c_over_R = c_over_R;
ipc0_refined5.x_cp     = x_cp;
ipc0_refined5.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined5.generate_solution;


% Refine mesh
% % Now, make an even finer mesh
px2_refined6 = refine_vector(px2_refined5);
py2_refined6 = refine_vector(py2_refined5);

% Naca0012, 8160 elements
ipc0_refined6          = inviscid_panel_case_multielement(px2_refined6 , py2_refined6 , 1, length(px2_refined6));
% Describe Straight Inflow
ipc0_refined6.u_inf    = 1    ;
ipc0_refined6.rho      = 1    ;
ipc0_refined6.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined6.c_over_R = c_over_R;
ipc0_refined6.x_cp     = x_cp;
ipc0_refined6.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined6.generate_solution;


% Refine mesh
% % Now, make an even finer mesh
px2_refined7 = refine_vector(px2_refined6);
py2_refined7 = refine_vector(py2_refined6);

% Naca0012, 16320 elements
ipc0_refined7          = inviscid_panel_case_multielement(px2_refined7 , py2_refined7 , 1, length(px2_refined7));
% Describe Straight Inflow
ipc0_refined7.u_inf    = 1    ;
ipc0_refined7.rho      = 1    ;
ipc0_refined7.alpha    = 5*pi/180;
% Describe Still-Air Rotation
ipc0_refined7.c_over_R = c_over_R;
ipc0_refined7.x_cp     = x_cp;
ipc0_refined7.rotation = rotation;
% Force, Solve and Postprocess
ipc0_refined7.generate_solution;

N_succession   = [                70,               140,               280,               560,              1020,              2040,              4080,              8160,             16320];
F_y_succession = [ipc0_coarsened.F_y, ipc0_refined0.F_y, ipc0_refined1.F_y, ipc0_refined2.F_y, ipc0_refined3.F_y, ipc0_refined4.F_y, ipc0_refined5.F_y, ipc0_refined6.F_y, ipc0_refined7.F_y];


