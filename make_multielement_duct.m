function [px_duct , py_duct, i_start, i_end] = make_multielement_duct(segment1, segment2)
% Makes multi-element duct geometry for ActiHS panel solver, to be used
% as input for ipc_panel_case_multielement class.
%
% Typical inputs:
% segment1.airfoil_file = 'donqio_reconstruct.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
% segment1.px_c         =  1/4;                      % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
% segment1.py_c         =  0;                        % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
% segment1.x_duct       = -0.20;                     % [m   ] Translation of airfoil rotation point from origin (x direction)
% segment1.y_duct       =  3+ 0.68;                  % [m   ] Translation of airfoil rotation point from origin (y direction)
% segment1.c_duct       =  1.00;                     % [m   ] Chord of of duct airfoil, as ratio of input file  (assumed to have chord 1)
% segment1.alpha_duct   = -5*pi()/180;               % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
% segment1.chord_ratio  = 1;                         % [--  ] Ratio for scaling airfoil chord
% segment1.flip_airfoil = true;                      % [bool] Flip airfoil upside down
% % Geometry of second duct segment
% segment2.airfoil_file2  = 'airfoils/naca0012.air'; % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
% segment2.px_c2          = 1/4;                     % [--  ] Airfoil Rotation Point (x coordinate, in airfoil units)
% segment2.py_c2          = 0;                       % [--  ] Airfoil Rotation Point (y coordinate, in airfoil units)
% segment2.x_duct2        = 1.2;                     % [m   ] Translation of airfoil rotation point from origin (x direction)
% segment2.y_duct2        = 3+ 3/4+0.1;              % [m   ] Translation of airfoil rotation point from origin (y direction)
% segment2.alpha_duct2    = -20*pi()/180;            % [rad ] Duct Airfoils Geometric Angle of Attack (in radians)
% segment2.chord_ratio2   = 0.5;                     % [--  ] Ratio for scaling airfoil chord
% segment2.flip_airfoil2  = true;                    % [bool] Flip airfoil upside down


% Extract inputs for first duct segment
airfoil_file  = segment1.airfoil_file ;
px_c          = segment1.px_c         ;
py_c          = segment1.py_c         ;
x_duct        = segment1.x_duct       ;
y_duct        = segment1.y_duct       ;
alpha_duct    = segment1.alpha_duct   ;
chord_ratio   = segment1.chord_ratio  ;
flip_airfoil  = segment1.flip_airfoil ;
% Extract inputs for first duct segment
airfoil_file2 = segment2.airfoil_file2;
px_c2         = segment2.px_c2        ;
py_c2         = segment2.py_c2        ;
x_duct2       = segment2.x_duct2      ;
y_duct2       = segment2.y_duct2      ;
alpha_duct2   = segment2.alpha_duct2  ;
chord_ratio2  = segment2.chord_ratio2 ;
flip_airfoil2 = segment2.flip_airfoil2;

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

end

