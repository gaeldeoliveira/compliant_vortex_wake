function [px_duct , py_duct, i_start , i_end] = make_duct_from_airfoil(px_air , py_air ,alpha_duct, c_duct, x_duct, y_duct, px_c, py_c)
% [px_duct , py_duct, i_start , i_end] = make_duct_from_airfoil(alpha_duct, c_duct, x_duct, y_duct, px_c, py_c)
%   Make duct as multi-element collection of airfoils, given:
%       px_air, py_air          - airfoil chance
%       alpha_duct              - duct angle in radians
%       c_duct                  - duct chord
%       x_duct                  - duct streamwise    translation
%       y_duct                  - duct cross-stream  translation
%       py_c                    - airfoil rotation point (x over chord)
%       px_c                    - airfoil rotation point (y over chord)
%   Output a duct:
%       px_duct, py_duct        - collection of coordinates 
%       i_start, i_end          - segmentation indices


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

end