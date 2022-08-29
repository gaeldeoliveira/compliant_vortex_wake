
% % Script for developping Donqi duct morphing technique
% % Load airfoil
% airfoil_file   = 'airfoils/donqio.air';
% % Flip airfoil upside down (boolean! true means it flips)
% flip_airfoil   = true;
%
%
% % Load Airfoil Geometry
% coord          = load(airfoil_file);
% % Extract single airfoil coordinates
% px_raw         = coord(:,1);
% py_raw         = coord(:,2);


function [px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = morph_airfoil_coordinates(px_raw, py_raw, a12, flip_airfoil)

% Make actual airfoil from raw coordinates
if not(flip_airfoil)
    % Keep airfoil as is
    px_air     = px_raw;
    py_air     = py_raw;
else
    % Flip airfoil and maintain CC point ordering
    px_air     = flip(  px_raw);
    py_air     = flip( -py_raw);
end

% Now make a set of 3rd order (2nd degree) Bernstein polynomials
b02 = @(t) 1 * t.^0 .* (1 - t).^2;
b12 = @(t) 2 * t.^1 .* (1 - t).^1;
b22 = @(t) 1 * t.^2 .* (1 - t).^0;
% Or make a general bernstein function (d = degree, k = 0...d index, t = 0 to stanwise stance
% bkd = @(t, k, d) t.^k .* (1 - t).^(d-k);

% Now make a shape function
a02   = 1.0;                               % Stays fixed to preserve LE
% a12_a = 1.0; a12_b = 0.80; a12_c = 0.60; % Given as input!
a22   = 1.0;                               % Stays fixed to preserve TE

% Shape Function
s2  = @(t, a02, a12, a22) a02*b02(t) + a12*b12(t) + a22*b22(t);

% % Now deform airfoil coordinates
%       px_air
%       py_air

% Scale coordinates to unit chord
px_air_scaled  = (px_air - min(px_air)) / (max(px_air) - min(px_air));
py_air_scaled  =  py_air                / (max(px_air) - min(px_air));

% Separate top from bottom side
[~ , i_LE]     = min(px_air_scaled);
% Get top side
px_top_scaled  = flipud(px_air_scaled(1:i_LE));
py_top_scaled  = flipud(py_air_scaled(1:i_LE));
% Get bottom side
px_bot_scaled  = px_air_scaled(i_LE:end);
py_bot_scaled  = py_air_scaled(i_LE:end);

% And morph top side
px_top_morphed = px_top_scaled;
py_top_morphed = s2(px_top_scaled, a02, a12, a22) .* py_top_scaled;


% Now compose airfoil coordinates back
% Morphed set (main result/output)
px_air_morphed = [flipud(px_top_morphed) ; px_bot_scaled(2:end)];
py_air_morphed = [flipud(py_top_morphed) ; py_bot_scaled(2:end)];
% Scaled set (secondary result/for comparison, used at intermediate step)

end

% plot(px_air_morphed, py_air_morphed);






