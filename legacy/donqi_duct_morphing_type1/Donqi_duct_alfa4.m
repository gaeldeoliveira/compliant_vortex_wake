

% Script for developping Donqi duct morphing technique
% Load airfoil
airfoil_file   = 'airfoils/donqio.air';
% Flip airfoil upside down (boolean! true means it flips)
flip_airfoil   = true;           
% Morphing factor
a12           = 0.5;

% Load Airfoil Geometry
coord          = load(airfoil_file);
% Extract single airfoil coordinates
px_raw         = coord(:,1);
py_raw         = coord(:,2);

% Morph coordinates
[px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = ...
            morph_airfoil_coordinates(px_raw, py_raw,  a12, flip_airfoil);
[px_air_morphed_a, py_air_morphed_a] = morph_airfoil_coordinates(px_raw, py_raw,  1, flip_airfoil);
[px_air_morphed_b, py_air_morphed_b] = morph_airfoil_coordinates(px_raw, py_raw,  0, flip_airfoil);
[px_air_morphed_c, py_air_morphed_c] = morph_airfoil_coordinates(px_raw, py_raw, -1, flip_airfoil);
[px_air_morphed_d, py_air_morphed_d] = morph_airfoil_coordinates(px_raw, py_raw, -2, flip_airfoil);
[px_air_morphed_e, py_air_morphed_e] = morph_airfoil_coordinates(px_raw, py_raw, -3, flip_airfoil);
[px_air_morphed_f, py_air_morphed_f] = morph_airfoil_coordinates(px_raw, py_raw, -4, flip_airfoil);


% Plot morphing results!
figure(1)
plot(px_air_scaled   , py_air_scaled    ,   '.-'); hold on   ;
plot(px_air_morphed_a, py_air_morphed_a); axis equal; grid on ;
plot(px_air_morphed_b, py_air_morphed_b); axis equal; grid on ;
plot(px_air_morphed_c, py_air_morphed_c); axis equal; grid on ;
plot(px_air_morphed_d, py_air_morphed_d); axis equal; grid on ;
plot(px_air_morphed_e, py_air_morphed_e); axis equal; grid on ;
plot(px_air_morphed_f, py_air_morphed_f); axis equal; grid on ;
xlabel('x/c'); ylabel('y/c'); title('Morphing of Donqi duct shape');
legend('Original shape (scaled)', ...
        'Morphed shape a_{12} =  1', ...
        'Morphed shape a_{12} =  0', ...
        'Morphed shape a_{12} = -1', ...
        'Morphed shape a_{12} = -2', ...
        'Morphed shape a_{12} = -3', ...
        'Morphed shape a_{12} = -4'  );
print -dpdf donqi_duct_morphing.pdf










