%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Basic Multiobjective Optimization Example Case
%
%           Aerodynamic Goal: L/D over a range of Cls
%           Structural  Goal: Shell stiffness (Rody Kemp, under Gael)
%
%           Parallel Execution Enabled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all
 %clear all
 
%% Inputs
BPO_upper = 16;  % Number of parameters (order=degree+1) for description of airfoil upper (suction ) side shape
BPO_lower = 16;  % Number of parameters (order=degree+1) for description of airfoil lower (pressure) side shape

 
%% System Work necessary to start

% Add sources to matlab path
fs = filesep();                   % Folder separator is OS dependent
addpath([cd fs 'src_kirikou' ]);  % Add optimizer source code folder
addpath([cd fs 'src_optiflow']);  % Add optimizer source code folder


%% Create Shape and Parametrization handling objects
% Create Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , BPO_upper);
p_lower = parametrization( 'cst_lower'  , BPO_lower);
 
% Create Shape Definition Objects with Parametrization object instances
N_dummy_parameters = 0;                                                        % Number of parameters that are not shape definition ones
SD = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'cst88');    % 'cst88' is any arbitrary name you like
 
% Make Shape Fit Object for constraint suggestion , using previously
% defined Shape Definition object
SF = shape_fit_cst(SD, 'fitcst88');                            % 'fitcst88' is any arbitrary name you like

%% Import shape of Donqi airfoil
% Now load airfoil and fit it to CST parametrization
airfoil_file = 'airfoils/donqio.air';

% Load and fit airfoil 
[x, imported_coordinates, regenerated_coordinates, parameters]= get_parameters_from_file(SF , airfoil_file);

% Reconstruct coordinates with default method!
[tx, tz] = SD.generate_coordinates(200, x);

% Check adimensionalization and default shape reconstruction
figure(185);
grid on; hold on;
plot(imported_coordinates.raw.tx_coordinates, imported_coordinates.raw.tz_coordinates);
plot(imported_coordinates.adimensionalized.tx_coordinates, imported_coordinates.adimensionalized.tz_coordinates);
plot(tx, tz);
legend('Raw Coordinates', 'Adimensionalized shape', 'Reconstructed shape');

% Reconstruct coordinates with various fitting methods!
[tx_linear          , tz_linear         ] = SD.generate_coordinates(200, parameters.linear.full         );
[tx_non_linear      , tz_non_linear     ] = SD.generate_coordinates(200, parameters.non_linear.full     );
[tx_full_non_linear , tz_full_non_linear] = SD.generate_coordinates(200, parameters.full_non_linear.full);

% Compare shape reconstruction methods
figure(187);
grid on; hold on;
plot(imported_coordinates.adimensionalized.tx_coordinates, imported_coordinates.adimensionalized.tz_coordinates, 'k');
plot(tx_linear          , tz_linear         );
plot(tx_non_linear      , tz_non_linear     );
plot(tx_full_non_linear , tz_full_non_linear);
legend('Adimensional shape', 'Linear reconstruction', 'Hybrid non-linear reconstruction', 'Full non-linear reconstruction');

% Verification of shape reconstruction 
figure(189);
grid on; hold on;
plot(imported_coordinates.adimensionalized.tx_coordinates, imported_coordinates.adimensionalized.tz_coordinates, 'k');
plot(tx_full_non_linear , tz_full_non_linear);
legend('Adimensional shape', 'Full non-linear reconstruction');

%% Import shape of Donqi airfoil


% % Save coordinates
% coordinates_full_non_linear = [tx_full_non_linear , tz_full_non_linear];
% %save('donqio_reconstruct.air' , 'coordinates_full_non_linear', '-ascii', '-double');
% save('donqio_reconstruct.air' , 'coordinates_full_non_linear', '-ascii');
% 
% [tx_full_non_linear , tz_full_non_linear] = SD.generate_coordinates(200, parameters.full_non_linear.full);






