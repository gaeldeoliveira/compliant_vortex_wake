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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   File    : Fusing method (2nd generation morphing) of airfoil shapes   %
%             for Donqi duct!                                             %
%   Purpose : Method development                                          %
%   Authors : Gael de Oliveira and Vinit Dighe                            %
%                                                                         %
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
plot(tx, tz, '--');
legend('Raw Coordinates', 'Adimensionalized shape', 'Reconstructed shape');
title(['DonQi airfoil import/scaling (\theta ='  , num2str(imported_coordinates.processing_parameters.theta), ...
                         '  \Delta x =', num2str(imported_coordinates.processing_parameters.translation_factor), ...
                         '   scale ='  , num2str(imported_coordinates.processing_parameters.scale_factor), ' )']);

% Reconstruct coordinates with various fitting methods!
[tx_linear          , tz_linear         ] = SD.generate_coordinates(200, parameters.linear.full         );
[tx_non_linear      , tz_non_linear     ] = SD.generate_coordinates(200, parameters.non_linear.full     );
[tx_full_non_linear , tz_full_non_linear] = SD.generate_coordinates(200, parameters.full_non_linear.full);
print -dpdf fig/185.pdf

% Compare shape reconstruction methods
figure(187);
grid on; hold on;
plot(imported_coordinates.adimensionalized.tx_coordinates, imported_coordinates.adimensionalized.tz_coordinates, 'k');
plot(tx_linear          , tz_linear         );
plot(tx_non_linear      , tz_non_linear     );
plot(tx_full_non_linear , tz_full_non_linear);
legend('Adimensional shape', 'Linear reconstruction', 'Hybrid non-linear reconstruction', 'Full non-linear reconstruction');
title(['Opt for fully non-linear polynomial fit (BPO_upper = ' , num2str(BPO_upper) , ...
                                              '  BPO_lower = ' , num2str(BPO_lower)]);
print -dpdf fig/187.pdf
% Verification of shape reconstruction 
figure(189);
grid on; hold on;
plot(imported_coordinates.adimensionalized.tx_coordinates, imported_coordinates.adimensionalized.tz_coordinates);
plot(tx_full_non_linear , tz_full_non_linear, '--');
legend('Adimensional shape', 'Full non-linear reconstruction');
title(['Opt for fully non-linear polynomial fit (BPO_upper = ' , num2str(BPO_upper) , ...
                                              '  BPO_lower = ' , num2str(BPO_lower)]);
print -dpdf fig/189.pdf

%% Import Donqi duct, make symmetric and flat variants

% % Select reference polynomial description for donqi airfoil
donqi_parameters_full  = parameters.full_non_linear.full;
donqi_parameters_upper = parameters.full_non_linear.upper;
donqi_parameters_lower = parameters.full_non_linear.lower;
% Identify trailing edge thickness
donqi_parameters_zte   = donqi_parameters_full(BPO_upper+BPO_lower+1);
% Donqi airfoil parameter rebuild (== donqi_parameters_full)
donqi_parameters_rfull = [donqi_parameters_upper, donqi_parameters_lower, donqi_parameters_zte];

% % Construct polynomial description of a symmetric version of donqi airfoil
symqi_parameters_upper = donqi_parameters_upper;
symqi_parameters_lower = donqi_parameters_upper;
% Set trailing edge thickness ( to may origignal donqi airfoil)
symqi_parameters_zte   = donqi_parameters_zte;
% Symqi airfoil parameter full vector build 
symqi_parameters_rfull = [symqi_parameters_upper, symqi_parameters_lower, symqi_parameters_zte];

% % Construct polynomial description of a flat plate version of donqi airfoil
fltqi_parameters_upper = donqi_parameters_upper;
fltqi_parameters_lower = zeros(size(donqi_parameters_lower));
fltqi_parameters_lower(1    ) =   donqi_parameters_lower(1);
fltqi_parameters_lower(2:end) = - donqi_parameters_upper(2:end);
% Set trailing edge thickness ( to may origignal donqi airfoil)
fltqi_parameters_zte   = donqi_parameters_zte;
% Symqi airfoil parameter full vector build 
fltqi_parameters_rfull = [fltqi_parameters_upper, fltqi_parameters_lower, fltqi_parameters_zte];

% % Regenerate coordinates for all three variants
[tx_donqi , tz_donqi] = SD.generate_coordinates(200, donqi_parameters_rfull);
[tx_symqi , tz_symqi] = SD.generate_coordinates(200, symqi_parameters_rfull);
[tx_fltqi , tz_fltqi] = SD.generate_coordinates(200, fltqi_parameters_rfull);
% Plot them
figure(191); 
subplot(311); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_donqi , tz_donqi); title('Donqi duct: Original                 ');
subplot(312); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_symqi , tz_symqi); title('Donqi duct: Symmetric variant');
subplot(313); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_fltqi , tz_fltqi); title('Donqi duct: Thin plate variant ');
print -dpdf fig/191.pdf
% % Write coordinates of all three variants to file!
% Save donqi coordinates
donqi_coordinates     = [tx_donqi , tz_donqi];
save('donqi_reconstructed.air' , 'donqi_coordinates', '-ascii');
% Save symqi coordinates
symqi_coordinates     = [tx_symqi , tz_symqi];
save('symqi_reconstructed.air' , 'symqi_coordinates', '-ascii');
% Save fltqi coordinates
fltqi_coordinates     = [tx_fltqi , tz_fltqi];
save('fltqi_reconstructed.air' , 'fltqi_coordinates', '-ascii');

%% Make interpolating polynomials
% Parabolic interpolant in eta
%   eta = -1 => p=fltqi_parameters_rfull;
%   eta =  0 => p=donqi_parameters_rfull;
%   eta =  1 => p=symqi_parameters_rfull;

% % Make (meta)polynomial coefficients
a_i = 0.5 * (symqi_parameters_rfull + fltqi_parameters_rfull) - donqi_parameters_rfull;
b_i = 0.5 * (symqi_parameters_rfull - fltqi_parameters_rfull)                         ;
c_i =                                                           donqi_parameters_rfull;

% % Evaluate (meta)polynomial to get intermediate airfoil shape parameters
% Make evaluator function (not vectorized, no need)
intqi_parameters_fun = @(eta) a_i * eta^2 + b_i * eta + c_i;
% Test evaluator (all should be true! to a few eps!)
fltqi_parameters_rfull == intqi_parameters_fun(-1) %#ok<NOPTS,EQEFF>
donqi_parameters_rfull == intqi_parameters_fun( 0) %#ok<NOPTS,EQEFF>
symqi_parameters_rfull == intqi_parameters_fun( 1) %#ok<NOPTS,EQEFF>

% % And generate some intermediate airfoils 
% Getting thicker
intqi_parameters_rfull_p02 = intqi_parameters_fun( 0.2);
intqi_parameters_rfull_p04 = intqi_parameters_fun( 0.4);
intqi_parameters_rfull_p06 = intqi_parameters_fun( 0.6);
intqi_parameters_rfull_p08 = intqi_parameters_fun( 0.8);
% Getting thinner
intqi_parameters_rfull_m02 = intqi_parameters_fun(-0.2);
intqi_parameters_rfull_m04 = intqi_parameters_fun(-0.4);
intqi_parameters_rfull_m06 = intqi_parameters_fun(-0.6);
intqi_parameters_rfull_m08 = intqi_parameters_fun(-0.8);

% % Make corresponding coordinates
% Getting thicker
[tx_intqi_p02 , tz_intqi_p02] = SD.generate_coordinates(200, intqi_parameters_rfull_p02);
[tx_intqi_p04 , tz_intqi_p04] = SD.generate_coordinates(200, intqi_parameters_rfull_p04);
[tx_intqi_p06 , tz_intqi_p06] = SD.generate_coordinates(200, intqi_parameters_rfull_p06);
[tx_intqi_p08 , tz_intqi_p08] = SD.generate_coordinates(200, intqi_parameters_rfull_p08);
% Getting thinner
[tx_intqi_m02 , tz_intqi_m02] = SD.generate_coordinates(200, intqi_parameters_rfull_m02);
[tx_intqi_m04 , tz_intqi_m04] = SD.generate_coordinates(200, intqi_parameters_rfull_m04);
[tx_intqi_m06 , tz_intqi_m06] = SD.generate_coordinates(200, intqi_parameters_rfull_m06);
[tx_intqi_m08 , tz_intqi_m08] = SD.generate_coordinates(200, intqi_parameters_rfull_m08);

% % And plot
% Getting thicker
figure(193); 
grid on; hold on; axis equal;
plot(tx_donqi     , tz_donqi    , 'k');
plot(tx_intqi_p02 , tz_intqi_p02);
plot(tx_intqi_p04 , tz_intqi_p04);
plot(tx_intqi_p06 , tz_intqi_p06);
plot(tx_intqi_p08 , tz_intqi_p08);
plot(tx_symqi     , tz_symqi    , 'k');
legend('Original Donqi', '\eta = 0.2',  '\eta = 0.4', '\eta = 0.6', '\eta = 0.8', 'Symqi (symmetric Donqi)')
print -dpdf fig/193.pdf
% Getting thinner
figure(195); 
grid on; hold on; axis equal;
plot(tx_donqi     , tz_donqi    , 'k');
plot(tx_intqi_m02 , tz_intqi_m02);
plot(tx_intqi_m04 , tz_intqi_m04);
plot(tx_intqi_m06 , tz_intqi_m06);
plot(tx_intqi_m08 , tz_intqi_m08);
plot(tx_fltqi     , tz_fltqi    , 'k');
legend('Original Donqi', '\eta =-0.2',  '\eta =-0.4', '\eta =-0.6', '\eta =-0.8', 'Fltqi (flat Donqi)')
print -dpdf fig/195.pdf

% % Make second round of plots
% Getting thicker
figure(197); 
subplot(311); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_donqi , tz_donqi); title('Donqi duct: Original                 ');
subplot(312); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_intqi_p04 , tz_intqi_p04); title('Intqi duct: \eta = 0.4                ');
subplot(313); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_symqi , tz_symqi); title('Symqi duct: Symmetric variant');
print -dpdf fig/197.pdf
% Getting thicker
figure(199); 
subplot(311); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_donqi , tz_donqi); title('Donqi duct: Original                 ');
subplot(312); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_intqi_m04 , tz_intqi_m04); title('Intqi duct: \eta =-0.4                ');
subplot(313); grid on; hold on; axis equal; axis([-.25 1.25 -.16 .16])
plot(tx_fltqi , tz_fltqi); title('Donqi duct: Thin plate variant ');
print -dpdf fig/199.pdf




 







