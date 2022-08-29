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
%   Purpose : Interface development                                       %
%   Authors : Gael de Oliveira and Vinit Dighe                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %close all
 %clear all
 
%% Envisioned interfaces
% [px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = fuse_airfoil_coordinates(BPO_upper, BPO_lower, airfoil_file, eta_airfoil, flip_airfoil);
 
%% Inputs
BPO_upper    = 16;                      % Number of parameters (order=degree+1) for description of airfoil upper (suction ) side shape
BPO_lower    = 16;                      % Number of parameters (order=degree+1) for description of airfoil lower (pressure) side shape
airfoil_file = 'airfoils/donqio.air';   % Filename of reference airfoil
eta_airfoil  = 0 ;                      % Fusing parameter 
flip_airfoil = true;                    % % If true airfoil is flipped upside down (boolean!) (post-morphing)
 
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

%% Import reference (Donqi) airfoil shape 

% % Load airfoil and fit it to CST parametrization
[x, imported_coordinates, regenerated_coordinates, parameters]= get_parameters_from_file(SF , airfoil_file);

% % Select reference polynomial description for donqi airfoil
donqi_parameters_full  = parameters.full_non_linear.full;
donqi_parameters_upper = parameters.full_non_linear.upper;
donqi_parameters_lower = parameters.full_non_linear.lower;
% Identify trailing edge thickness
donqi_parameters_zte   = donqi_parameters_full(BPO_upper+BPO_lower+1);
% Donqi airfoil parameter rebuild (== donqi_parameters_full)
donqi_parameters_rfull = [donqi_parameters_upper, donqi_parameters_lower, donqi_parameters_zte];

%% Make symmetric and flat variants of reference (Donqi) airfoil shape

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
% fltqi_parameters_rfull == intqi_parameters_fun(-1) %#ok<NOPTS,EQEFF>
% donqi_parameters_rfull == intqi_parameters_fun( 0) %#ok<NOPTS,EQEFF>
% symqi_parameters_rfull == intqi_parameters_fun( 1) %#ok<NOPTS,EQEFF>

% % Generate fused airfoil
intqi_parameters_rfull = intqi_parameters_fun(eta_airfoil);

% % Make corresponding coordinates
[tx_intqi , tz_intqi] = SD.generate_coordinates(200, intqi_parameters_rfull);

%% Now flip coordinates, if applicable
% Make actual airfoil from raw coordinates

% For (reconstructed, scaled/rotated/shifted) reference (DonQi) airfoil 
if not(flip_airfoil)
    % Keep airfoil as is
    px_air_scaled = tx_donqi;
    py_air_scaled = tz_donqi;
else
    % Flip airfoil and maintain CC point ordering
    px_air_scaled = flip(  tx_donqi);
    py_air_scaled = flip( -tz_donqi);
end


% For fused (morphed, 2nd generation) (IntQi) airfoil 
if not(flip_airfoil)
    % Keep airfoil as is
    px_air_morphed = tx_intqi;
    py_air_morphed = tz_intqi;
else
    % Flip airfoil and maintain CC point ordering
    px_air_morphed = flip(  tx_intqi);
    py_air_morphed = flip( -tz_intqi);
end








