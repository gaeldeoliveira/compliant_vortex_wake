%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   File    : Load Donqi duct airfoil, fuse it (2nd generation morphing)  %
%             and export airfoil to file (no computations)                %
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

%% Inputs
% Geometry
% Airfoil:  definition
airfoil_file = 'airfoils/donqio.air';   % Filename of reference airfoil
BPO_upper    = 16;                      % Number of parameters (order=degree+1) for description of airfoil upper (suction ) side shape
BPO_lower    = 16;                      % Number of parameters (order=degree+1) for description of airfoil lower (pressure) side shape
eta_airfoil  = -.2 ;                    % Fusing parameter 
flip_airfoil = false;                   % % If true airfoil is flipped upside down (boolean!) (post-morphing)

%% Add sources to matlab path
fs = filesep();                   % Folder separator is OS dependent
addpath([cd fs 'src_kirikou' ]);  % Add optimizer source code folder
addpath([cd fs 'src_optiflow']);  % Add optimizer source code folder

%%  Morph airfoil shape
% Scale, morph and flip airfoil coordinates
[px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = ...
            fuse_airfoil_coordinates(BPO_upper, BPO_lower, airfoil_file, eta_airfoil, flip_airfoil);
        
%%  Export results to airfoil file

% Bundle airfoil coordinates into single array
array_air_morphed = [px_air_morphed, py_air_morphed]; %#ok<NASGU>
array_air_scaled = [px_air_scaled  , py_air_scaled ] ;
                
% Build airfoil string file (for morphed)
if eta_airfoil >= 0
    airfoil_file_morphed = [airfoil_file(1:end-4) , '_eta_p', num2str(abs(eta_airfoil)) '.air'];
else
    airfoil_file_morphed = [airfoil_file(1:end-4) , '_eta_m', num2str(abs(eta_airfoil)) '.air'];
end
% Build airfoil string file (for scaled)
airfoil_file_scaled = [airfoil_file(1:end-4) , '_scaled.air'];
                
% Save to file in right format
save( airfoil_file_morphed , 'array_air_morphed','-ascii', '-double');
save( airfoil_file_scaled  , 'array_air_scaled','-ascii', '-double');

% Plot to make sure all is well
% Bundle airfoil coordinates into single array
plot(px_air_morphed, py_air_morphed); hold on;
plot(px_air_scaled  , py_air_scaled); grid on; axis equal;
legend(['Fused airfoil, eta=' num2str(eta_airfoil)], 'Scaled Donqi');


        

