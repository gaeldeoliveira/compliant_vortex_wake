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
%   Purpose : Interface test                                              %
%   Authors : Gael de Oliveira and Vinit Dighe                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all
 

%% Inputs
BPO_upper    = 16;                      % Number of parameters (order=degree+1) for description of airfoil upper (suction ) side shape
BPO_lower    = 16;                      % Number of parameters (order=degree+1) for description of airfoil lower (pressure) side shape
airfoil_file = 'airfoils/donqio.air';   % Filename of reference airfoil
eta_airfoil  = -.2 ;                    % Fusing parameter 
flip_airfoil = true;                    % % If true airfoil is flipped upside down (boolean!) (post-morphing)

%% Add sources to matlab path
fs = filesep();                   % Folder separator is OS dependent
addpath([cd fs 'src_kirikou' ]);  % Add optimizer source code folder
addpath([cd fs 'src_optiflow']);  % Add optimizer source code folder

%% Test interfaces
[px_air_morphed, py_air_morphed, px_air_scaled, py_air_scaled] = fuse_airfoil_coordinates(BPO_upper, BPO_lower, airfoil_file, eta_airfoil, flip_airfoil)

%% Plot results
figure(201); hold on; grid on; axis equal;
plot(px_air_scaled , py_air_scaled );
plot(px_air_morphed, py_air_morphed, '--');
legend('DonQi Airfoil', ['\eta = ' , num2str(eta_airfoil)]);
title('Test of fusing (2nd generation morphing) function');
print -dpdf fig/201.pdf




