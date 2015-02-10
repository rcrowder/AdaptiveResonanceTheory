function [x_t2,z_t2,B_t2]=NDTransientCells(x_t1,z_t1,I_t1,B1,K2,Tb)
%
% NDTransientCells :: Function to simulate the change-sensitive receptors
% or the non-directional transient cells in the Retina
%
%% Inputs
%
% x_t1 :: Stage 1 neuronal activities at time instant t1
%
% z_t1 :: Habituative transmitter activities at time instant t1
% 
% Tb :: Output threshold of non-directional transient cells
%
% I_t1 :: Luminance increments at time instant t1 in the stimulus
%
%% Outputs
% 
% x_t2 :: Stage 1 neuronal activities at time instant t2
%
% z_t2 :: Habituative transmitter activities at time instant t2
%
% B_t2 :: Non-directional transient cell activities at time instant t2
%
%% Reference
% Grossberg, S., & Rudd, M.E. (1992). Cortical dynamics of visual motion perception: Short-range and long-range apparent motion. Psychological Review, 99, 78-121.
%
%% Author
% Praveen K. Pilly (advaitp@gmail.com)
%
%% License policy
% Written by Praveen K. Pilly, Department of Cognitive and Neural Systems, Boston University
% Copyright 2009, Trustees of Boston University
%
% Permission to use, copy, modify, distribute, and sell this software and its documentation for any purpose is hereby granted
% without fee, provided that the above copyright notice and this permission notice appear in all copies, derivative works and
% associated documentation, and that neither the name of Boston University nor that of the author(s) be used in advertising or
% publicity pertaining to the distribution or sale of the software without specific, prior written permission. Neither Boston
% University nor its agents make any representations about the suitability of this software for any purpose. It is provided "as
% is" without warranty of any kind, either express or implied. Neither Boston University nor the author indemnify any
% infringement of copyright, patent, trademark, or trade secret resulting from the use, modification, distribution or sale of
% this software.
%
%% Last modified
% Sept 8, 2009

% fixed time step of numerical integration (Forward Euler's method)
dt=0.001; % (in sec)

%% Parameters
A1=1; % Response speed scaling parameter for Stage 1 neurons
% B1=10; % Passive decay rate scaling parameter for Stage 1 neurons
Tb=0.1;
A2=1; % Response speed scaling parameter for habituative transmitters
% K2=50; % Habituation rate scaling parameter for habituative transmitters

x_t2=x_t1+dt*A1*(-B1*x_t1+(1-x_t1).*I_t1);
z_t2=z_t1+dt*A2*(1-z_t1-K2*x_t1.*z_t1);

b_t2=x_t2.*z_t2;
B_t2=max(b_t2-Tb,0); % thresholding and half-wave rectification

return
