% Template file for DREAM for model V0 against MCPA. This template has to
% be adapted for the different model versions. 

%% Problem settings defined by user

DREAMPar.d = 18;                         % Dimension of the problem
DREAMPar.N = 3;                          % Number of Markov chains
DREAMPar.T = 100000;                     % Number of generations
DREAMPar.lik = 12;                       % Model output is simulation: Gaussian likelihood function

%% Provide information parameter space and initial sampling

load('24D_ModelA.mat')                                  % Prior values of the model parameters

Par_info.initial        = 'normal';                     % Latin hypercube sampling
Par_info.boundhandling  = 'reflect';                    % Explicit boundary handling
Par_info.mu             = x;                            % if 'normal', define mean of distribution
Par_info.cov            = 2.5 * eye(DREAMPar.d);        % if 'normal', define covariance matrix
Par_info.min            = [-4 1 -10 -4 -14 -8 -10 -5 -7 -5 -5 -4 0.1 0.1 -2 0.8 0.1 3.5];                                                     
Par_info.max            = [3 10 6 5 -8 4 -2 2 -2 2 2 4 0.9 0.9 0 1 0.9 5]; 

%% Define name of function (.m file) for posterior exploration

Func_name = 'template_calibrationMCPA';

%% Define the measured streamflow data

% The data can be found in file Data_file, Pesticide_data.

load('data_MCPA.txt');                     % Data

Meas_info.Y         = (data_MCPA(:,2));    % Mean values
Meas_info.Sigma     = (data_MCPA(:,3));    % Standard deviation

%% Optional settings

options.save = 'yes';                       % Save workspace DREAM during run
options.parallel = 'yes';                   % Run in parallel
options.modout = 'yes';                     % Return model (function) simulations of samples (yes/no)?
% options.restart = 'yes';                  % Return model (function) simulations of samples (yes/no)?

%% Run the DREAM algorithm

[chain,output,fx,Z] = DREAM_ZS(Func_name,DREAMPar,Par_info,Meas_info,options);