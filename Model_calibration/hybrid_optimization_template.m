%% HYBRIC OPTIMIZATION METHOD %%

% This file is a template file to run the hybrid optimization method. The
% presented example is only for the model version V0 and it has to be
% adapted for the other models. The files that have to be run here are the 
% template_calibration24D.m files adapted for the specific model version. 
% The first guess needs to be adapted as well for the different versions. 
% This file can be used for calibration against MCPA and 2,4-D data. 


ObjectiveFunction = @D24_DORyo; % Calling the Function

load('best_dormancy.mat')           % First guess value of the parameters for version V0
p0 = p;

lb = [-4 1 -10 -4 -14 -8 -10 -5 -7 -5 -5 -4 0.1 0.1 -2 0.8 0.1 3.5];        % Min values
ub = [3 10 6 5 -8 4 -2 2 -2 2 2 4 0.9 0.9 0 1 0.9 5];                       % Max values

hybridopts  = optimoptions('fmincon','OptimalityTolerance',1e-10,'MaxFunEvals',100000);
options     = optimoptions('simulannealbnd','HybridFcn',{'fmincon',hybridopts});%,'MaxIter',100);
[x,fval]    = simulannealbnd(ObjectiveFunction,p0,lb,ub,options);
