%% LOCAL SENSITIVITY ANALYSIS %%

% This file is used for local sensitivity analysis and it is a template
% that has to be adapted for the other models. We are presenting here the
% example for the model version V0 and only with 24_D data. 

%% CALLING THE DATA %%

% The data can be found in file Data_file, Pesticide_data.

load('data_24D.txt')                          % File containing the data
meas_data               = data_24D(:,2);      % Data vector
days_meas               = data_24D(:,1);      % Time of each data point

%% CALLING PARAMETER VALUES %%

load('best_24D_V0.mat') ;     % Best values obtained from the hybrid optimization. Data found in Data_file, Best_Outputs_24D     
p = x;
N           = 18;               % Number of parameters
parameters  = [10^p(1);p(2);10^p(3);10^p(4);10^p(5);10^p(6);10^p(7);10^p(8);10^p(9);10^p(10);...
               10^p(11);10^p(12);p(13);p(14);10^p(15);p(16);p(17);10^p(18)]';   % Parameter values without log-transformations; therefore, the model file has to change to eliminate the exponential parameters.
parametersT = parameters;       % Renaming the values vector

%% CALLING THE MODEL %%

[output]  = template_calibration24D(parameters);        % Output of the model (vector file of outputs. Same lenght as data points)

%% SENSITIVITY ANALYSIS %%

pert_size   = 0.01;                                     % Size of the perturbation
S           = zeros(length(output),N);                  % Sensitivity matrix
absS        = zeros(length(output),N);                  % Log-likelihood function

for j=1:N
    pert_pars    = parameters;
    pert_pars(j) = parameters(j)*(1+pert_size);
    
    [outputper]= D24_DORyo_nolog(pert_pars);            % Output of the model perturbing each parameter 0.01
 
    for i=1:length(output)                              % Number of observations. In my case, I have 54 observation, so I will get 54 sensitivity values per parameter
                if outputper(i) > 0
            S   (i,j) = (outputper(i)-output(i))*(output(i))^(-1)*pert_size^(-1);
            absS(i,j) = (outputper(i)-output(i))*(pert_size*parameters(j))^(-1);
                else
                    S(i,j)=0;
                    absS(i,j)=0;
                end
    end
end

%% SENSITIVITY SCORE %% 

% Overall sensitivity measure. Parameter with the highest overall sensitivity measure is the most
% identifiable parameter.

SC          = length(S)*rms(S).^2';         % Overall sensitivity measure (OSM)
X_labels    = {'f_T','n_H','K_G','mu_max','f_1','K_M','C_T','a_a','a_i',...
               'k_r','k_d','a_s','Y_P','aCO2','K_F_P','n_F_P','a_r','DNT_T_'}; 
Final       = [X_labels;num2cell(SC')]';    % It contains the OSM and the name of the variables

%% IDENTIFIABILITY SCORE %%

Z                  = S; 
X                  = [];
[max0, max_index0] = max(length(Z)*rms(Z).^2');     % Parameter with the highest overall sensitivity 

IPOI               = zeros(N,1);                    % Indices_of_parameters_in_order_of_most_identifiable
ISPO               = zeros(N,1);                    % Ident_scores_of_parameters_in_order_of_most_identifiable
IPOI(1)            = max_index0;
ISPO(1)            = max0;

% Scores for the other paramters 

previous_max_index = max_index0;

for ii = 2:N
    X    = [X Z(:,previous_max_index)];
    Zhat = X*inv(X'*X)*X'*Z;
    R    = Z - Zhat;
    [maxP, max_indexP] = max(length(R)*rms(R).^2');
    IPOI(ii)            = max_indexP;
    ISPO(ii)            = maxP;
    previous_max_index = max_indexP;
end

for j=1:N    
    IS(IPOI(j))=ISPO(j);    
end

Final          = [X_labels;num2cell(parametersT); num2cell(SC') ;num2cell(IS)]'; % It contains the OSM and the name of the variables

%% PRACTICAL OR POSTERIOR IDENTIFIABILITY %%

meas_std    = model2_15(:,3);               % Vector containing the standard deviation of the data points
W           = diag(meas_std.^2,0);          % Covariance matrix
FIM         = absS'*inv(W)*absS;            % Fisher information matrix
FIM(11,:)    = [];
FIM(7,:)    = [];
FIM(:,11)    = [];
FIM(:,7)    = [];
CRLB        = diag(inv(FIM));               % Cramer-Rao bound

%confidence interval:

CI              = 1.96*sqrt(CRLB);          % Confidence interval
parameters(11)  = [];
parameters(7)   = [];
PE              = (CI./parameters')*100;    % Percent_error

x               = CI;
b               = [x(1:6)' 0 x(7:end)'];    % This modification is only for the version V0 because we got two parameters with 0 sensitivity that have to be discarted to get the other values
CI              = [b(1:10) 0 b(11:end)];
y               = PE;
b1              = [y(1:6)' 0 y(7:end)'];
PE              = [b1(1:10) 0 b1(11:end)];

Final(:,5)      = num2cell(CI);
Final(:,6)      = num2cell(PE);

%% 95% INTERVAL CONFIDENCE %%

M               = length(meas_data);                        % Number of observation
RES             =(sum(((output-meas_data)./meas_std).^2));  % Square of residuals for the best obtained parameter combination
F95             = finv(0.95, N, M-N);                       % Fishers distribution
C95             = N/(M-N)*RES*F95;                          % Confidence region

JTJ             = transpose(absS)*absS;                     % Fisher Information Matrix (J^T.T)
JTJ(11,:)       = [];
JTJ(7,:)        = [];
JTJ(:,11)       = [];
JTJ(:,7)        = [];

cov             = inv(JTJ);                                 % Covariance Matrix
JTJ_diag        = diag(JTJ);
JTJ_inv_diag    = diag(cov);
Delta_D         = C95./sqrt(JTJ_diag);                      % Dependent confidence interval
percent_Delta_D = 100*Delta_D./parameters';
Delta_I         = C95.*sqrt(diag(cov));                     % Independent confidence interval
percent_Delta_I = 100*Delta_I./parameters';

x1              = Delta_D;
b2              = [x1(1:6)' 0 x1(7:end)'];
Delta_D         = [b2(1:10) 0 b2(11:end)];
y1              = percent_Delta_D;
b3              = [y1(1:6)' 0 y1(7:end)'];
percent_Delta_D = [b3(1:10) 0 b3(11:end)];

x2              = Delta_I;
b4              = [x2(1:6)' 0 x2(7:end)'];
Delta_I         = [b4(1:10) 0 b4(11:end)];
y2              = percent_Delta_I;
b5              = [y2(1:6)' 0 y2(7:end)'];
percent_Delta_I = [b5(1:10) 0 b5(11:end)];
% 
Final(:,7)      = num2cell(Delta_D);
Final(:,8)      = num2cell(percent_Delta_D);

Final(:,9)      = num2cell(Delta_I);
Final(:,10)     = num2cell(percent_Delta_I);

%% CORRELATION MATRIX %%

for i=1:length(parameters)
    for j=1:length(parameters)
        cor(i,j)=cov(i,j)/(sqrt(cov(i,i)*cov(j,j)));
    end
end

X_labels(11) =[];
X_labels(7) =[];

X_labels1{1}     = [];
X_labels2 = [X_labels1 X_labels];

out     = [X_labels;num2cell(cor)];
corfin  = [X_labels2' out];

% Plot %

figure('Name','Correlation')

xlabels = {'fT','nH','KG','mu_max','f1','KM','aa','ai','kr',...
           'as','YP','aCO2','KFP','nFP','ar','DNA-T'};
ylabels = {'fT','nH','KG','mu_max','f1','KM','aa','ai','kr',...
           'as','YP','aCO2','KFP','nFP','ar','DNA-T'};
       
       map = [1 0 0
    0.94 0.41 0.17
    0 1 0
    1 1 0
    1 1 1
    1 1 0
    0 1 0
    0.94 0.41 0.17
    1 0 0];

g = colormap(map);       
h = heatmap(xlabels,ylabels,cor,'Colormap',g);
h.ColorLimits = [-1 1];

%% FINAL TABLE %%

headings = {'Parameter';'Best Fit Value';'Sensitivity';'Identifiability';...
            'Cramer-Rao boundary';'Percentage Error';'Dependent CI';'Percentage Error';...
            'Independent CI';'Percentage Error'}';
        
Summary = [headings;Final]; 
