%% COMPARISON %%

% We calculated here the information criteria AICc and BIC. We present the
% example of the full model version V0 and the data has to be adapted for
% the other models. The data of the SSE results in the file Data_file,
% SSE_Outputs. 

%% AICc and BIC %%

load('SSE_V0.mat')

m           = 18;       % Number of model parameters 
SSEdor      = RES;      % SSE
n           = 54;       % Number of observations

AIC_0       =  2*m+n*log(SSEdor/n)+(2*m*(m+1)/(n-m-1));

BIC_0       = -2*log(SSEdor)+2*log(n)*m;