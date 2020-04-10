% This script provides a basic application example                             
% of the Elementary Effects Test. Useful to get started with the EET.          
%                                                                              
% METHOD                                                                       
%                                                                              
% This script provides an example of application of the Elementary Effects     
% Test (EET) or 'method of Morris' (Morris, 1991; Saltelli et al., 2008).      
%                                                                              
% The EET is a One-At-the-Time method for global Sensitivity Analysis.         
% It computes two indices for each input:                                      
% i) the mean (mi) of the EEs, which measures the total effect of an input     
% over the output;                                                             
% ii) the standard deviation (sigma) of the EEs, which measures the degree     
% of interactions with the other inputs.                                       
% Both sensitivity indices are relative measures, i.e. their value does not    
% have any specific meaning per se but it can only be used in pair-wise        
% comparison (e.g. if input x(1) has higher mean EEs than input x(3) than      
% x(1) is more influential than x(3)).                                         
%                                                                              
% For an application example in the environmental domain, see for instance     
% Nguyen and de Kok (2007).                                                    
%                                                                              
% MODEL AND STUDY AREA                                                         
%                                                                              
% The model under study is the rainfall-runoff model Hymod                     
% (see help of function hymod_sim.m for more details)                          
% applied to the Leaf catchment in Mississipi, USA                             
% (Sorooshian et al., 1983).                                                   
% The inputs subject to SA are the 5 model parameters, and the scalar          
% output for SA is a metric of model performance.                              
%                                                                              
% INDEX                                                                        
%                                                                              
% Steps:                                                                       
% 1. Add paths to required directories                                         
% 2. Load data and set-up the HBV model                                        
% 3. Sample inputs space                                                       
% 4. Run the model against input samples                                       
% 5. Compute the elementary effects                                            
% 6. Example of how to repeat computions after adding up new                   
%    input/output samples.                                                     
%                                                                              
% REFERENCES                                                                   
%                                                                              
% Morris, M.D. (1991), Factorial sampling plans for preliminary                
% computational experiments, Technometrics, 33(2).                             
%                                                                              
% Nguyen, T.G. and de Kok, J.L. (2007). Systematic testing of an integrated    
% systems model for coastal zone management using sensitivity and              
% uncertainty analyses. Env. Mod. & Soft., 22, 1572-1587.                      
%                                                                              
% Saltelli, A., et al. (2008) Global Sensitivity Analysis, The Primer,         
% Wiley.                                                                       
%                                                                              
% Sorooshian, S., Gupta, V., Fulton, J. (1983). Evaluation of maximum          
% likelihood parameter estimation techniques for conceptual rainfall-runoff    
% models: Influence of calibration data variability and length on model        
% credibility. Water Resour. Res., 19, 251-259.                                
                                                                               
% This script prepared by Francesca Pianosi and Fanny Sarrazin                 
% University of Bristol, 2014                                                  
% mail to: francesca.pianosi@bristol.ac.uk                                     
                                                                               
%% Step 1 (add paths)    

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab    
% current directory to the SAFE directory. Otherwise, you may define           
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:            
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';                      
                                                                               
% Set current directory to 'my_dir' and add path to sub-folders:               
cd(my_dir)                                                              
addpath(genpath(my_dir))
                                                                               
%% Step 2 (setup the Hymod model)                                              

load('with_hybrid_15.mat')
                                                                               
% Number of uncertain parameters subject to SA:                                
M    = 18 ;                                                                     
% Parameter ranges (from literature):                                          
xmin = [-4 1 -10 -4 -14 -8 -10 -5 -7 -5 -5 -4 0.1 0.1 -2 0.8 0.1 3.5];                                                     
xmax = [3 10 6 5 -8 4 -2 2 -2 2 2 4 0.9 0.9 0 1 0.9 5];                                                     
  
DistrFun  = 'norm' ; % 
DistrPar  = { [x(1) 2.5]; [x(2) 1.5]; [x(3) 2.5]; [x(4) 1.5]; [x(5) 1.5]; [x(6) 2.5]; [x(7) 1.5]; [x(8) 1.75]; ...
              [x(9) 1.5]; [x(10) 2.5]; [x(11) 1.5]; [x(12) 2]; [x(13) 0.1];[x(14) 0.07];[x(15) 0.35];...
              [x(16) 0.025]; [x(17) 0.09]; [x(18) 0.5]} ; % Parameter ranges

% Name of parameters (will be used to costumize plots):                        
X_labels = {'f_T','n_H','K_G','mu_max','f_1','K_M','C_T','a_a','a_i',...
            'k_r','k_d','a_s','Y_P','aCO2','K_F_P','n_F_P','a_r','DNT_T_°'} ;                                    

out      = 13;  % Number of trajectory outputs (see "MCPA_dormax.m")

% Define output:                                                               
myfun = 'V0_dormax' ;                                                          
                                                                               
%% Step 3 (sample inputs space)                                                
                                                                               
r = 50 ; % Number of Elementary Effects                                       
% [notice that the final number of model evaluations will be equal to          
% r*(M+1)]                                                                     
                                                                               
% option 1: use the sampling method originally proposed by Morris (1991):      
% L = 6  ; % number of levels in the uniform grid                              
% design_type  = 'trajectory'; % (note used here but required later)           
% X = Morris_sampling(r,xmin,xmax,L); % (r*(M+1),M)                            
                                                                               
% option 2: Latin Hypercube sampling strategy                                  
SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';                                                        
% other options for design type:                                               
%design_type  = 'trajectory';   

%%
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);   

save('param.mat','X');           % Mean of the elementary effects

%% Step 4 (run the model)                                                      
Y = model_evaluation(myfun,X) ; %,rain,evap,flow) ; % size (r*(M+1),1)              

save('outputs.mat','Y');           % Mean of the elementary effects

%% Step 5 (Computation of the Elementary effects)
                                   
mi          = zeros(out,M); 
sigma       = zeros(out,M);
mi_sd       = zeros(out,M);
sigma_sd    = zeros(out,M);
mi_lb       = zeros(out,M);
sigma_lb    = zeros(out,M);
mi_ub       = zeros(out,M);
sigma_ub    = zeros(out,M);
EE          = zeros(out*r,M);
a           = [0 1 2 3 4 5 6 7 8 9 10 11 12]; % out - 1
rr          = [r/5:r/5:r] ;                                                           
m_r         = zeros(length(rr)*out,M);
s_r         = zeros(length(rr)*out,M);
m_lb_r      = zeros(length(rr)*out,M);
m_ub_r      = zeros(length(rr)*out,M);

for j = 1:out
    
% Compute Elementary Effects:                                                  
[ mi(j,:), sigma(j,:) ] = EET_indices(r,xmin,xmax,X,Y(:,j),design_type); 

end 

save('mean_EE.mat','mi');           % Mean of the elementary effects
save('standard_EE.mat','sigma');    % Standard deviation of the elementary effects

% for j = 1:9
% % Plot results in the plane (mean(EE),std(EE)):                                
% EET_plot(mi(j,:), sigma(j,:),X_labels )
% end

% Use bootstrapping to derive confidence bounds:        

for j = 1:out

Nboot=15;                                                                     
[mi(j,:),sigma(j,:),EE(((r*a(j))+1):((r*a(j))+r),1:M),mi_sd(j,:),sigma_sd(j,:),mi_lb(j,:),...
sigma_lb(j,:),mi_ub(j,:),sigma_ub(j,:)] = EET_indices(r,xmin,xmax,X,Y(:,j),design_type,Nboot);

% Plot bootstrapping results in the plane (mean(EE),std(EE)):                  
% EET_plot(mi(j,:),sigma(j,:),X_labels,mi_lb(j,:),mi_ub(j,:),sigma_lb(j,:),sigma_ub(j,:))
end

save('mean_EE.mat','mi');           % Mean of the elementary effects
save('standard_EE.mat','sigma');    % Standard deviation of the elementary effects
save('matrix_EE.mat','EE');         % Matrix of elementary effects
save('standard_mi.mat','mi_sd');    % Standard deviation of 'mi' across Nboot resamples
save('standard_si.mat','sigma_sd'); % Standard deviation of 'sigma' across Nboot resamples
save('lower_mi.mat','mi_lb');       % Lower bound of 'mi'
save('lower_si.mat','sigma_lb');    % Lower bound of 'sigma'
save('upper_mi.mat','mi_ub');       % Upper bound of 'mi'
save('upper_si.mat','sigma_ub');    % Upper bound of 'sigma'

% Repeat computations using a decreasing number of samples so as to assess     
% if convergence was reached within the available dataset:

for j = 1:out
m_r(((length(rr)*a(j))+1):((length(rr)*a(j))+length(rr)),1:M) = EET_convergence(EE(((r*a(j))+1):((r*a(j))+r),1:M),rr); 
end

% save('mean_EEs.mat','m_r');         % Mean EEs at different sampling size

for j = 1:out

    % Repeat convergence analysis using bootstrapping:                             
Nboot = 15;                                                                   
[m_r(((length(rr)*a(j))+1):((length(rr)*a(j))+length(rr)),1:M),s_r(((length(rr)*a(j))+1):((length(rr)*a(j))+length(rr)),1:M),...
 m_lb_r(((length(rr)*a(j))+1):((length(rr)*a(j))+length(rr)),1:M),m_ub_r(((length(rr)*a(j))+1):((length(rr)*a(j))+length(rr)),1:M)] = ...
 EET_convergence(EE(((r*a(j))+1):((r*a(j))+r),1:M),rr,Nboot); 

end  

save('mean_EEs.mat','m_r');         % Mean EEs at different sampling size
save('Standard_EEs.mat','s_r');     % Std of EEs at different sampling size
save('lower_EEs.mat','m_lb_r');     % Lower bound of mean EEs at different sampling size
save('upper_EEs.mat','m_ub_r');     % Upper bound of mean EEs
