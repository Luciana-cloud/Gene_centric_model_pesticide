function vf_ = version_V2(t,x,p1,q)

% The dead bacteria pool was set to zero (CdB= 0). To account for cellular 
% decay, a fraction153of the decaying active and inactive bacteria was set 
% to directly contribute to the CO2 pool. 

%% STOCKS

DNA_a      = x(1); % Active DNA [mmol C g-1]
DNA_i      = x(2); % Inactive DNA [mmol C g-1]
P          = x(3); % MPCA in solution [mmol C cm-3]

%% PARAMETERS %%

% BACTERIAL PARAMETERS

f_T    = 10^p1(1);      % Conversion factor [transcripts/gene]
n_H    = p1(2);         % Hill coeficient  [1]
K_G    = 10^p1(3);      % Activation coefficient of MCPA for gene expression [mmol C cm-3]
mu_max = 10^p1(4);      % Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
f_1    = 10^p1(5);      % f_CG - Conversion factor (mmol of C per # of gene)[mmol C gene(-1)]
K_M    = 10^p1(6);      % Michaelis menten constant (mmol C cm-3)
C_T    = 10^p1(7);      % Threshold MCPA concentration for microbial growth (mmol C cm^(-3))
a_a    = 10^p1(8);      % Decay rate of active?DNA?_tfdA [d-1]
a_i    = 10^p1(9);      % Decay rate of inactive?DNA?_tfdA [d-1]
k_r    = 10^p1(10);     % Coefficient rate of reactivation of inactive DNA_tfdA [d-1]
k_d    = 10^p1(11);     % Coefficient rate of deactivation of inactive DNA_tfdA [d-1]
Y_P    = p1(12);        % Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
aCO2   = p1(13);        % Fraction of the degraded relic DNA that contributes to CO2 pool [1]

% SORPTION PARAMETERS

K_F_P  = 10^p1(14);     % Freundlich coeff of MCPA sorption isotherm (mmol MCPA g^-1 soil/(mmol MCPA cm^-3)^nF_MCPA)
n_F_P  = p1(15);        % Freundlich exponent of MCPA sorption isotherm (1)

%% VALUES OF CONSTANT %%

n_M   = q(1);           % n_M - Conversion factor from mmol MCPA to mmol C  [mmol C mmol MCPA^-1]
n_S   = q(2);           % n_S  - Switch function parameter
th_V  = q(3);           % th_V - Average volumetric soil water content (cm^3 cm^-3)
rho_B = q(4);           % rho_B - Bulk density of soils (g cm^-3)

%% BIOKINETIC FUNCTIONS %%

% Substrate (MCPA)  dependent specific growth rate of bacterial pesticide
% degraders bacteria  [d-1]
mu_P = mu_max*P^(n_H+1)*(K_G^(n_H)+P^(n_H))^(-1)*(K_M+P)^(-1);

% switch-function to express genes [1]
tau_P = (1+exp(C_T^(-1)*n_S^(-1)*(C_T-P)))^(-1);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(4,1);

% Active DNA_tfdA  - [mmol C g-1]
vf_(1) = DNA_i*tau_P*k_r-DNA_a*(1-tau_P)*k_d+DNA_a*(mu_P-a_a);

% Inactive DNA_tfdA - [mmol C g-1]
vf_(2) = DNA_a*(1-tau_P)*k_d-DNA_i*tau_P*k_r-DNA_i*a_i;

% Solution phase MCPA concentration - mmmol C cm^(-3) soil solution)
vf_(3) = -(mu_P*DNA_a*(th_V^(-1)*rho_B)*(Y_P^(-1)))*...
          (1+th_V^(-1)*rho_B*P^(n_F_P-1)*K_F_P*n_F_P)^(-1);

% Total CO2 [mmol C g-1]
vf_(4)= DNA_a*mu_P*(1-Y_P)*Y_P^(-1)+(DNA_i*a_i+DNA_a*a_a)*aCO2;

end