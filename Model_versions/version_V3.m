function vf_ = version_V3(t,x,p4,q)

% The inactive bacteria pool and the dead bacteria pool were set to zero 
% The parameterKMwas also set to zero (KM= 0) to neglect pesticide-dependent
% growth and keep pesticide-dependent gene expression. 

%% STOCKS

DNA_a      = x(1);      % Active DNA [mmol C g-1]
P          = x(2);      % MPCA in solution [mmol C cm-3]

%% PARAMETERS %%

% BACTERIAL PARAMETERS

f_T    = 10^p4(1);      % Conversion factor [transcripts/gene]
n_H    = p4(2);         % Hill coeficient  [1]
K_G    = 10^p4(3);      % Activation coefficient of MCPA for gene expression [mmol C cm-3]
mu_max = 10^p4(4);      % Maximal specific growth rate of bacterial pesticide degraders bacteria A [d-1]
f_1    = 10^p4(5);      % f_CG - Conversion factor (mmol of C per # of gene)[mmol C gene(-1)]
a_a    = 10^p4(6);      % Decay rate of active?DNA?_tfdA [d-1]
Y_P    = p4(7);         % Substrate (MCPA)uptake efficiency of by bacterial pesticide degraders [1]
aCO2   = p4(8);         % Fraction of the degraded relic DNA that contributes to CO2 pool [1]

% SORPTION PARAMETERS

K_F_P  = 10^p4(9);      % Freundlich coeff of MCPA sorption isotherm (mmol MCPA g^-1 soil/(mmol MCPA cm^-3)^nF_MCPA)
n_F_P  = p4(10);        % Freundlich exponent of MCPA sorption isotherm (1)

%% VALUES OF CONSTANT %%

n_M   = q(1);           % n_M - Conversion factor from mmol MCPA to mmol C  [mmol C mmol MCPA^-1]
th_V  = q(3);           % th_V - Average volumetric soil water content (cm^3 cm^-3)
rho_B = q(4);           % rho_B - Bulk density of soils (g cm^-3)

%% BIOKINETIC FUNCTIONS %%

% Substrate (MCPA)  dependent specific growth rate of bacterial pesticide
% degraders bacteria  [d-1]
mu_P = mu_max*P^(n_H)*(K_G^(n_H)+P^(n_H))^(-1);

%% ---- vector fields = right hand sides of ODE system ---- %%

vf_ = zeros(3,1);

% Active DNA_tfdA  - [mmol C g-1]
vf_(1) = DNA_a*(mu_P-a_a);

% Solution phase MCPA concentration - mmmol C cm^(-3) soil solution)
vf_(2) = -(mu_P*DNA_a*(th_V^(-1)*rho_B)*(Y_P^(-1)))*...
          (1+th_V^(-1)*rho_B*P^(n_F_P-1)*K_F_P*n_F_P)^(-1);

% Total CO2 [mmol C g-1]
vf_(3)= DNA_a*mu_P*(1-Y_P)*Y_P^(-1)+DNA_a*a_a*aCO2; % +DNA_s*a_s*aCO2;

end