function [mineralA,DNA_totA,transcriptsA,time,DNA_a,DNA_i,DNA_s,P,CO2] = template_simulations24D(p)

% This file is a template to run model simulations using the calibrated
% outputs. This is a template example using the version model V0 and it has
% to be adapted for the other model versions.

%% DEFINE FIXED PARAMETERS %%

q(1) = 8;       % n_M - Conversion factor from mmol 2,4D to mmol C  [mmol C mmol 2,4D^-1]
q(2) = 0.1;     % n_S  - Switch function parameter
q(3) = 0.35;    % th_V - Average volumetric soil water content (cm^3 cm^-3)
q(4) = 1.1;     % rho_B - Bulk density of soils (g cm^-3)

%% TIME %%

base   = 9600*0.3+1;                        % Number of output times
t_base1 = linspace(0,100*0.3,base);         % Simulation period and vector of output times (d)
t_base1(end)=[];
t_base2 = linspace(30,35,2400*0.3+1);       % Simulation period and vector of output times (d)
t_basef = [t_base1 t_base2];
t       = [193; 577; 769; 865; 1021; 1137; 1225; 1341; 1411; 1515; 1599;...
           1687; 1889; 2177; 2371; 2371; 2387; 2411; 2464; 2466; 2515; 2563; ...
           2603; 2663; 2751; 2855; 3406];
       
t_base  = t_basef(t);                       % Respike at time 24.6875 days t_base(15)

tspan =  [t_base(1) t_base(15)];            % Time for first degradation phase
tspan1 = [t_base(16) t_base(27)];           % Time for second degradation phase

%% OPTIONS FOR THE SOLVER %%

    abstol = 1e-9;
    reltol = 1e-7;
    o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent); %,'Jacobian',@MCPA1_jac);%,'NonNegative',1:7);

%% INITIAL CONDITIONS FOR FIRST PESTICIDE APPLICATION %%
        
        K_F     = 10^p(15);
        n_F     = p(16);
        C_tot   = 7.4e-04;          % Total initial C conc. mmol g-1 soil
        C_L0    = 0;
        th_V    = q(3);             % Average volumetric soil water content (cm^3 cm^-3)
        rho_B   = q(4);             % Bulk density of soils (g cm^-3)
        f_1     = 10^p(5);          % f_CG - Conversion factor (mmol of C per # of gene)[mmol C gene(-1)]
        
        fun=@(C_L0) MCPA_init(C_L0,C_tot,K_F,n_F,th_V,rho_B);
        
        %  switch off solver progress information on solver progress of fsolve
        f_opts = optimoptions('fsolve','Display','none');
        
        % calculate initial MCPA in solution phase
        C_L0 = fsolve(fun,C_L0,f_opts); % C in solution
        
        % DNA %
        
        DNA_T     = 10^p(18)*f_1;       % total initial gene abundance [mmol C g-1]
        a         = p(17);              % Parameter that determines the percentage of relic DNA in soils
        
c(1) = 0;                       % Active DNA _tfdA  [mmolC g-1]: DNA^ac
c(2) = DNA_T*(1-a);             % Inactive DNA _tfdA  [mmolC g-1]: DNA^ac
c(3) = DNA_T*a;                 % Stabilized DNA in soils in solution [mmolC g-1] 
c(4) = C_L0;                    % MCPA in solution (mass per volume of water) (mmol C cm^(-3))
c(5) = 0;                       % Total CO2 [mmolC g-1]

     tic
     [ty,cu] = ode15s(@MCPA1_c_g,tspan,c,o_opts,p',q);     

%% INITIAL CONDITIONS FOR SECOND PESTICIDE APPLICATION %%

C_totA = 7.4e-04;    % Total initial 2,4 conc. mmol g-1 soil
C_L0  = 0;
        
fun=@(C_L0) MCPA_init(C_L0,C_totA,K_F,n_F,th_V,rho_B);
        
%  switch off solver progress information on solver progress of fsolve
f_opts = optimoptions('fsolve','Display','none');
        
% calculate initial MCPA in solution phase
C_L0 = fsolve(fun,C_L0,f_opts); % MCPA in solution

c(1) = cu(end,1);              % Active DNA _tfdA  [gene copies cm-3]: DNA^ac     
c(2) = cu(end,2);              % Active DNA _tfdA  [gene copies cm-3]: DNA^ac
c(3) = cu(end,3);              % Stabilized DNA in soils in solution [gene copies cm-3] 
c(4) = cu(end,4)+C_L0;         % MCPA in solution (mass per volume of water) (mmol MCPA cm^(-3))
c(5) = 0;                      % Total CO2 [mmolCg-1] 

     [ty2,cu1] = ode15s(@MCPA1_c_g,tspan1,c,o_opts,p',q);     

%% COMPUTING FULL TIME SERIES OF THE MODEL OUTPUTS %%

    DNA_a = vertcat(cu(:,1),cu1(:,1));  % [mmol C g-1]
    DNA_i = vertcat(cu(:,2),cu1(:,2));  % [mmol C g-1]
    DNA_s = vertcat(cu(:,3),cu1(:,3));  % [mmol C g-1]  
    P     = vertcat(cu(:,4),cu1(:,4));  % [mmol C cm-3]
    CO2   = vertcat(cu(:,5),cu1(:,5));  % [mmol C g-1]
    time  = vertcat(ty,ty2);

%% COMPUTING FULL TIME SERIES OF THE SIMULATED VARIABLES TO MATCH TO THE DATA %%

        Initial = 7.4e-04;          % [mmol C g-1]
        n_H     = p(2);             % [-]
        K_G     = 10^p(3);          % [mmol C cm^-3]
        f_T     = 10^p(1);          % [transcripts/gene]
        RNA     = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1)); % [transcripts/gene]

    mineralA (:,1)       = CO2.*100/Initial;                % [%]
    transcriptsA (:,1)   = RNA.*(DNA_a*f_1^(-1));           % [transcripts g-1]
    DNA_totA (:,1)       = (DNA_a+DNA_s+DNA_i)*f_1^(-1);    % [gene g-1]
        
end
