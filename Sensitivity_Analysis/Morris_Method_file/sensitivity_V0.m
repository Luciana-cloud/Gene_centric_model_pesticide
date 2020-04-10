function [mineralA,transcriptsA,DNA_totA,DNA_a,DNA_i,DNA_s,RNA,P,CO2,rate,ty,SSE] = sensitivity_V0(p)

%% DEFINE FIXED PARAMETERS %%

q(1) = 8;       % n_M - Conversion factor from mmol 2,4D to mmol C  [mmol C mmol 2,4D^-1]
q(2) = 0.1;     % n_S  - Switch function parameter
q(3) = 0.35;    % th_V - Average volumetric soil water content (cm^3 cm^-3)
q(4) = 1.1;     % rho_B - Bulk density of soils (g cm^-3)

%% TIME %%

base   = 9600*0.3+1;                    % Number of output times
t_base1 = linspace(0,100*0.3,base);     % Nimulation period and vector of output times (d)
t_base1(end)=[];
t_base2 = linspace(30,35,2400*0.3+1);   % Nimulation period and vector of output times (d)
t_basef = [t_base1 t_base2];
t       = [193; 577; 769; 865; 1021; 1137; 1225; 1341; 1411; 1515; 1599;...
           1687; 1889; 2177; 2371; 2371; 2387; 2411; 2464; 2466; 2515; 2563; ...
           2603; 2663; 2751; 2855; 3406];
       
t_base  = t_basef(t);                   % Respike at time 24.6875 days t_base(16)

% For calibration:

tspan1   =  [t_base(1) t_base(2) t_base(3) t_base(4) t_base(5) t_base(6) t_base(7)...
          t_base(8) t_base(9) t_base(10) t_base(11) t_base(12) t_base(13) t_base(14) t_base(15)];

% For smoother simulations:

tspan2  =  [t_base(1) t_base(15)];


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
c(2) = (DNA_T*(1-a));             % Inactive DNA _tfdA  [mmolC g-1]: DNA^ac
c(3) = (DNA_T*a);                 % Stabilized DNA in soils in solution [mmolC g-1] 
c(4) = (C_L0);                    % MCPA in solution (mass per volume of water) (mmol C cm^(-3))
c(5) = 0;                       % Total CO2 [mmolC g-1]

 try
     warning off
     tic
     [ty,cu]    = ode15s(@version_V0,tspan2,c,o_opts,p',q);       %  Smooth simulations    
     [tyC,cuC]  = ode15s(@version_V0,tspan1,c,o_opts,p',q);       %  Calibration
     
 catch ME
     warning off
     tyC = t_base(1:15);
     cuC = ones(length(t_base(1:15)),length(c))*1e+99;     
 end
if length(cuC) < length(t_base(1:15))
    cuC = ones(length(t_base(1:15)),length(c))*1e+99;
end

%% COMPUTING FULL TIME SERIES OF THE MODEL OUTPUTS %%

        Initial = 7.4e-04;      % [mmol C g-1]
        n_H     = p(2);         % [-]
        K_G     = 10^p(3);      % [mmol C cm^-3]
        f_T     = 10^p(1);      % [transcripts/gene]
        mu_max  = 10^p(4);
        K_M     = 10^p(6);
        Y_P     = p(13);
        a_s     = p(12);
        aCO2    = p(14);
        
% Smooth simulations:
    
    DNA_a = ((cu(:,1)));        % [gene g-1]
    DNA_i = ((cu(:,2)));        % [gene C g-1]
    DNA_s = ((cu(:,3)));        % [gene C g-1]  
    P     = ((cu(:,4)));        % [mmol C cm-3]
    CO2   = ((cu(:,5)));        % [mmol C g-1]

         RNA     = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1)); % [transcripts/gene]
         mu_P    = mu_max*P.^(n_H+1).*(K_G^(n_H)+P.^(n_H)).^(-1).*(K_M+P).^(-1);

%% COMPUTING FULL TIME SERIES OF THE SIMULATED VARIABLES TO MATCH TO THE DATA %%

    mineralA (:,1)       = CO2.*100/Initial;
    transcriptsA (:,1)   = RNA.*(DNA_a)*f_1^(-1);          % [transcripts g-1]
    DNA_totA (:,1)       = (DNA_a+DNA_s+DNA_i)*f_1^(-1);   % [gene g-1]
    
    % Calibration:
    
    DNA_aC = ((cuC(:,1)));  % [gene g-1]
    DNA_iC = ((cuC(:,2)));  % [gene C g-1]
    DNA_sC = ((cuC(:,3)));  % [gene C g-1]  
    PC     = ((cuC(:,4)));  % [mmol C cm-3]
    CO2C   = ((cuC(:,5)));  % [mmol C g-1]
    
    RNAC     = f_T.*(PC.^(n_H)).*((K_G^(n_H)+PC.^(n_H)).^(-1)); % [transcripts/gene]
    rate = mu_P.*DNA_a*(1-Y_P)*Y_P^(-1)+DNA_s*a_s*aCO2;    
    
% compute full time series of measured parameters from model output

    mineralAC (:,1)       = CO2C.*100/Initial;
    transcriptsAC (:,1)   = RNAC.*(DNA_aC)*f_1^(-1);          % [transcripts g-1]
    DNA_totAC (:,1)       = (DNA_aC+DNA_sC+DNA_iC)*f_1^(-1);   % [gene g-1]
    
%% CALCULATION OF THE SSE %%

output    = [mineralAC(1:15,1); (transcriptsAC([7 9 10 12 13 15],1));...
        (DNA_totAC([1 3 7 9 10 12 13 15],1))];         
    
mod_data2 = [output(1:15);output(16:21);output(22:end)];
%          
load('data_24D.txt') % The data can be found in file Data_file, Pesticide_data.
% 
meas_data2 = [data_24D(1:15,2);(data_24D(27:32,2));(data_24D(40:47,2))];

stad_data2 = [data_24D(1:15,3);(data_24D(27:32,3));(data_24D(40:47,3))];
res  =(sum(((mod_data2-meas_data2)./stad_data2).^2))/length(meas_data2); % square of residuals
SSE = res^0.5;

end