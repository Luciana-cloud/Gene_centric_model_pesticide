function [final,output] = template_calibration24D(p)

% This file is used for running the model calibration against 2,4-D data.
% The objective of this file is to calculate the SSE and plot together the
% distribution of the residuals. This file is only a template and it has to
% be adecuate for the particular model version. Here, we present the
% solution for the model version V0. This file can be used for manual and
% automatized calibration. 

%% DEFINE FIXED PARAMETERS %%

q(1) = 8;           % n_M - Conversion factor from mmol 2,4D to mmol C  [mmol C mmol 2,4D^-1]
q(2) = 0.1;         % n_S  - Switch function parameter
q(3) = 0.35;        % th_V - Average volumetric soil water content (cm^3 cm^-3)
q(4) = 1.1;         % rho_B - Bulk density of soils (g cm^-3)

%% TIME %%

base   = 9600*0.3+1;                    % Number of output times
t_base1 = linspace(0,100*0.3,base);     % Simulation period and vector of output times (d)
t_base1(end)=[];
t_base2 = linspace(30,35,2400*0.3+1);   % Simulation period and vector of output times (d)
t_basef = [t_base1 t_base2];
t       = [193; 577; 769; 865; 1021; 1137; 1225; 1341; 1411; 1515; 1599;...
           1687; 1889; 2177; 2371; 2371; 2387; 2411; 2464; 2466; 2515; 2563; ...
           2603; 2663; 2751; 2855; 3406];
       
t_base  = t_basef(t);                   % Respike at time 24.6875 days t_base(16)

tspan =  [t_base(1) t_base(2) t_base(3) t_base(4) t_base(5) t_base(6) t_base(7)...
          t_base(8) t_base(9) t_base(10) t_base(11) t_base(12) t_base(13) t_base(14) t_base(15)];   % Time for first degradation phase
tspan1 = [t_base(16) t_base(17) t_base(18) t_base(19) t_base(20) t_base(21) ...
          t_base(22) t_base(23) t_base(24) t_base(25) t_base(26) t_base(27)];                       % Time for second degradation phase
      
%% OPTIONS FOR THE SOLVER %%

    abstol = 1e-9;
    reltol = 1e-7;
    o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent); %,'Jacobian',@MCPA1_jac);%,'NonNegative',1:7);

%% INITIAL CONDITIONS FOR FIRST PESTICIDE APPLICATION %%
        
        K_F     = 10^p(15);
        n_F     = p(16);
        C_tot   = 7.4e-04;      % Total initial C conc. mmol g-1 soil
        C_L0    = 0;
        th_V    = q(3);         % Average volumetric soil water content (cm^3 cm^-3)
        rho_B   = q(4);         % Bulk density of soils (g cm^-3)
        f_1     = 10^p(5);      % f_CG - Conversion factor (mmol of C per # of gene)[mmol C gene(-1)]
        
        fun=@(C_L0) MCPA_init(C_L0,C_tot,K_F,n_F,th_V,rho_B);
        
        %  switch off solver progress information on solver progress of fsolve
        f_opts = optimoptions('fsolve','Display','none');
        
        % calculate initial MCPA in solution phase
        C_L0 = fsolve(fun,C_L0,f_opts); % C in solution
        
        % DNA %
        
        DNA_T     = 10^p(18)*f_1;       % total initial gene abundance [mmol C g-1]
        a         = p(17);              % Parameter that determines the percentage of relic DNA in soils
        
c(1) = 0;                       % Active DNA _tfdA  [mmolC g-1]
c(2) = DNA_T*(1-a);             % Inactive DNA _tfdA  [mmolC g-1]
c(3) = DNA_T*a;                 % Dead DNA in soils in solution [mmolC g-1] 
c(4) = C_L0;                    % MCPA in solution (mass per volume of water) [mmol C cm^-3]
c(5) = 0;                       % Total CO2 [mmolC g-1]

 try
     warning off
     tic
     [ty,cu] = ode15s(@MCPA1_c_g,tspan,c,o_opts,p',q);     

 catch ME
     warning off
     ty = t_base(1:15);
     cu = ones(length(t_base(1:15)),length(c))*1e+99;     
 end
if length(cu) < length(t_base(1:15))
    cu = ones(length(t_base(1:15)),5)*1e+99;
end

if isreal(cu)==0
    cu = ones(length(t_base(1:15)),5)*1e+99;    
end

%% INITIAL CONDITIONS FOR SECOND PESTICIDE APPLICATION %%

C_totA = 7.4e-04;    % Total initial 2,4 conc. mmol g-1 soil
C_L0  = 0;
        
fun=@(C_L0) MCPA_init(C_L0,C_totA,K_F,n_F,th_V,rho_B);
        
%  switch off solver progress information on solver progress of fsolve
f_opts = optimoptions('fsolve','Display','none');
        
% calculate initial MCPA in solution phase
C_L0 = fsolve(fun,C_L0,f_opts); % MCPA in solution

c(1) = cu(end,1);              % Active DNA _tfdA  [mmolC g-1]     
c(2) = cu(end,2);              % Inactive DNA _tfdA  [mmolC g-1]
c(3) = cu(end,3);              % Dead DNA _tfdA  [mmolC g-1] 
c(4) = cu(end,4)+C_L0;         % MCPA in solution (mass per volume of water) [mmol C cm^-3]
c(5) = 0;                      % Total CO2 [mmolCg-1] 

 try
     warning off
     tic
     [ty,cu1] = ode15s(@MCPA1_c_g,tspan1,c,o_opts,p',q);     

 catch ME
     warning off
     ty = t_base(16:end);
     cu1 = ones(length(t_base(16:end)),length(c))*1e+99;     
 end
 
if length(cu1) < length(t_base(16:end))
    cu1 = ones(length(t_base(16:end)),5)*1e+99;
end

if isreal(cu1)==0
    cu1 = ones(length(t_base(16:end)),5)*1e+99;    
end

%% COMPUTING FULL TIME SERIES OF THE MODEL OUTPUTS %%

    DNA_a = vertcat(cu(:,1),cu1(:,1));  % [mmol C g-1]
    DNA_i = vertcat(cu(:,2),cu1(:,2));  % [mmol C g-1]
    DNA_s = vertcat(cu(:,3),cu1(:,3));  % [mmol C g-1]  
    P     = vertcat(cu(:,4),cu1(:,4));  % [mmol C cm-3]
    CO2   = vertcat(cu(:,5),cu1(:,5));  % [mmol C g-1]

%% COMPUTING FULL TIME SERIES OF THE SIMULATED VARIABLES TO MATCH TO THE DATA %%

        Initial = 7.4e-04;          % Intial 2,4-D applied [mmol C g-1]
        n_H     = p(2);             % [-]
        K_G     = 10^p(3);          % [mmol C cm^-3]
        f_T     = 10^p(1);          % [transcripts/gene]
        RNA     = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1)); % [transcripts/gene]

    mineralA (:,1)       = CO2.*100/Initial;                % [%]
    transcriptsA (:,1)   = RNA.*(DNA_a*f_1^(-1));           % [transcripts g-1]
    DNA_totA (:,1)       = (DNA_a+DNA_s+DNA_i)*f_1^(-1);    % [gene g-1]
    
%% CALCULATION OF THE SSE %%

output = [mineralA(1:15,1);mineralA(17:end,1); (transcriptsA([7 9 10 12 13 15 18 20 21 22 23 25 27],1));...
        (DNA_totA([1 3 7 9 10 12 13 15 18 20 21 22 23 25 27],1))];      % Output to be compared to data       
    
mod_data2 = [output(1:26);output(27:39);output(40:54)];                 % Output to be compared to data 

% Calling Data %

% The data can be found in file Data_file, Pesticide_data.

load('data_24D.txt')

meas_data2 = [data_24D(1:26,2);(data_24D(27:54,2))];  % Mean values
stad_data2 = [data_24D(1:26,3);(data_24D(27:54,3))];  % Standard deviation

% Plotting residuals %

t1 = linspace(1,54,54);
resi = (((mod_data2-meas_data2)./stad_data2).^2);
plot(t1(1:26),resi(1:26),'-o','Color','r')
hold on
plot(t1(27:39),resi(27:39),'-o','Color','b')
plot(t1(40:54),resi(40:54),'-o','Color','k')
title('Residuals after manual calibration for model with dormancy')
xlabel('Individual values','FontSize',12);
hold off

% Calculating SSE %

res  =(sum(((mod_data2-meas_data2)./stad_data2).^2))/length(model2_15); % square of residuals
final = res^0.5

end