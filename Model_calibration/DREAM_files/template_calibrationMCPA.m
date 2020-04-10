function [output] = template_calibrationMCPA(p)

%% DEFINE FIXED PARAMETERS %%

q(1) = 9;               % n_M - Conversion factor from mmol MCPA to mmol C  [mmol C mmol MCPA^-1]
q(2) = 0.1;             % n_S  - Switch function parameter
q(3) = 0.35;            % th_V - Average volumetric soil water content (cm^3 cm^-3)
q(4) = 1.1;             % rho_B - Bulk density of soils (g cm^-3)

%% TIME %%

base    = 4800*0.4+1;                   % Number of output times
t_base  = linspace(0,100*0.8,base);     % Simulation period and vector of output times (d)
t       = [1;25;97; 169; 265; 346; 419; 529; 649; 792; 793; 794; 796;802;...
           816;829;865;937;1033;1105;1273];
tspan   = [t_base(t(1)) t_base(t(2)) t_base(t(3)) t_base(t(4)) t_base(t(5)) ...
           t_base(t(6)) t_base(t(7)) t_base(t(8)) t_base(t(9)) t_base(t(10))];      % Time for first degradation phase
tspan1  = [t_base(t(11)) t_base(t(12)) t_base(t(13)) t_base(t(14)) t_base(t(15)) ...
           t_base(t(16)) t_base(t(17)) t_base(t(18)) t_base(t(19)) t_base(t(20)) ...
           t_base(t(21)) 67];                                                       % Time for first degradation phase

%% OPTIONS FOR THE SOLVER %%

abstol = 1e-9;
reltol = 1e-7;
o_opts = odeset('AbsTol',abstol,'RelTol',reltol,'Events',@stopevent,'NonNegative',1:5); % 'Events',@stopevent);%,'NonNegative',1:7);
  
%% INITIAL CONDITIONS FOR FIRST PESTICIDE APPLICATION %%
        
K_F   = 10^p(15);
n_F   = p(16);
C_tot = 9e-04;    % Total initial MCPA conc. mmol g-1 soil
C_L0  = 0;
th_V  = q(3);     % Average volumetric soil water content (cm^3 cm^-3)
rho_B = q(4);     % Bulk density of soils (g cm^-3)
f_1   = 10^p(5);   % f_CG - Conversion factor (mmol of C per # of gene)[mmol C gene(-1)]

fun=@(C_L0) MCPA_init(C_L0,C_tot,K_F,n_F,th_V,rho_B);
        
%  switch off solver progress information on solver progress of fsolve
f_opts = optimoptions('fsolve','Display','none');
        
% calculate initial MCPA in solution phase
C_L0 = fsolve(fun,C_L0,f_opts); % MCPA in solution
        
% DNA %
                
DNA_T     = 10^p(18)*f_1;    % total initíal gene abundance [gene copies cm-3]
a         = p(17);    % Parameter that determines the percentage of relic DNA in soils
     
c(1) = 0;                       % Active DNA _tfdA  [mmolC g-1]: DNA^ac
c(2) = DNA_T*(1-a);             % Inactive DNA _tfdA  [mmolC g-1]: DNA^ac
c(3) = DNA_T*a;                 % Stabilized DNA in soils in solution [mmolC g-1] 
c(4) = C_L0;                  % MCPA in solution (mass per volume of water) (mmol C cm^(-3))
c(5) = 0;                       % Total CO2 [mmolC g-1]   
       
 try
     warning off
     tic
     [ty,cu] = ode15s(@version_V0,tspan,c,o_opts,p',q);     

 catch ME
     warning off
     ty = length(tspan);
     cu = ones(length(tspan),length(c))*1e+99;     
 end
if length(cu) < length(tspan)
    cu = ones(length(tspan),length(c))*1e+99;
end

if isreal(cu)==0
    cu = ones(length(tspan),length(c))*1e+99;    
end

%% INITIAL CONDITIONS FOR SECOND PESTICIDE APPLICATION %%

C_totA  = 9e-04;    % Total initial 2,4 conc. mmol g-1 soil
C_L0    = 0;
        
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

 try
     warning off
     tic
     [ty1,cu1] = ode15s(@version_V0,tspan1,c,o_opts,p',q);     

 catch ME
     warning off
     ty1 = length(tspan1);
     cu1 = ones(length(tspan1),length(c))*1e+99;     
 end
if length(cu1) < length(tspan1)
    cu1 = ones(length(tspan1),length(c))*1e+99;
end

if isreal(cu1)==0
    cu1 = ones(length(tspan1),length(c))*1e+99;    
end

%% COMPUTING FULL TIME SERIES OF THE MODEL OUTPUTS %%

    DNA_a = vertcat(cu(:,1),cu1(:,1));  % [mmol C g-1]
    DNA_i = vertcat(cu(:,2),cu1(:,2));  % [mmol C g-1]
    DNA_s = vertcat(cu(:,3),cu1(:,3));  % [mmol C g-1]  
    P     = vertcat(cu(:,4),cu1(:,4));  % [mmol C cm-3]
    CO2   = vertcat(cu(:,5),cu1(:,5));  % [mmol C g-1]
    time  = vertcat(ty,ty1);
    
%% COMPUTING FULL TIME SERIES OF THE SIMULATED VARIABLES TO MATCH TO THE DATA %%

        Initial = 9e-04;            % [mmol C g-1]
        n_H     = p(2);             % [-]
        K_G     = 10^p(3);          % [mmol MCPA cm^-3]
        f_T     = 10^p(1);          % [transcripts/gene]
        RNA     = f_T.*(P.^(n_H)).*((K_G^(n_H)+P.^(n_H)).^(-1));    % [transcripts/gene]
        
        mineralA (:,1)       = CO2.*100/Initial;
        transcriptsA (:,1)   = RNA.*(DNA_a*f_1^(-1));               % [transcripts g-1]
        DNA_totA (:,1)       = (DNA_a+DNA_s+DNA_i)*f_1^(-1);        % [gene g-1]
    
%% OPTIMIZATION %%

output  = [mineralA([2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 22]);...
         (transcriptsA([4 7 9 10 12 13 14 15 16 17 18 19 20 21]));...
         (DNA_totA([1 2 4 7 9 10 12 13 14 15 16 17 18 19 20 21]))];

end