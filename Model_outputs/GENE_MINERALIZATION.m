%% RATE OF MINERALIZATION VS GENE EXPRESSION %%

% This file is a template for the figure 4 of the main manuscript. We are
% presenting some examples of the model versions V0, V3 and V4 and for
% other different runs, the file should be adapted. The data for the
% current graphs can reproduce from the file
% TEMPLATE_SIMULATIONS_INTERPOLATIONS, and the data can be found Data_file,
% DREAM_outputs.

%% MODEL VERSIONS V0 %%

% Calling Data %

% Time %
load('time_V0.mat')

% DNA active %
load('DNAa_V0.mat')

% DNA relic %
load('DNAs_V0.mat')

% Growth Rate %
load('mu_P_V0.mat')

% Transcripts %
load('transcripts_V0.mat')

% Parameters %
load('V0_24D.mat')

% Y_P %
Y_P_A  = [Da(:,13,1);Da(:,13,2);Da(:,13,3)];

% a_s %
a_s_A  = [Da(:,12,1);Da(:,12,2);Da(:,12,3)];

% aCO2 %
aCO2_A = [Da(:,14,1);Da(:,14,2);Da(:,14,3)];

%% Transcripts vs. DNA_a*mu_P*(1-Y_P)*Y_P^(-1) +DNA_s*a_s*aCO2 %

figure('Name','V0')

pos = 226076; % 226076 - circular; 576507 - not complete; linear - 57982, linear 2 - 92699
plot(transcripts_A{pos}(1:117),DNAa_A{pos}(1:117).*mu_P_A{pos}(1:117)*(1-Y_P_A(pos))*Y_P_A(pos)^(-1)+DNAs_A{pos}(1:117)*10^a_s_A(pos)*aCO2_A(pos),'b','LineWidth',2.75)  % First pesticide application
hold on
plot(transcripts_A{pos}(118:end),DNAa_A{pos}(118:end).*mu_P_A{pos}(118:end)*(1-Y_P_A(pos))*Y_P_A(pos)^(-1)+DNAs_A{pos}(118:end)*10^a_s_A(pos)*aCO2_A(pos),'r','LineWidth',2.75)  % Second pesticide application
legend('1° Application','2° Application','Location','northwest')
hold off
box on
set(gca,'FontSize',30)
ylim([0 5e-4])
xlim([0 6e4])

%% MODEL VERSIONS V3 %%

% Time %
load('time_V3.mat')

% Active %
load('DNAa_V3.mat')

% Growth Rate %
load('mu_P_V3.mat')

% Transcripts %
load('transcripts_V3.mat')

% Parameters %
load('V3_24D.mat')

% Y_P %
Y_P_D  = [D4a(:,7,1);D4a(:,7,2);D4a(:,7,3)];

% aCO_2 %
a_a_D  = [D4a(:,6,1);D4a(:,6,2);D4a(:,6,3)];

% a_s %
aCO2_D = [D4a(:,8,1);D4a(:,8,2);D4a(:,8,3)];

%% Transcripts vs. DNA_a*mu_P*(1-Y_P)*Y_P^(-1)+DNA_a*a_a*aCO2 %

figure('Name','V3')

pos = 87821; % linear - 151214/ circular - 87821
plot(transcripts_D{pos}(1:74),DNAa_D{pos}(1:74).*mu_P_D{pos}(1:74)*(1-Y_P_D(pos))*Y_P_D(pos)^(-1)+DNAa_D{pos}(1:74)*10^a_a_D(pos)*aCO2_D(pos),'b','LineWidth',2.75)         % First pesticide application
hold on
plot(transcripts_D{pos}(75:end),DNAa_D{pos}(75:end).*mu_P_D{pos}(75:end)*(1-Y_P_D(pos))*Y_P_D(pos)^(-1)+DNAa_D{pos}(75:end)*10^a_a_D(pos)*aCO2_D(pos),'r','LineWidth',2.75) % Second pesticide application
box on
set(gca,'FontSize',30)
ylim([0 5e-4])
xlim([0 6e4])

%% MODEL VERSIONS V4 %%

% Time %
load('time_V4.mat')

% Active %
load('DNAa_V4E.mat')

% Growth Rate %
load('mu_P_V4.mat')

% Transcripts %
load('transcripts_V4.mat')

% Parameters %
load('out_UNREL_24D.mat')

% Y_P %
Y_P_E  = [D(:,6,1);D(:,6,2);D(:,6,3)];

% aCO_2 %
aCO2_E = [D(:,7,1);D(:,7,2);D(:,7,3)];

% a_s %
a_a_E  = [D(:,5,1);D(:,5,2);D(:,5,3)];

%% Transcripts vs. DNA_a*mu_P*(1-Y_P)*Y_P^(-1)+DNA_a*a_a*aCO2 %

pos = 87821; % circular - 117161//
plot(transcripts_E{pos}(1:68),DNAa_E{pos}(1:68).*mu_P_E{pos}(1:68)*(1-Y_P_E(pos))*Y_P_E(pos)^(-1)+DNAa_E{pos}(1:68)*10^a_a_E(pos)*aCO2_E(pos),'b','LineWidth',2.75)          % First pesticide application
hold on
plot(transcripts_E{pos}(69:end),DNAa_E{pos}(69:end).*mu_P_E{pos}(69:end)*(1-Y_P_E(pos))*Y_P_E(pos)^(-1)+DNAa_E{pos}(69:end)*10^a_a_E(pos)*aCO2_E(pos),'r','LineWidth',2.75)  % Second pesticide application
box on
set(gca,'FontSize',30)
ylim([0 5e-4])
xlim([0 6e4])