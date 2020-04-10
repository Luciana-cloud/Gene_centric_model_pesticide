%% MODEL DYNAMCS %%

% This file is intented to reproduce the figures S2 to S5 from the
% Supporting Information. The outputs for this file can be found in the
% file Data_file, Outputs_D_24_D and Outputs_D_MCPA.

%% FOR 2,4-D DYNAMICS %%

% CALLING DATA %

% Measured Time %

load('time_fin_2,4.mat')

%% Model V0 %
% Active DNA %
load('YacV0_5d.mat')
load('YacV0_50d.mat')
load('YacV0_95d.mat')
% Inactive DNA %
load('YiV0_5d.mat')
load('YiV0_50d.mat')
load('YiV0_95d.mat')
% Relic DNA %
load('YsV0_5d.mat')
load('YsV0_50d.mat')
load('YsV0_95d.mat')
% Total Pesticide %
load('YPTV0_5d.mat')
load('YPTV0_50d.mat')
load('YPTV0_95d.mat')

%% Model V3 %
% DNA %
load('YacV3_5d.mat')
load('YacV3_50d.mat')
load('YacV3_95d.mat')
% PT %
load('YPV3_50d.mat')
load('YPV3_95d.mat')
load('YPV3_5d.mat')

%% Model V4 %
% DNA %
load('YaV4_5d.mat')
load('YaV4_50d.mat')
load('YaV4_95d.mat')
% PT %
load('YPTV4_5d.mat')
load('YPTV4_50d.mat')
load('YPTV4_95d.mat')

%% Model V4P %
% DNA %
load('YaV4P_5d.mat')
load('YaV4P_50d.mat')
load('YaV4P_95d.mat')
% PT %
load('YPTV4P_5d.mat')
load('YPTV4P_50d.mat')
load('YPTV4P_95d.mat')

%% Model V1 %
% DNA %
load('YaV1_5d.mat')
load('YaV1_50d.mat')
load('YaV1_95d.mat')
load('YsV1_5d.mat')
load('YsV1_50d.mat')
load('YsV1_95d.mat')
% PT %
load('PT_V1_5d.mat')
load('PT_V1_50d.mat')
load('PT_V1_95d.mat')

%% Model V2 %
% DNA %
load('YaV2_5d.mat')
load('YaV2_50d.mat')
load('YaV2_95d.mat')
load('YiV2_5d.mat')
load('YiV2_50d.mat')
load('YiV2_95d.mat')
% PT %
load('PT_V2_5d.mat')
load('PT_V2_50d.mat')
load('PT_V2_95d.mat')

%% PLOTTING %%

%% Model V0 %%

% DNA pools %

figure('Name','Validation')
plot(time_fin,YacV0_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YacV0_95d,'color','none','LineWidth',1.25)
plot(time_fin,YiV0_5d,'color','none','LineWidth',1.25)
plot(time_fin,YiV0_95d,'color','none','LineWidth',1.25)
plot(time_fin,YsV0_5d,'color','none','LineWidth',1.25)
plot(time_fin,YsV0_95d,'color','none','LineWidth',1.25)
A = fill_betweenN(time_fin,YacV0_5d,YacV0_95d); % Active DNA
B = fill_betweenN1(time_fin,YiV0_5d,YiV0_95d);% Inactive DNA
C = fill_betweenN3(time_fin,YsV0_5d,YsV0_95d);  % Relic DNA
plot(time_fin,YacV0_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75)
plot(time_fin,YiV0_50d,'color',[90/256 180/256 172/256],'LineWidth',2.75)
plot(time_fin,YsV0_50d,'color',[173/256 221/256 142/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20)
hold off
box on 
legend([A B C],'Active','Inactive','Dead','Location','northwest')
      
% Total Pesticide in Soil %

figure('Name','Validation')
plot(time_fin,YPTV0_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV0_95d,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,YPTV0_5d,YPTV0_95d); % Total 
plot(time_fin,YPTV0_50d,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20)
hold off
box on 

%% Model V3 %%

% Total DNA %   

plot(time_fin,YacV3_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YacV3_95d,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YacD_5d,YacD_95d);
plot(time_fin,YacV3_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20)
hold off
box on 

% Total Pesticide in Soil %

figure('Name','Validation')
plot(time_fin,YPTV3_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV3_95d,'color','none','LineWidth',1.25)
fill_betweenM(time_fin,YPTV3_5d,YPTV3_95d); % Total 
plot(time_fin,YPTV3_50d,'color','m','LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20)
hold off
box on 

%% Model V4 %%

% DNA pools %  

plot(time_fin,YaV4_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV4_95d,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV4_5d,YaV4_95d);
plot(time_fin,YaV4_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20,'xTickLabel',[])

% Total pesticide in soil %

figure('Name','Validation')
hold on
plot(time_fin,YPTV4_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV4_95d,'color','none','LineWidth',1.25)
fill_betweenM(time_fin,YPTV4_5d,YPTV4_95d); % Total 
plot(time_fin,YPTV4_50d,'color','m','LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20,'xTickLabel',[],'yTickLabel',[])

%% Model V4p %%

% DNA pools %

plot(time_fin,YaV4p_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV4p_95d,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV4p_5d,YaV4p_95d);
plot(time_fin,YaV4p_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20)

% Total pesticide in soils %

figure('Name','Validation')
hold on
plot(time_fin,YPTV4p_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV4p_95d,'color','none','LineWidth',1.25)
fill_betweenM(time_fin,YPTV4p_5d,YPTV4p_95d); % Total 
plot(time_fin,YPTV4p_50d,'color','m','LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20,'xTickLabel',[],'yTickLabel',[])

%% Model V1 %%

% DNA % 

figure('Name','Model B-DNA')
plot(time_fin,YaV1_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV1_95d,'color','none','LineWidth',1.25)
plot(time_fin,YsV1_5d,'color','none','LineWidth',1.25)
plot(time_fin,YsV1_95d,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV1_5d,YaV1_95d); % blue 
fill_betweenN3(time_fin,YsV1_5d,YsV1_95d); % residual
A=plot(time_fin,YaV1_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75);
B=plot(time_fin,YsV1_50d,'color',[173/256 221/256 142/256],'LineWidth',2.75);
ylim([0 1.5e-3])
hold off
box on
set(gca,'FontSize',20,'xTickLabel',[]) 
legend ([A B],'Active','Relic')

% Total pesticide in soil %

figure('Name','Model B-P')
plot(time_fin,PT_V1_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,PT_V1_95d,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,PT_V1_5d,PT_V1_95d);
plot(time_fin,PT_V1_50d,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20)

%% Model V2 %%

% Total DNA %

figure('Name','Model C-DNA')
plot(time_fin,YaV2_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV2_95d,'color','none','LineWidth',1.25)
plot(time_fin,YiV2_5d,'color','none','LineWidth',1.25)
plot(time_fin,YiV2_95d,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV2_5d,YaV2_95d); % blue 
fill_betweenN1(time_fin,YiV2_5d,YiV2_95d); % red 
A = plot(time_fin,YaV2_50d,'color',[216/256 179/256 101/256],'LineWidth',2.75);
B = plot(time_fin,YiV2_50d,'color',[90/256 180/256 172/256],'LineWidth',2.75);
ylim([0 1.5e-3])
hold off
box on
set(gca,'FontSize',20,'xTickLabel',[])

% Total pesticide in soil %

figure('Name','Model C-P')
plot(time_fin,PT_V2_5d,'color','none','LineWidth',1.25)
hold on
plot(time_fin,PT_V2_95d,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,PT_V2_5d,PT_V2_95d);
plot(time_fin,PT_V2_50d,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20,'yTickLabel',[])

%% FOR MCPA DYNAMICS %%

% CALLING DATA %

% Measured Time %

load('time_fin.mat')

%% Model V0 %%
% DNA %
load('YacV0_5.mat')
load('YacV0_50.mat')
load('YacV0_95.mat')
load('YinV0_5.mat')
load('YinV0_50.mat')
load('YinV0_95.mat')
load('YsV0_5.mat')
load('YsV0_50.mat')
load('YsV0_95.mat')

% PT %
load('YPTV0_5.mat')
load('YPTV0_50.mat')
load('YPTV0_95.mat')

%% Model V3 %

% DNA %
load('YaV3_5.mat')
load('YaV3_50.mat')
load('YaV3_95.mat')

% PT %
load('YPTV3_5.mat')
load('YPTV3_50.mat')
load('YPTV3_95.mat')

%% Model V4 %%

% DNA %
load('YaV4_5.mat')
load('YaV4_50.mat')
load('YaV4_95.mat')

% PT %
load('YPTV4_5.mat')
load('YPTV4_50.mat')
load('YPTV4_95.mat')

%% Model V4p %%

% DNA %
load('YaV4p_5.mat')
load('YaV4p_50.mat')
load('YaV4p_95.mat')

% PT %
load('YPTV4p_5.mat')
load('YPTV4p_50.mat')
load('YPTV4p_95.mat')

%% PLOTTING %%

%% Model V0 %%

% Total DNA %

figure('Name','Validation')
plot(time_fin,YacV0_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YacV0_95,'color','none','LineWidth',1.25)
plot(time_fin,YinV0_5,'color','none','LineWidth',1.25)
plot(time_fin,YinV0_95,'color','none','LineWidth',1.25)
plot(time_fin,YsV0_5,'color','none','LineWidth',1.25)
plot(time_fin,YsV0_95,'color','none','LineWidth',1.25)
A = fill_betweenN(time_fin,YacV0_5,YacV0_95); % Active DNA
B = fill_betweenN1(time_fin,YinV0_5,YinV0_95);% Inactive DNA
C = fill_betweenN3(time_fin,YsV0_5,YsV0_95);  % Relic DNA
plot(time_fin,YacV0_50,'color',[216/256 179/256 101/256],'LineWidth',2.75)
plot(time_fin,YinV0_50,'color',[90/256 180/256 172/256],'LineWidth',2.75)
plot(time_fin,YsV0_50,'color',[173/256 221/256 142/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
box on 
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])
       
% Total pesticide in soil %

figure('Name','Validation')
plot(time_fin,YPTV0_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV0_95,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,YPTV0_5,YPTV0_95); % Total 
plot(time_fin,YPTV0_50,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off

%% Model V3 %

% Total DNA %

plot(time_fin,YaV3_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV3_95,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV3_5,YaV3_95);
plot(time_fin,YaV3_50,'color',[216/256 179/256 101/256],'LineWidth',2.45)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off

% Total pesticide in soil %

figure('Name','Validation')
plot(time_fin,YPTV3_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV3_95,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,YPTV3_5,YPTV3_95); % Total 
plot(time_fin,YPTV3_50,'color',[27/256 158/256 119/256],'LineWidth',2.45)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off

%% Model V4 %%

% Total DNA %   

plot(time_fin,YaV4_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV4_95,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV4_5,YaV4_95);
plot(time_fin,YaV4_50,'color',[216/256 179/256 101/256],'LineWidth',2.75)
ylim([0 1.5e-3])
hold off
set(gca,'FontSize',20,'yTickLabel',[])

% Total pesticide in soil %

figure('Name','Validation')
plot(time_fin,YPTV4_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV4_95,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,YPTV4_5,YPTV4_95); % Total 
plot(time_fin,YPTV4_50,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off

%% Model V4p %%

% Total DNA %    

figure('Name','Validation')
plot(time_fin,YaV4p_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YaV4p_95,'color','none','LineWidth',1.25)
fill_betweenN(time_fin,YaV4p_5,YaV4p_95);
plot(time_fin,YaV4p_50,'color',[216/256 179/256 101/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off

% Total pesticide in soil %

figure('Name','Validation')
plot(time_fin,YPTV4p_5,'color','none','LineWidth',1.25)
hold on
plot(time_fin,YPTV4p_95,'color','none','LineWidth',1.25)
fill_betweenN2(time_fin,YPTV4p_5,YPTV4p_95); % Total 
plot(time_fin,YPTV4p_50,'color',[27/256 158/256 119/256],'LineWidth',2.75)
ylim([0 1.5e-3])
set(gca,'FontSize',20,'yTickLabel',[])
hold off