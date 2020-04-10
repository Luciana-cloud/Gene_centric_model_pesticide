%% CALIBRATION AND VALIDATION OUTPUTS %%

% This file includes the figures 2 and 3 of the manuscript and the figure
% S2 of the Supporting Information. The data for this file can be found in
% the file Data_file, Outputs_F_24_D and Outputs_F_MCPA, and in the file
% Data_file, Pesticide_Data

%% 2.4-D %%

% CALLING DATA %

load('data_24D.txt')
meas_data               = data_24D(:,2);
meas_data(meas_data==0) = 5e-4;
meas_std                = data_24D(:,3);
days_meas               = data_24D(:,1);

% TIME %

load('time_fin_2,4.mat')

%% MINERALIZATION %%

% Model V0 %
load('sdV0_M_2.mat')% Standard deviation
load('YV0_M_2.mat') % Mean Value
load('YV0_24D.mat')% Median
% Model V3 %
load('sdV3_M_2.mat')% Standard deviation
load('YV3_M_2.mat') % Mean Value
load('YV3_24D.mat')% Median
% Model V4P %
load('sdV4P_M_2.mat')% Standard deviation
load('YV4P_M_2.mat') % Mean Value
load('YV4P_24D.mat')% Median
% Model V4 %
load('sdV4_M_2.mat')% Standard deviation
load('YV4_M_2.mat') % Mean Value
load('YV4_24D.mat')% Median

%% Plot %%

figure('Name','Validation')
B=plot(time_fin,YV3_24D,'color',[227/256 74/256 51/256],'LineWidth',2);   % RED
hold on
C=plot(time_fin,YV4_24D,'color',[49/256 163/256 84/256],'LineWidth',2);% GREEN
D=plot(time_fin,YV4p_24D,'color',[250/256 159/256 181/256],'LineWidth',2);  % MAGENTA
A=plot(time_fin,YV0_24D,'color','black','LineWidth',3);    % BLACK
E = errorbar(days_meas(1:26),meas_data(1:26),meas_std(1:26),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 189/256],'LineWidth',0.75,'Color',[49/256 130/256 189/256]);
hold off
legend([A B C D E],'Full','Opt','Bio 1','Bio 2','Observations','Location','northwest')
box on
ylim([0 100])
set(gca,'FontSize',20)


figure('Name','Validation')
D = scatter(meas_data(1:26),YV3_M_2,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(1:26),YV3_M_2,meas_std(1:26), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:26),YV3_M_2,sdV3_M_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(1:26),YV4_M_2,50,[49/256 163/256 84/256],'filled');
eb(1) = errorbar(meas_data(1:26),YV4_M_2,meas_std(1:26), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:26),YV4_M_2,sdV4_M_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[49/256 163/256 84/256], 'LineWidth', 0.75)
F = scatter(meas_data(1:26),YV4P_M_2,50,[250/256 159/256 181/256],'filled');
eb(1) = errorbar(meas_data(1:26),YV4p_M_2,meas_std(1:26), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:26),YV4p_M_2,sdV4p_M_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[250/256 159/256 181/256], 'LineWidth', 0.75)
A = scatter(meas_data(1:26),YV0_M_2,50,'black','filled');
eb(1) = errorbar(meas_data(1:26),YV0_M_2,meas_std(1:26), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:26),YV0_M_2,sdV0_M_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2)
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1.0;
ylim([0 100])
box on
hold off                                        
set(gca,'FontSize',20,'yTickLabel',[])

%% Transcripts in soil %%

% Model V0 %
load('sdV0_T_2.mat')% Standard deviation
load('YV0_T_2.mat') % Mean Value
load('YV0_T_24D.mat')% Median
% Model V3 %
load('sdV3_T_2.mat')% Standard deviation
load('YV3_T_2.mat') % Mean Value
load('YV3_T_24D.mat')% Median
% Model V4 %
load('sdV4_T_2.mat')% Standard deviation
load('YV4_T_2.mat') % Mean Value
load('YV4_T_24D.mat')% Median

%% Plot %%

figure('Name','Validation')
B=plot(time_fin,YV3_T_24D,'color',[227/256 74/256 51/256],'LineWidth',2);   % red
hold on
C=plot(time_fin,YV4_T_24D,'color',[49/256 163/256 84/256],'LineWidth',2);% green
A=plot(time_fin,YV0_T_24D,'color','black','LineWidth',3);    % blacl
D = errorbar(days_meas(27:39),meas_data(27:39),meas_std(27:39),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 189/256],'LineWidth',0.75,'Color',[49/256 130/256 189/256]);
hold off
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])


figure('Name','Validation')
D = scatter(meas_data(27:39),YV3_T_2,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(27:39),YV3_T_2,meas_std(27:39), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(27:39),YV3_T_2,sdV3_T_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(27:39),YV4_T_2,50,[49/256 163/256 84/256],'filled');
eb(1) = errorbar(meas_data(27:39),YV4_T_2,meas_std(27:39), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(27:39),YV4_T_2,sdV4_T_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [49/256 163/256 84/256], 'LineWidth', 0.75)
A = scatter(meas_data(27:39),YV0_T_2,50,'black','filled');
eb(1) = errorbar(meas_data(27:39),YV0_T_2,meas_std(27:39), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(27:39),YV0_T_2,sdV0_T_2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2)
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1.0;
ylim([0 inf])
hold off
box on
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])
       
%% Genes in soil %%

% Model V0 %
load('sdV0_G.mat')% Standard deviation
load('YV0_G.mat') % Mean Value
load('YV0_G_24D.mat')% Median
% Model V3 %
load('sdV3_G.mat')% Standard deviation
load('YV3_G.mat') % Mean Value
load('YV3_G_24D.mat')% Median
% Model V4 %
load('sdV4_G.mat')% Standard deviation
load('YV4_G.mat') % Mean Value
load('YV4_G_24D.mat')% Median
% Model V4p %
load('sdV4p_G.mat')% Standard deviation
load('YV4p_G.mat') % Mean Value
load('YV4p_G_24D.mat')% Median

%% Plot %%

figure('Name','Validation')
B=plot(time_fin,YV3_G_24D,'color',[227/256 74/256 51/256],'LineWidth',2);   % red
hold on
C=plot(time_fin,YV4_G_24D,'color',[49/256 163/256 84/256],'LineWidth',2);% green
D=plot(time_fin,YV4p_G_24D,'color',[250/256 159/256 181/256],'LineWidth',2);  % magenta
A=plot(time_fin,YV0_G_24D,'color','black','LineWidth',3);    % blue
E = errorbar(days_meas(40:54),meas_data(40:54),meas_std(40:54),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 189/256],'LineWidth',0.75,'Color','b');
ylim([0 2.5e7])
hold off
box on
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])


figure('Name','Validation')
D = scatter(meas_data(40:end),YV3_G,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(40:end),YV3_G,meas_std(40:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(40:end),YV3_G,sdV3_G, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(40:end),YV4_G,50,[49/256 163/256 84/256],'filled');
eb(1) = errorbar(meas_data(40:end),YV4_G,meas_std(40:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(40:end),YV4_G,sdV4_G, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [49/256 163/256 84/256], 'LineWidth', 0.75)
F = scatter(meas_data(40:end),YV4p_G,50,[250/256 159/256 181/256],'filled');
eb(1) = errorbar(meas_data(40:end),YV4p_G,meas_std(40:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(40:end),YV4p_G,sdV4p_G, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[250/256 159/256 181/256], 'LineWidth', 0.75)
A = scatter(meas_data(40:end),YV0_G,50,'black','filled');
eb(1) = errorbar(meas_data(40:end),YV0_G,meas_std(40:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(40:end),YV0_G,sdV0_G, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2)
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1.0;
ylim([0 10e6])
xlim([0 inf])
hold off
box on
set(gca,'FontSize',20,'yTickLabel',[])

%% MCPA %%

% CALLING DATA %

load('data_MCPA.txt')
meas_data = data_MCPA(:,2);
meas_data(meas_data==0) = 5e-4;
meas_std = data_MCPA(:,3);
days_meas = data_MCPA(:,1);

%% Mineralization %%

% Model V0 %
load('sdV0.mat')% Standard deviation
load('YV0.mat') % Mean Value
load('YV0_m.mat')% Median
% Model V3 %
load('sdV3.mat')% Standard deviation
load('YV3.mat') % Mean Value
load('YV3_m.mat')% Median
% Model V4p %
load('sdV4p.mat')% Standard deviation
load('YV4p.mat') % Mean Value
load('YV4p_m.mat')% Median
% Model V4 %
load('sdV4.mat')% Standard deviation
load('YV4.mat') % Mean Value
load('YV4_m.mat')% Median

%% Plotting %%

figure('Name','Validation')
D=plot(time_fin,YV3_m,'color',[227/256 74/256 51/256],'LineWidth',2);  % MAGENTA
hold on
B=plot(time_fin,YV4_m,'color',[49/256 163/256 84/256],'LineWidth',2);   % RED
C=plot(time_fin,YV4p_m,'color',[250/256 159/256 181/256],'LineWidth',2);% GREEN
A=plot(time_fin,YV0_m,'color','black','LineWidth',3);    % BLACK
E = errorbar(days_meas(1:20),meas_data(1:20),meas_std(1:20),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 189/256],'LineWidth',0.75,'Color',[49/256 130/256 189/256]);
ylim([0 100])
hold off
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])


figure('Name','Validation')
D = scatter(meas_data(1:20),YV3,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(1:20),YV3,meas_std(1:20), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:20),YV3,sdV3, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(1:20),YV4,50,[49/256 163/256 84/256],'filled');
eb(1) = errorbar(meas_data(1:20),YV4,meas_std(1:20), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:20),YV4,sdV4, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [49/256 163/256 84/256], 'LineWidth', 0.75)
F = scatter(meas_data(1:20),YF,50,[250/256 159/256 181/256],'filled');
eb(1) = errorbar(meas_data(1:20),YV4p,meas_std(1:20), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:20),YV4p,sdV4p, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[250/256 159/256 181/256], 'LineWidth', 0.75)
A = scatter(meas_data(1:20),YV0,50,'black','filled');
hold on
eb(1) = errorbar(meas_data(1:20),YV0,meas_std(1:20), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(1:20),YV0,sdV0, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2) 
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1;
box on
xlim([0 70])
ylim([0 100])
hold off
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])

%% Transcripts in soil %%

% Model V0 %
load('sdV0_T.mat')% Standard deviation
load('YV0_T.mat') % Mean Value
load('YV0_mT.mat')% Median
% Model V3 %
load('sdV3_T.mat')% Standard deviation
load('YV3_T.mat') % Mean Value
load('YV3_mT.mat')% Median
% Model V4 %
load('sdV4_T.mat')% Standard deviation
load('YV4_T.mat') % Mean Value
load('YV4_mT.mat')% Median

%% Plotting %%

figure('Name','Validation')
B=plot(time_fin,YV3_mT,'color',[227/256 74/256 51/256],'LineWidth',2);   % Red
hold on
C=plot(time_fin,YV4_mT,'color',[49/256 163/256 51/256],'LineWidth',2);% Green
A=plot(time_fin,YV0_mT,'color','black','LineWidth',3);    % Blue
D = errorbar(days_meas(21:34),meas_data(21:34),meas_std(21:34),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 189/256],'LineWidth',0.75,'Color',[49/256 130/256 189/256]);
hold off
set(gca,'FontSize',20)
   

figure('Name','Validation')
D = scatter(meas_data(21:34),YV3_T,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(21:34),YV3_T,meas_std(21:34), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(21:34),YV3_T,sdV3_T, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(21:34),YE_T,50,[49/256 163/256 51/256],'filled');
eb(1) = errorbar(meas_data(21:34),YV4_T,meas_std(21:34), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(21:34),YV4_T,sdV4_T, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[49/256 163/256 51/256], 'LineWidth', 0.75)
A = scatter(meas_data(21:34),YV0_T,50,'black','filled');
hold on
eb(1) = errorbar(meas_data(21:34),YV0_T,meas_std(21:34), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(21:34),YV0_T,sdV0_T, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2)
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1;
box on
hold off
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])

%% Genes in soil %%

% Model V0 %
load('sdV0_d.mat')% Standard deviation
load('YV0_d.mat') % Mean Value
load('YV0_mG.mat')% Median
% Model V3 %
load('sdV3_d.mat')% Standard deviation
load('YV3_d.mat') % Mean Value
load('YV3_mG.mat')% Median
% Model V4 %
load('sdV4_d.mat')% Standard deviation
load('YV4_d.mat') % Mean Value
load('YV4_mG.mat')% Median
% Model V4p %
load('sdV4p_d.mat')% Standard deviation
load('YV4p_d.mat') % Mean Value
load('YV4p_mG.mat')% Median

%% Plot %%

figure('Name','Validation')
B=plot(time_fin,YV3_mG,'color',[227/256 74/256 51/256],'LineWidth',2);   % Red
hold on
C=plot(time_fin,YV4_mG,'color',[49/256 163/256 84/256],'LineWidth',2);% Green
D=plot(time_fin,YV4p_mG,'color',[250/256 159/256 181/256],'LineWidth',2);  % Magenta
A=plot(time_fin,YV0_mG,'color','black','LineWidth',3);    % Blue
E = errorbar(days_meas(35:end),meas_data(35:end),meas_std(35:end),'.','MarkerSize',20,...
    'MarkerEdgeColor',[49/256 130/256 181/256],'LineWidth',0.75,'Color',[49/256 130/256 181/256]);
ylim([0 2.5e7])
hold off
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])
   

figure('Name','Validation')
D = scatter(meas_data(35:end),YV3_d,50,[227/256 74/256 51/256],'filled');
hold on
eb(1) = errorbar(meas_data(35:end),YV3_d,meas_std(35:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(35:end),YV3_d,sdV3_d, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [227/256 74/256 51/256], 'LineWidth', 0.75)
E = scatter(meas_data(35:end),YE_d,50,[49/256 163/256 84/256],'filled');
eb(1) = errorbar(meas_data(35:end),YV4_d,meas_std(35:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(35:end),YV4_d,sdV4_d, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [49/256 163/256 84/256], 'LineWidth', 0.75)
F = scatter(meas_data(35:end),YF_d,50,[250/256 159/256 181/256],'filled');
eb(1) = errorbar(meas_data(35:end),YV4p_d,meas_std(35:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(35:end),YV4p_d,sdV4p_d, 'vertical', 'LineStyle', 'none');
set(eb, 'color',[250/256 159/256 181/256], 'LineWidth', 0.75)
A = scatter(meas_data(35:end),YV0_d,50,'black','filled');
hold on
eb(1) = errorbar(meas_data(35:end),YV0_d,meas_std(35:end), 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(meas_data(35:end),YV0_d,sdV0_d, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'black', 'LineWidth', 2)
A1              = refline(1,0);
A1.Color        = [189/256 189/256 189/256];
A1.LineStyle    = '--';
A1.LineWidth    = 1;
ylim([0 2.5e7])
xlim([0 inf])
box on
hold off    
set(gca,'FontSize',20,'yTickLabel',[],'xTickLabel',[])
