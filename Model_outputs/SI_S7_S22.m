%% POSTPROCESSING DREAM %%

% This file is a template for post-processing the DREAM outputs. We are
% presenting here the example of the full model V0 DREAM output against
% 2,4-D, and it has to be adapted for the other model versions. The 
% figures here correspond to the Supporting Information S7 to S22.

%% CALLING DATA %%

% We are calling the selected parameter esemble from DREAM that meet the
% criteria of Rˆ-diagnostic (5). The corresponding outputs can be found in
% Data_file, DREAM_outputs. 

load('V0_24D.mat')

%% MARGINAL DISTRIBUTIONS %%
          
figure(1)
X_1 = D5(:,1,:);
histogram(X_1);
xlim([-4 3])
xlabel('f_T - log');
ylabel('log');
saveas(gcf,'X_1eM','epsc')


figure(2)
X_2 = D2(:,2,:);
histogram(X_2);
% xlim([1 10])
xlabel('n_H');
saveas(gcf,'X_2cM','epsc')

figure(2)
X_1 = D6(:,1,:);
histogram(X_1);
% xlim([-10 6])
xlabel('K_G - log');
% ylabel('log');
saveas(gcf,'X_1Mf','epsc')

figure(3)
X_2 = D6(:,2,:);
histogram(X_2);
% xlim([-4 5])
xlabel('\mu_{max} - log');
% ylabel('log');
saveas(gcf,'X_2Mf','epsc')

figure(4)
X_3 = D6(:,3,:);
histogram(X_3);
% xlim([-14 -8])
xlabel('f_1 - log');
% ylabel('log');
saveas(gcf,'X_3Mf','epsc')

figure(6)
X_6 = D2(:,6,:);
histogram(X_6);
% xlim([-8 4])
xlabel('K_M - log');
% ylabel('log');
saveas(gcf,'X_6cM','epsc')

figure(7)
X_7 = D2(:,7,:);
histogram(X_7);
% xlim([-10 -2])
xlabel('C_T - log');
% ylabel('log');
saveas(gcf,'X_7cM','epsc')

figure(5)
X_4 = D6(:,4,:);
histogram(X_4);
xlim([-5 2])
xlabel('a_a - log');
ylabel('log');
saveas(gcf,'X_4Mf','epsc')

figure(9)
X_9 = D2(:,9,:);
histogram(X_9);
% xlim([-7 -2])
xlabel('a_i - log');
% ylabel('log');
saveas(gcf,'X_9cM','epsc')

% figure(3)

figure(10)
X_10 = D2(:,10,:);
histogram(X_10);
% xlim([-5 2])
xlabel('k_r - log');
% ylabel('log');
saveas(gcf,'X_10cM','epsc')

figure(11)
X_11 = D2(:,11,:);
histogram(X_11);
% xlim([-5 2])
xlabel('k_d - log');
% ylabel('log');
saveas(gcf,'X_11cM','epsc')

figure(12)
X_12 = D(:,12,:);
histogram(X_12);
% xlim([-4 4])
xlabel('a_s - log');
% ylabel('log');
saveas(gcf,'X_12M','epsc')

figure(6)
X_5 = D6(:,5,:);
histogram(X_5);
% xlim([0.1 0.9])
xlabel('Y_P');
saveas(gcf,'X_5Mf','epsc')

figure(7)
X_6 = D6(:,6,:);
histogram(X_6);
% xlim([0.1 0.9])
xlabel('a_{CO2}');
saveas(gcf,'X_6Mf','epsc')

figure(8)
X_7 = D6(:,7,:);
histogram(X_7);
% xlim([-2 0])
xlabel('K_{FP} - log');
% ylabel('log');
saveas(gcf,'X_7Mf','epsc')

figure(9)
X_8 = D6(:,8,:);
histogram(X_8);
% xlim([0.8 1])
xlabel('n_{FP}');
saveas(gcf,'X_8Mf','epsc')