%% MASTER OUTPUTS %%

% This file is tent to run all the models and save all the outputs per
% model in order to make easy plotting the results for the publication.
% This is a template example using the Full model V0 for 24D outputs.

%% MODEL WITH DORMANCY - MODEL V0 %%

% Call the accepted parameters from DREAM %

load('V0_24D.mat')

% Empty variables %

mineral_V0        = {};
DNA_tot_V0        = {};
transcripts_V0    = {};
time_V0           = {};
DNAa_V0           = {};
DNAi_V0           = {};
DNAs_V0           = {};
P_V0              = {}; 
CO2_V0            = {};

% Running the model %

for i = 1:length(D)
    for j = 1:3
        [mineral_V0{end+1},DNA_tot_V0{end+1},transcripts_V0{end+1},...
         time_V0{end+1},DNAa_V0{end+1},DNAi_V0{end+1},DNAs_V0{end+1},...
         P_V0{end+1},CO2_V0{end+1}] = template_simulations24D(D(i,:,j)); 
    end
    for iii = 1:10
        if (i/length(D))*100 == iii*10
            fprintf('\nProgress: ' + string(iii) + '0 %% in Model 1\n\n')
        else
        end
    end
end

save('mineral_V0.mat','mineral_V0')
save('DNA_tot_V0.mat','DNA_tot_V0')
save('transcripts_V0.mat','transcripts_V0')
save('time_V0.mat','time_V0')
save('DNAa_V0.mat','DNAa_V0')
save('DNAi_V0.mat','DNAi_V0')
save('DNAs_V0.mat','DNAs_V0')
save('P_V0.mat','P_V0')
save('CO2_V0.mat','CO2_V0')

%% INTERPOLATION OF OUTPUTS %%

% This part of the file is to get the outputs for the figure. We needed to
% interpolate in order to fit to the measurements, and we selected the
% mean, median and 5 to 95% confidence interval. These are only templates
% for active bacteria of V0 model variant, but we are attaching to complete
% outputs to reproduce the figures of the manuscript. 

%% DNA pools %

% Active Bacteria %

Time_A  = [time_V0];
DNAac   = [DNAa_V0];
vq1     = {};

for i=1:length(Time_A)
    [x, index]  = unique(cell2mat(Time_A(i))); 
    DNAacA      = cell2mat(DNAac(i));
    vq1{end+1}  = (interp1(x,DNAacA(index),time_fin,'linear'))';
end

A1ac   = cell2mat(vq1); % Active

YacV0_5  = prctile(A1ac,5,2);
YacV0_50 = prctile(A1ac,50,2);
YacV0_95 = prctile(A1ac,95,2);

save('YacV0_5.mat','YacV0_5')
save('YacV0_50.mat','YacV0_50')
save('YacV0_95.mat','YacV0_95')

% Inactive Bacteria %

Time_A  = [time_V0];
DNAin   = [DNAi_V0];
vq2     = {};

for i=1:length(Time_A)
    [x, index]  = unique(cell2mat(Time_A(i))); 
    DNAinA      = cell2mat(DNAin(i));
    vq2{end+1}  = (interp1(x,DNAinA(index),time_fin,'linear'))';
end

A1in   = cell2mat(vq2); % Inactive
YinV0_5  = prctile(A1in,5,2);
YinV0_50 = prctile(A1in,50,2);
YinV0_95 = prctile(A1in,95,2);

save('YinV0_5.mat','YinV0_5')
save('YinV0_50.mat','YinV0_50')
save('YinV0_95.mat','YinV0_95')

% Dead Bacteria %

Time_A  = [time_V0];
DNAs    = [DNAs_V0];
vq3     = {};

for i=1:length(Time_A)
    [x, index]  = unique(cell2mat(Time_A(i))); 
    DNAsA       = cell2mat(DNAs(i));
    vq3{end+1}  = (interp1(x,DNAsA(index),time_fin,'linear'))';
end

A1s    = cell2mat(vq3); % Relic
YsV0_5   = prctile(A1s,5,2);
YsV0_50  = prctile(A1s,50,2);
YsV0_95  = prctile(A1s,95,2);

save('YsV0_5.mat','YsV0_5')
save('YsV0_50.mat','YsV0_50')
save('YsV0_95.mat','YsV0_95')