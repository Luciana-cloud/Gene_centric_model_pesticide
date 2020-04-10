%% NORMALIZATION %%

% This file is used to normalize the outputs from the Morris method as
% product of the file: "workflow_eet_V0.m"

out      = 14;  % Number of trajectory outputs (see "V0_dormax.m")

% Creating empty matrices for the number of trajectory outputs:

mi_max      = zeros(out,1);
sigma_max   = zeros(out,1);

% Selecting the maximum value of each trajectory output:

for i = 1:out
mi_max(i) = max(mi(i,:));
sigma_max(i) = max(sigma(i,:));
end

% Normalize:

mi_nor = mi./mi_max;
sigma_nor = sigma./sigma_max;

% mi_nor(12,:)    = [];
% sigma_nor(12,:) = [];

% Determine the euclidean distance between mi and sigma nor. Here I get only 
% value which is not what I really want:

normf = sqrt(sum((mi_nor - sigma_nor) .^ 2));

% Determine the euclidean distance between mi and sigma nor for each
% trajectory output:

normA = zeros(out,20);

for j = 1:out
    for k = 1:20
normA(j,k) = sqrt(sum((mi_nor(j,k) - sigma_nor(j,k)) .^ 2));
    end
end

% Plot the distances in a heat map:

ylabels = {'Minimum AT','Maximum HY','Maximum DEA',...
           'Maximum DIA','Maximum CA','Maximum CO2',...
           'Minimum AT in relation to total AT','Maximum HY in relation to total AT','Maximum DEA+DIA in relation to total AT','Minimum DIA in relation to total AT',...
           'Maximum CA in relation to total AT','Maximum CO2 in relation to total AT',...
           'Minimum HY in relation to total AT','Minimum CA in relation to total AT'};
xlabels = {'dAT_HY','dAT_DD','dHY_CYA','dDEA_CYA','dDIA_CYA','dCYA_CO2','KF_AT','KF_HY','KF_DEA','KF_DIA',...
           'KF_CYA','nF_AT','nF_HY','nF_DEA','nF_DIA','nF_CYA','Ki','c_NO3','f3','ATo'};
       
       map = [1 1 1
    0 1 1
    0 1 0
    1 1 0
    1 0 0];

       g = colormap(hot);
       h = heatmap(xlabels,ylabels,normA,'Colormap',g);
       set(h,'Colormap',flipud(g))

