function [muV, stdV] = VolumeExperiment(MassTot,rho)
% volume glycerol
g = 1e-3; mL = 1e-6; mg = 1e-6;
Mass = (MassTot.data(:,3)-MassTot.data(:,2))*mg;
VolG = Mass/rho ;
muV = mean(VolG); stdV = std(VolG);
end