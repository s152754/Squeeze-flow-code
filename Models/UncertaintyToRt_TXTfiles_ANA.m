clear all; close all; clc;
n = input('Experiment:')
% write uncertainty to radius in txt files
directoryIN = "C:\Users\s152754\OneDrive - TU Eindhoven\Documents\PhD\Project\Squeezeflow_Matlab\Squeeze-flow-code\ParametricUncertainties\";
directoryOUT = "C:\Users\s152754\OneDrive - TU Eindhoven\Documents\PhD\Project\Squeezeflow_Matlab\Squeeze-flow-code\Models\";

directoryExp = "C:\Users\s152754\OneDrive - TU Eindhoven\Documents\PhD\Project\Squeezeflow_Matlab\Squeeze-flow-code\Experiments\output\";

switch n
    case "GlyV01F"
        filename = "PU_GlyV01F_laplace3000.txt";    
        IR = importdata(directoryExp+"RoutExp_GlyV01F.txt");
        tend = IR.data(end,1);
    case "GlyV01F025"
        filename = "PU_GlyV01F025_laplace3000.txt";
        IR = importdata(directoryExp+"RoutExp_GlyV01F025.txt");
        tend = IR.data(end,1);
    case "GlyV01F05"
        filename = "PU_GlyV01F05_laplace3000.txt";  
        IR = importdata(directoryExp+"RoutExp_GlyV01F05.txt");
        tend = IR.data(end,2);
    case "GlyV02F025"
        filename = "PU_GlyV02F025_laplace3000.txt";  
        IR = importdata(directoryExp+"RoutExp_GlyV02F025.txt");
        tend = IR.data(end,1);
    case "ShaV03F025"
        filename = "MC_Shampoo_V03_F025_3000_nosurf_NewtonianData2_mu0.txt";          
    otherwise
        disp('does not work')
end
tic
MC      = importdata(directoryIN+filename);
dt      = 0.5;

% dependent variables
t = 0:dt:tend;       
for j = 1:length(MC.data(:,1))
    F       = MC.data(j,1);
    R0      = MC.data(j,2);
    V       = MC.data(j,3);
    eta     = MC.data(j,4);
    H0 = V./(R0.^2*pi);
    Hn = H0; Rn = R0;
    for i = 1:length(t)
      HN= H0*(1+((8*H0^2*F*t(i))/(3*pi*eta*R0^4)))^(-1/4);
      RN(i,j) = sqrt(V/(pi*HN));
    end
end
muRN    = mean(RN,2);
stdRN   = std(RN,0,2);

outputname = char(filename);
outputname = outputname(1:end-4);
outputname = string(outputname) + "_out.txt";
writematrix([t',RN,muRN,stdRN],directoryOUT + outputname); %writematrix([t',RN],outputname);%
toc