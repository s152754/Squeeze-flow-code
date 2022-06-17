clear all; close all; warning off ; clc
x = input('Experiment: ');
% compare the videos of one experiment
mm = 1e-3; g = 1e-3; mL = 1e-6;

%% cases Glycerol
directoryExp = "C:\Users\s152754\OneDrive - TU Eindhoven\Documents\PhD\Project\Squeezeflow_Matlab\Squeeze-flow-code\Experiments\output\";

switch x
    case "GlyV01F"
        filename = "PU_GlyV01F_";
        Material = "Glycerol";
        
        IR = importdata(directoryExp+"RoutExp_GlyV01F.txt");
        muR0 = IR.data(1,end-1); stdR0 = IR.data(1,end);

        MassTot = importdata('MassGlycerol.txt');
        murho = 1.26*g/mL;
        [muV, stdV] = VolumeExperiment(MassTot,murho);

        ForceTot = importdata("ForceF.txt");
        [muF, stdF] = Force(ForceTot);

        [mueta, stdeta] = viscosity(Material);

    case "GlyV01F025"
        filename = "PU_GlyV01F025_";
        Material = "Glycerol";

        IR = importdata(directoryExp+"RoutExp_GlyV01F025.txt");
        muR0 = IR.data(1,end-1); stdR0 = IR.data(1,end);

        MassTot = importdata('MassGlycerol.txt');
        murho = 1.26*g/mL;
        [muV, stdV] = VolumeExperiment(MassTot,murho);

        ForceTot = importdata("ForceF025.txt");
        [muF, stdF] = Force(ForceTot);
        
        [mueta, stdeta] = viscosity(Material);

    case "GlyV01F05"
        filename = "PU_GlyV01F05_";
        Material = "Glycerol";

        IR = importdata(directoryExp+"RoutExp_GlyV01F05.txt");
        muR0 = IR.data(1,end-1); stdR0 = IR.data(1,end);

        MassTot = importdata('MassGlycerol.txt');
        murho = 1.26*g/mL;
        [muV, stdV] = VolumeExperiment(MassTot,murho);

        ForceTot = importdata("ForceF05.txt");
        [muF, stdF] = Force(ForceTot);
        
        [mueta, stdeta] = viscosity(Material);

    case "GlyV02F025"
        filename = "PU_GlyV02F025_";
        Material = "Glycerol";

        IR = importdata(directoryExp+"RoutExp_GlyV02F025.txt");
        muR0 = IR.data(1,end-1); stdR0 = IR.data(1,end);

        MassTot = importdata('MassGlycerolV02.txt');
        murho = 1.26*g/mL;
        [muV, stdV] = VolumeExperiment(MassTot,murho);

        ForceTot = importdata("ForceF025.txt");
        [muF, stdF] = Force(ForceTot);

        [mueta, stdeta] = viscosity(Material);
    case "ShaV03F025"
        filename = "PU_ShaV03F025_";
        Material = "Shampoo";

        MassTot = importdata('MassShampoo.txt');
        [murho, stdrho] = densityShampoo();
        [muV, stdV] = VolumeExperiment(MassTot,murho);

        IR = importdata('InitialRadius/Radius_ShaV03F025.txt');
        IR.data(5,:) = NaN; % wrongful data
        muR0 = nanmean(IR.data(:,2)); stdR0 = nanstd(IR.data(:,2));

        ForceTot = importdata("ForceF025.txt");
        [muF, stdF] = Force(ForceTot);

        [mueta, stdeta] = viscosity(Material);
     
    otherwise
        disp('does not work')
end

% surface tension
if strcmp(Material,"Glycerol")
    musurf  = mean([42.53 47.39 46.17]*1e-3);
    stdsurf = std([42.53 47.39 46.17]*1e-3);
    factlow = 1; factup = 2;
elseif strcmp(Material,"Shampoo")
    musurf  = mean([24.22 26.25 26.63]*1e-3);
    stdsurf = std([24.22 26.25 26.63]*1e-3);
    factlow = 1; factup = 2;
end


%% create distrution of uncertainties
nsamples = 3000;
syms mu stdn
mulog = log(mu)-0.5*log(1+(stdn/mu)^2);
stdlog = sqrt(log(1+(stdn/mu)^2));

% welke input file wil je maken --> cases maken
Fdist   = makedist('normal','mu',muF,'sigma',stdF);                   
R0dist  = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[muR0 stdR0])),'sigma',eval(subs(stdlog,[mu stdn],[muR0 stdR0])));             
Vdist   = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[muV stdV])),'sigma',eval(subs(stdlog,[mu stdn],[muV stdV])));                 

F       = random(Fdist,nsamples,1);
R0      = random(R0dist,nsamples,1);
V       = random(Vdist,nsamples,1);

m = input('model: ');
switch m
    case "Analytical"
        etadist     = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[mueta stdeta])),'sigma',eval(subs(stdlog,[mu stdn],[mueta stdeta]))); 
        
        eta         = random(etadist,nsamples,1);
        T           = table(F,R0,V,eta);
        filename2   = "analytical";

    case "LaPlace"
        surfdist    = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[musurf stdsurf])),'sigma',eval(subs(stdlog,[mu stdn],[musurf stdsurf])));
        factdist    = makedist('uniform','lower',4,'upper',8);
        etadist     = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[mueta stdeta])),'sigma',eval(subs(stdlog,[mu stdn],[mueta stdeta]))); 
        
        eta         = random(etadist,nsamples,1);
        surf        = random(surfdist,nsamples,1);
        fact        = random(factdist,nsamples,1);
        T           = table(F,R0,V,eta,surf,fact);
        filename2   = "laplace";
    
    case "Carreau"
        ndist       = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[mun stdnn])),'sigma',eval(subs(stdlog,[mu stdn],[mun stdnn])));
        lambdadist  = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[mulambda stdlambda])),'sigma',eval(subs(stdlog,[mu stdn],[mulambda stdlambda])));
        eta0dist    = makedist('lognormal','mu',eval(subs(mulog,[mu stdn],[mueta0 stdeta0])),'sigma',eval(subs(stdlog,[mu stdn],[mueta0 stdeta0])));

        n           = random(ndist,nsamples,1);
        l           = random(lambdadist,nsamples,1);
        eta0        = random(eta0dist,nsamples,1);
        T           = table(F,R0,V,n,l,eta0);
        filename2   = "carreau_";
end
filename3           = num2str(nsamples); 


%% Saving the uncertainty input parameter file
directory = 'C:\Users\s152754\OneDrive - TU Eindhoven\Documents\PhD\Project\Squeezeflow_Matlab\Squeeze-flow-code\ParametricUncertainties';
table_path = fullfile(directory, filename+filename2+filename3+".txt");
writetable(T,table_path);
rng('default');  % For reproducibility

%% plots report
%{
mu = muF; stds = stdF; dist = Fdist;
figure()
x = linspace(mu-4*stds,mu+4*stds,100);
y = pdf(dist,x);
figure()
plot(x,y,'b','LineWidth',1)
xlabel('$F$ [N]','Interpreter','latex')
ylabel('$P_F$ [N$^{-1}$]','Interpreter','latex')
xlim([mu-4*stds mu+4*stds])
ax = gca; 
ax.FontSize = 18;

lower = 4; upper = 8; dist=factdist;
figure()
x = linspace(lower-5,upper+5,200);
y = pdf(dist,x);
figure()
plot(x,y,'b','LineWidth',1)
xlabel('$\alpha$ [-]','Interpreter','latex')
ylabel('$P_{\alpha}$ [-]','Interpreter','latex')
xlim([lower-5 upper+5])
ylim([0 0.3])
ax = gca; 
ax.FontSize = 18;
%}