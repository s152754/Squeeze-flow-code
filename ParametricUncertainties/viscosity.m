function [mu_eta, std_eta] = viscosity(material)

% this is viscosity for newtonian model%

mea = 10;
if strcmp(material,"Glycerol")
    for i = 1:mea
        read = readtable(sprintf("rheometerGly/SRST_TXT_S%01.0f.txt",i+1));
        array = table2array(read(:,2:3));
        data{i} = rmmissing(array);
    
        if i==1
            rads(i,:) = data{i}(:,2);
            eta_size = size(data{i}(:,1));
            eta = zeros(mea,eta_size(1,1));
        end
    
        eta(i,:) = data{i}(:,1);
        name(i) = "dataset " + num2str(i);
%         loglog(rads,eta(i,:),'o','MarkerSize', 10)
%         xlim([0 100])
%         ylim([0.0 1.5])
%         xlabel('$\dot{\gamma}$ [1/s]','Interpreter','latex')
%         ylabel('$\eta$ [Pa*s]','Interpreter','latex')
%         ax = gca; 
%         ax.FontSize = 18;
%         hold on
    end
%         legend(name)
%         figure
   % post-processing
mu_s1 = mean(eta());
std_s1 = std(eta);

% gemiddelde over alle data punten
mu_eta = mean(eta,'all');
std_eta = std(eta,0,'all');

x_mu = logspace(0,2,2)
y_mu = linspace(mu_eta,mu_eta,2)

figure()
errorbar(rads,mu_s1,std_s1,'b','LineWidth',1)
hold on
plot(x_mu,y_mu,'r','LineWidth',1)
plot(x_mu,y_mu+std_eta,'--r','LineWidth',1)
plot(x_mu,y_mu-std_eta,'--r','LineWidth',1)
set(gca,'YScale','log','XScale','log')
ax = gca; 
ax.FontSize = 18;
ylabel('$\eta$ [Pa$\cdot$s]','Interpreter','latex')
xlabel('$\dot{\gamma}$ [s$^{-1}$]','Interpreter','latex')
legend({'$\mu_e \pm 2\sigma_e$','$\mu_{\rm{all}}$','$\mu_{\rm{all}} \pm 2\sigma_{\rm{all}}$'},'Location','northeast','Interpreter','latex')

    
elseif strcmp(material,"Shampoo")
    % create matrices
    etainf      = zeros(1,mea);
    eta0        = zeros(1,mea);
    lambda      = zeros(1,mea);
    a           = zeros(1,mea);
    n           = zeros(1,mea);
    
    for i = 1:mea
        read = readtable(['shearrheometer\Shampoo_21\shampoo_21_SRST_TXT_S' num2str(i+2) '.txt']);
        array = table2array(read(:,2:3));
        data{i} = rmmissing(array);
    
        eta = data{i}(:,1);
        rad = data{i}(:,2);
%         % mu0 determination
        etas(i,:) = eta;
        etas
        mu0(i) = mean(eta(1:5));
% 
%         % power law
        type = fittype('a*x^b');
        options = fitoptions(type);
        options.upper = [1000 10];
        options.Lower = [0 -10];
        [F, gof,opt]=fit(rad(11:16),eta(11:16),type, options);
        RMSE(i) = gof.rmse;
        xx = logspace(1,2,100);
        
        shearrate = xx;
        K(i) = F.a; n(i) = F.b;
        etashear(i,:) = K(i)*shearrate.^n(i);
        pow(i) = 1+n(i);

    end
    rads = rad;
%     % mu0
    mu_eta = mean(etas); std_eta = std(etas);
    directory = 'C:\Users\s152754\OneDrive - TU Eindhoven\Documents\Master\afstuderen\experiment\camera\MonteCarlo';
T = table(rads,mu_eta',std_eta');
table_path = fullfile(directory, "Rheometerdata_rad_mu_sig.txt");
writetable(T,table_path);
    mu0mu = mean(mu0)
    std0 = std(mu0)
    x_mu = linspace(0.1,0.7,2);
    y_mu = linspace(mu0mu,mu0mu,2);
    figure()
    errorbar(rads,mu_eta,2*std_eta,'b','LineWidth',1)
%     hold on
%     plot(x_mu,y_mu,'r','LineWidth',1)
    set(gca,'YScale','log','XScale','log')
    ax = gca; 
    ax.FontSize = 18;
    ylabel('$\eta$ [Pa*s]','Interpreter','latex')
    xlabel('$\dot{\gamma}$ [s$^{-1}$]','Interpreter','latex')
    legend({'$\mu_{0,e} \pm 2\sigma_{0,e}$','$\mu_{all}$'},'Location','northeast','Interpreter','latex')
for i = 1:mea
% carreau
        powf = pow(i)
        mu0f = mu0(i)
        type = fittype('0.001 + (mu0f-0.001) / (1+(a*x)^b)^((1-powf)/2)');
        options = fitoptions(type);
        options.upper = [1 10];
        options.Lower = [0 0];
        [F, gof,opt]=fit(rads,mu_eta',type, options);
        xx = logspace(-1,2,100);
        shearrate = xx;
        lambda(i) = F.a; yas(i) = F.b;
        etashear(i,:) = 0.001 + (mu0f-0.001)./(1+(lambda(i)*shearrate).^yas(i)).^((1-powf)/2);
end
mu_lambda = mean(lambda)
std_lambda = std(lambda)

%carreafit
figure()
    loglog(rads,etas(1,:),'-ob','LineWidth',1)
    hold on
    loglog(xx,etashear(1,:),'r','LineWidth',1)      
    set(gca,'YScale','log','XScale','log')
    ax = gca; 
    ax.FontSize = 18;
    ylabel('$\eta$ [Pa*s]','Interpreter','latex')
    xlabel('$\dot{\gamma}$ [s$^{-1}$]','Interpreter','latex')
    legend({'Measurement data','Carreau fit'},'Location','southwest','Interpreter','latex')

%     % n
%     mu_n = mean(1+n)
%     std_n = std(1+n)
%     figure()
%     loglog(rads,etas(1,:),'-ob','LineWidth',1)
%     hold on
%     loglog(xx,etashear(1,:),'r','LineWidth',1)      
%     set(gca,'YScale','log','XScale','log')
%     ax = gca; 
%     ax.FontSize = 18;
%     ylabel('$\eta$ [Pa*s]','Interpreter','latex')
%     xlabel('$\dot{\gamma}$ [s$^{-1}$]','Interpreter','latex')
%     legend({'Measurement data','Power-law fit'},'Location','northeast','Interpreter','latex')
end
end