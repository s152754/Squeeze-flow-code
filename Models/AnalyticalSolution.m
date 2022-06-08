function [pN,ranaN,HN,RN,t,gammadot] = AnalyticalSolution(eta,F,tend,V,H0,R0)
%% Newtonian
dr = 100;
pN = zeros(tend,dr);
HN = zeros(1,tend);
RN = zeros(1,tend);
hdotN = zeros(1,tend);
ranaN = zeros(tend,dr);
t = zeros(1,tend);
gammadot = zeros(1,tend);
for i = 1:tend*10
    t(i) = i/10;
    HN(i) = H0*(1+((8*H0^2*F*t(i))/(3*pi*eta*R0^4)))^(-1/4);
    RN(i) = sqrt(V/(pi*HN(i)));

    hdotN(i) = (2*HN(i)^3*F)/(3*pi*eta*RN(i)^4);
    gammadot(i) = (3*hdotN(i)*RN(i)/HN(i)^2);
    for j = 1:dr
        rN = j*RN(i)/dr;
        ranaN(i,j) = rN;
        pN(i,j) = 3*eta/HN(i)^3*hdotN(i)*(RN(i)^2-rN^2);
    end
end
xlin = linspace(0.1,600,2);
ylin = linspace(10,10,2)

figure()
loglog(t,gammadot,'b','LineWidth',1)
hold on
loglog(xlin,ylin,'--r','LineWidth',1)
ylabel('$\dot{\gamma}$ [s$^{-1}$]','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
legend({'Newtonian $\dot{\gamma}$','Shear thinnning limit'},'Location','northeast','Interpreter','latex')
ax = gca; 
ax.FontSize = 18;


figure()
plot(t,RN)
end
%% additional figures
%{
% dimension plates (after which the droplet fill fall of the plate)
yPMMA = [5*cm/mm 5*cm/mm];
xPMMA = linspace(t(1),t(end),2);

figure()
plot(t,RN/mm,'LineWidth',1)
hold on
plot(xPMMA, yPMMA,'LineWidth',1)
title('Analytical solution Newtonian fluid')
ylabel('Radius [mm]')
xlabel('Time [s]')
legend('Radius over time','Maximum radius on plate')

figure()
loglog(t,RN/mm)
title('Analytical solution Newtonian fluid')
ylabel('Radius [mm]')
xlabel('Time [s]')
%}

%% Power-law
%{
hdotPL = zeros(1,tend);
pPL = zeros(tend,dr);
HPL = zeros(1,tend);
n = 1; %0-1 shear thinning 1-... shear thickening (we kijken alleen naar shear thinning nu)
K = mu;
for i = 1:tend+1
    tPL(i) = i-1;
    HPL(i) = (abs(1/2*F*(1/n+1)^(-1)*((-1-n)/(1+2*n))*(1/(2*pi*K))^(1/n)*(n+1)^(1/n)*((-1-n)/(2*(n+3)))^(-1/n)*(V/pi)^(-1/2-3/(2*n))*(-3/2-5/(2*n))*tPL(i)+H0^(-(3/2-5/(2*n)))))^(-1/(3/2+5/(2*n)));
    RPL = sqrt(V./(pi*HPL(i)));
    
    hdotPL(i) = 1/2*F^(1/n)*(1/n+1)^(-1)*((-1/n-1)/(1/n+2))*HPL(i)^(1/n+2)*(1/(2*pi*K))^(1/n)*(n+1)^(1/n)*((-n-1)/(2*(n+3))*RPL^(n+3))^(-1/n);
    % wss fout in boundary
    for j = 1:dr
        rPL = j*RPL/dr;
        pPL(i,j) = (1/4*hdotPL(i))^n*(1/n+1)^n*((-1-1/n)/(1/n+2))*(1/2*HPL(i))^(n/(1/n+2))*K*1/(n+1)*(rPL^(n+2)-RPL^(n+2));
    end
end
figure()
plot(t,HN)
hold on
plot(tPL,HPL)
legend('Analytical Newtonian', 'Analytical Power-Law n=1')
%}
