clear all
format long

pi0 = 0.1;
pi1 = 0.15;
pi2 = 1;
pi3 = 2;

% check pi0 < 1
if (pi0 >= 1)
    'pi0 is'
    pi0
    'pi0 too big'
end

Lss = (1+pi0*pi2+0.5*pi1)/(2*pi0*pi3)*(sqrt(1+2*(1-pi0)*pi0*pi3/(1+pi0*pi2+0.5*pi1)^2)-1);
Lcrit = (1+pi0*pi2)/(2*pi0*pi3)*(sqrt(1+4*(1-pi0)*pi0*pi3/(1+pi0*pi2)^2)-1);

Lmincoarse = 0.01;
Ncoarse = 30;

L1coarse = linspace(Lmincoarse,1,Ncoarse);
L2coarse = linspace(Lmincoarse,1,Ncoarse);

[LL1coarse,LL2coarse] = meshgrid(L1coarse,L2coarse);

coarse_inds = find(LL1coarse+LL2coarse < 1);

JJ = 1+pi2*(LL1coarse+LL2coarse)+pi3*(LL1coarse.^2+LL2coarse.^2);
dL1coarse = (1-LL1coarse-LL2coarse)./JJ-pi0-pi1*LL1coarse./JJ;
dL2coarse = (1-LL1coarse-LL2coarse)./JJ-pi0-pi1*LL2coarse./JJ;

start_inds = 1:2;
end_inds = 3:Ncoarse;

LL1coarse_plt1 = [LL1coarse(start_inds,:) LL1coarse(:,start_inds)'];
LL2coarse_plt1 = [LL2coarse(start_inds,:) LL2coarse(:,start_inds)'];
dLL1coarse_plt1 = [dL1coarse(start_inds,:) dL1coarse(:,start_inds)'];
dLL2coarse_plt1 = [dL2coarse(start_inds,:) dL2coarse(:,start_inds)'];

LL1coarse_plt2 = LL1coarse(end_inds,end_inds);
LL2coarse_plt2 = LL2coarse(end_inds,end_inds);
dLL1coarse_plt2 = dL1coarse(end_inds,end_inds);
dLL2coarse_plt2 = dL2coarse(end_inds,end_inds);

coarse_inds_plt1 = find(LL1coarse_plt1+LL2coarse_plt1 < 1);
coarse_inds_plt2 = find(LL1coarse_plt2+LL2coarse_plt2 < 1);

Lmin = 0.01;
N = 100;

L1 = linspace(Lmin,1,N);
L2 = linspace(Lmin,1,N);

[LL1,LL2] = meshgrid(L1,L2);

fine_inds = find(LL1coarse+LL2coarse < 1);

JJ = 1+pi2*(LL1+LL2)+pi3*(LL1.^2+LL2.^2);
dL1 = (1-LL1-LL2)./JJ-pi0-pi1*LL1./JJ;
dL2 = (1-LL1-LL2)./JJ-pi0-pi1*LL2./JJ;

LY = 0.5*(((LL1-Lss)+(LL2-Lss)).^2+pi3/pi1*((LL1-Lss)-(LL2-Lss)).^2);

figure(1);clf;
quiver(LL1coarse_plt1(coarse_inds_plt1),LL2coarse_plt1(coarse_inds_plt1),...
    dLL1coarse_plt1(coarse_inds_plt1),dLL2coarse_plt1(coarse_inds_plt1),2.5,...
    'Color',[0    0.4470    0.7410]);
hold on
quiver(LL1coarse_plt2(coarse_inds_plt2),LL2coarse_plt2(coarse_inds_plt2),...
    dLL1coarse_plt2(coarse_inds_plt2),dLL2coarse_plt2(coarse_inds_plt2),2.5,...
    'Color',[0    0.4470    0.7410]);
% plot(L1(L1 < P),P-L1(L1 < P),'k');
plot(L1,zeros(N,1),'k','LineWidth',2);
plot(zeros(N,1),L2,'k','LineWidth',2);
plot(L1,1-L1,'k','LineWidth',2);

LY(find(LL1+LL2 > 1)) = NaN;
LY2(find(LL1+LL2 > 1)) = NaN;
LYtmp(find(LL1+LL2 > 1)) = NaN;
LYtmp2(find(LL1+LL2 > 1)) = NaN;

levels = [0.01 0.025 0.05 0.1 0.15];
levels = [0.1 0.5 2 3];
[C, h] = contour(LL1,LL2,LY,levels,'ShowText','on','LineWidth',1.5,'LineStyle','--');
clabel(C,h,'FontSize',20)
scatter(Lss,Lss,800,'k*')
scatter(Lcrit,0,100,'rd','MarkerFaceColor','r')
scatter(0,Lcrit,100,'rd','MarkerFaceColor','r')
xlim([-0.01 1]);
ylim([-0.01 1]);
%axis equal

% %solve ODE for boundary region
% params = [A D P];
% hfun = @(x,h) h_ode1(x,h,params);
% [L1out,hout] = ode45(hfun,[P^2 4*P^2],0);
% plot(sqrt(L1out),sqrt(hout),'k','LineWidth',2);
% 
% hfun = @(y,h) h_ode2(y,h,params);
% [L2out,hout] = ode45(hfun,[P^2 4*P^2],0);
% plot(sqrt(hout),sqrt(L2out),'k','LineWidth',2);

figure(2);clf;
pcolor(LL1,LL2,LY);
colorbar
set(gca,'colorscale','log')

%% add corner boundaries

%solve ODE for boundary region
params = [pi0 pi1 pi2 pi3];
dfun = @(t,y) myderiv_back(t,y,params);
L0 = [Lcrit; 0];
tspan = [0 -100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
figure(1)
plot(L1out,L2out,'r','LineWidth',2);

L0 = [0; Lcrit];
tspan = [0 -100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
plot(L1out,L2out,'r','LineWidth',2);
%% now apply severing

mult_fac = 1/(1-Lss);
pi0 = mult_fac*0.1;
pi1 = 0.15;
pi2 = 1/mult_fac;
pi3 = 2/mult_fac^2;

pi0 = 0.2;
pi1 = 0.15;
pi2 = 0.5;
pi3 = 0.5;


% check pi0 < 1
if (pi0 >= 1)
    'pi0 is'
    pi0
    'pi0 too big'
end

Lss2 = 1/mult_fac*(1+pi0*pi2+0.5*pi1)/(2*pi0*pi3)*(sqrt(1+2*(1-pi0)*pi0*pi3/(1+pi0*pi2+0.5*pi1)^2)-1);
Lcrit2 = 1/mult_fac*(1+pi0*pi2)/(2*pi0*pi3)*(sqrt(1+4*(1-pi0)*pi0*pi3/(1+pi0*pi2)^2)-1);

Lmincoarse = 0.01;
Ncoarse = 30;

L1coarse = linspace(Lmincoarse,1,Ncoarse);
L2coarse = linspace(Lmincoarse,1,Ncoarse);

[LL1coarse,LL2coarse] = meshgrid(L1coarse,L2coarse);

coarse_inds = find(LL1coarse+LL2coarse < 1);

JJ = 1+pi2*(mult_fac*LL1coarse+2*LL2coarse)+pi3*((mult_fac*LL1coarse).^2+(mult_fac*LL2coarse).^2);
dL1coarse = 1/mult_fac*((1-mult_fac*LL1coarse-mult_fac*LL2coarse)./JJ-pi0-pi1*mult_fac*LL1coarse./JJ);
dL2coarse = 1/mult_fac*((1-mult_fac*LL1coarse-mult_fac*LL2coarse)./JJ-pi0-pi1*mult_fac*LL2coarse./JJ);

start_inds = 1:2;
end_inds = 3:Ncoarse;

LL1coarse_plt1 = [LL1coarse(start_inds,:) LL1coarse(:,start_inds)'];
LL2coarse_plt1 = [LL2coarse(start_inds,:) LL2coarse(:,start_inds)'];
dLL1coarse_plt1 = [dL1coarse(start_inds,:) dL1coarse(:,start_inds)'];
dLL2coarse_plt1 = [dL2coarse(start_inds,:) dL2coarse(:,start_inds)'];

LL1coarse_plt2 = LL1coarse(end_inds,end_inds);
LL2coarse_plt2 = LL2coarse(end_inds,end_inds);
dLL1coarse_plt2 = dL1coarse(end_inds,end_inds);
dLL2coarse_plt2 = dL2coarse(end_inds,end_inds);

coarse_inds_plt1 = find(LL1coarse_plt1+LL2coarse_plt1 < 0.5);
coarse_inds_plt2 = find(LL1coarse_plt2+LL2coarse_plt2 < 0.5);

coarse_inds_plt3 = find((LL1coarse+LL2coarse >= 0.5) ...
    & (LL1coarse+LL2coarse < 1));

Lmin = 0.01;
N = 100;

L1 = linspace(Lmin,1,N);
L2 = linspace(Lmin,1,N);

[LL1,LL2] = meshgrid(L1,L2);

fine_inds = find(LL1coarse+LL2coarse < 1);

JJ = 1+pi2*(mult_fac*LL1+mult_fac*LL2)+pi3*((mult_fac*LL1).^2+(mult_fac*LL2).^2);
dL1 = 1/mult_fac*((1-mult_fac*LL1-mult_fac*LL2)./JJ-pi0-pi1*mult_fac*LL1./JJ);
dL2 = 1/mult_fac*((1-mult_fac*LL1-mult_fac*LL2)./JJ-pi0-pi1*mult_fac*LL2./JJ);

LY = 0.5*(((2*LL1-2*Lss2)+(2*LL2-2*Lss2)).^2+pi3/pi1*((2*LL1-2*Lss2)-(2*LL2-2*Lss2)).^2);

figure(3);clf;
% quiver(LL1coarse_plt1(coarse_inds_plt1),LL2coarse_plt1(coarse_inds_plt1),...
%     dLL1coarse_plt1(coarse_inds_plt1),dLL2coarse_plt1(coarse_inds_plt1),1.5,...
%     'Color',[0    0.4470    0.7410]);
% hold on
% quiver(LL1coarse_plt2(coarse_inds_plt2),LL2coarse_plt2(coarse_inds_plt2),...
%     dLL1coarse_plt2(coarse_inds_plt2),dLL2coarse_plt2(coarse_inds_plt2),2,...
%     'Color',[0    0.4470    0.7410]);
% quiver(LL1coarse(coarse_inds_plt3),LL2coarse(coarse_inds_plt3),...
%     dL1coarse(coarse_inds_plt3),dL2coarse(coarse_inds_plt3),1,...
%     'Color',[0    0.4470    0.7410]);
% plot(L1(L1 < P),P-L1(L1 < P),'k');
plot(L1,zeros(N,1),'k','LineWidth',2);
hold on
plot(zeros(N,1),L2,'k','LineWidth',2);
plot(L1,1-L1,'k','LineWidth',2);

LY(find(LL1+LL2 > 1)) = NaN;
LY2(find(LL1+LL2 > 1)) = NaN;
LYtmp(find(LL1+LL2 > 1)) = NaN;
LYtmp2(find(LL1+LL2 > 1)) = NaN;

levels = [0.01 0.025 0.05 0.1 0.15];
levels = [0.1 0.5 1];
[C, h] = contour(LL1,LL2,LY,levels,'ShowText','on','LineWidth',1.5,'LineStyle','--');
clabel(C,h,'FontSize',12)
scatter(Lss2,Lss2,800,'k*')
scatter(Lcrit2,0,100,'rd','MarkerFaceColor','r')
scatter(0,Lcrit2,100,'rd','MarkerFaceColor','r')
xlim([-0.01 1]);
ylim([-0.01 1]);
%axis equal

% %solve ODE for boundary region
% params = [A D P];
% hfun = @(x,h) h_ode1(x,h,params);
% [L1out,hout] = ode45(hfun,[P^2 4*P^2],0);
% plot(sqrt(L1out),sqrt(hout),'k','LineWidth',2);
% 
% hfun = @(y,h) h_ode2(y,h,params);
% [L2out,hout] = ode45(hfun,[P^2 4*P^2],0);
% plot(sqrt(hout),sqrt(L2out),'k','LineWidth',2);

figure(4);clf;
pcolor(LL1,LL2,LY);
colorbar
set(gca,'colorscale','log')

%% add corner boundaries

%solve ODE for boundary region
params = [pi0 pi1 pi2 pi3 mult_fac];
dfun = @(t,y) myderiv2_back(t,y,params);
L0 = [Lcrit2; 0];
tspan = [0 -100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
figure(3)
plot(L1out,L2out,'r','LineWidth',2);

L0 = [0; Lcrit2];
tspan = [0 -100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
plot(L1out,L2out,'r','LineWidth',2);

a1 = gca;
f6 = figure(6);clf;
a2 = copyobj(a1,f6);

%% add another severing trajectory for fun

%solve ODE for boundary region
params = [pi0 pi1 pi2 pi3 mult_fac];
dfun = @(t,y) myderiv2(t,y,params);
L0 = [0.7; 0];
L0 = [0.8; 0];
tspan = [0 100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
figure(6)
plot(L1out,L2out,'r','LineWidth',2);

figure(3)
plot(L1out,L2out,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);

figure(7);clf;
plot(tout,L1out,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(tout,L2out,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
% legend('L_1','L_2','Location','Best')
set(gca, 'fontsize', 18);
ylim([0 0.8])

figure(8);clf;
plot(tout,L1out,'Color',[0.4660 0.6740 0.1880],'LineWidth',4);
hold on
plot(tout,L2out,'Color',[0.4660 0.6740 0.1880],'LineWidth',4);
% legend('L_1','L_2','Location','Best')
set(gca, 'fontsize', 34);
xlim([0 20])
ylim([0 0.2])

figure(3);
%% add severing trajectory

%solve ODE for boundary region
params = [pi0 pi1 pi2 pi3 mult_fac];
dfun = @(t,y) myderiv2(t,y,params);
L0 = [Lss; 0];
% L0 = [1.5*Lss; 0];
tspan = [0 100];
[tout,Lout] = ode45(dfun,tspan,L0);
L1out = Lout(:,1); L2out = Lout(:,2);
figure(3)
plot(L1out,L2out,'Color',[0.4940 0.1840 0.5560],'LineWidth',2);

figure(5);clf;
plot(tout,L1out,'Color',[0.4940 0.1840 0.5560],'LineWidth',2);
hold on
plot(tout,L2out,'Color',[0.4940 0.1840 0.5560],'LineWidth',2);
% legend('L_1','L_2','Location','Best')
set(gca, 'fontsize', 18);
ylim([0 0.8])

%% define differential equation

function df = myderiv(t,LL,params)

pi0 = params(1);
pi1 = params(2);
pi2 = params(3);
pi3 = params(4);

L1 = LL(1);
L2 = LL(2);

JJ = 1+pi2*(L1+L2)+pi3*(L1.^2+L2.^2);
dL1 = (1-L1-L2)./JJ-pi0-pi1*L1./JJ;
dL2 = (1-L1-L2)./JJ-pi0-pi1*L2./JJ;

if (L1 <= eps && dL1 <= 0)
    dL1 = 0;
end
if (L2 <= eps && dL2 <= 0)
    dL2 = 0;
end

df = [dL1; dL2];
end

function df = myderiv_back(t,LL,params)

pi0 = params(1);
pi1 = params(2);
pi2 = params(3);
pi3 = params(4);

L1 = LL(1);
L2 = LL(2);

JJ = 1+pi2*(L1+L2)+pi3*(L1.^2+L2.^2);
dL1 = (1-L1-L2)./JJ-pi0-pi1*L1./JJ;
dL2 = (1-L1-L2)./JJ-pi0-pi1*L2./JJ;

if (L1+L2 >= 1-1.6e-2 && dL1+dL2 <= 0)
    dL1 = 0;
    dL2 = 0;
end

df = [dL1; dL2];
end

function df = myderiv2(t,LL,params)

pi0 = params(1);
pi1 = params(2);
pi2 = params(3);
pi3 = params(4);
mult_fac = params(5);

L1 = LL(1);
L2 = LL(2);

% JJ = 1+pi2*(L1+L2)+pi3*(L1.^2+L2.^2);
% dL1 = (1-L1-L2)./JJ-pi0-pi1*L1./JJ;
% dL2 = (1-L1-L2)./JJ-pi0-pi1*L2./JJ;

JJ = 1+pi2*mult_fac*(L1+L2)+pi3*mult_fac^2*(L1.^2+L2.^2);
dL1 = 1/mult_fac*((1-mult_fac*(L1+L2))./JJ-pi0-pi1*mult_fac*L1./JJ);
dL2 = 1/mult_fac*((1-mult_fac*(L1+L2))./JJ-pi0-pi1*mult_fac*L2./JJ);

if (L1 <= eps && dL1 <= 0)
    dL1 = 0;
end
if (L2 <= eps && dL2 <= 0)
    dL2 = 0;
end

df = [dL1; dL2];
end

function df = myderiv2_back(t,LL,params)

pi0 = params(1);
pi1 = params(2);
pi2 = params(3);
pi3 = params(4);
mult_fac = params(5);

L1 = LL(1);
L2 = LL(2);

% JJ = 1+pi2*(L1+L2)+pi3*(L1.^2+L2.^2);
% dL1 = (1-L1-L2)./JJ-pi0-pi1*L1./JJ;
% dL2 = (1-L1-L2)./JJ-pi0-pi1*L2./JJ;

JJ = 1+pi2*mult_fac*(L1+L2)+pi3*mult_fac^2*(L1.^2+L2.^2);
dL1 = 1/mult_fac*((1-mult_fac*(L1+L2))./JJ-pi0-pi1*mult_fac*L1./JJ);
dL2 = 1/mult_fac*((1-mult_fac*(L1+L2))./JJ-pi0-pi1*mult_fac*L2./JJ);

if (L1+L2 >= 1-1e-2 && dL1+dL2 <= 0)
    dL1 = 0;
    dL2 = 0;
end

df = [dL1; dL2];
end
