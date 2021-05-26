
% Play around with areas:

base = '/srv/ccrc/data03/z3500785/mom/';
model = 'MOM025_kb3seg';
outputs = 101:110;

load([base 'mat_data/' model sprintf('_output%03d_BaseVars.mat',outputs(1))]);

[xL,yL] = size(lon);
VG = zeros(TL,12);
HG = zeros(TL,12);

for output=outputs
    load([base 'mat_data/' model sprintf('_output%03d_VHza.mat',outputs(1))]);
    VG = VG + squeeze(nansum(V,1));
    HG = HG + squeeze(nansum(H,1));
end
VG = VG/length(outputs);
HG = HG/length(outputs);
clear V H;

VG = monmean(VG,2,ndays);
HG = monmean(HG,2,ndays);

gname = [base 'ocean_grid_mom025.nc'];

area = ncread(gname,'area_t');
ht = ncread(gname,'ht');

max(max(ht))

dz = 10;
z = -5500:dz:0;
A = zeros(size(z));

for zi=1:length(z)
    A(zi) = nansum(nansum(area(ht>-z(zi))));
end

Voz = cumsum(A*dz);%,'reverse');
Voz = Voz(end)-Voz;
V = cumsum(VG,'reverse');

Voz(1) = Voz(1)+1e12;
zT = zeros(size(V));
for Ti=1:TL
    zT(Ti) = interp1(Voz,z,V(Ti),'linear');
end

dzTdT = [0; diff(zT)/dT; 0];

AT = [A(1); -diff(V)/dT./dzTdT(2:end-1); A(end)];

load('MOM025_kb3seg_output101-110_GlobalBudget_Gterms.mat');

JSint = rho0*Cp*cumsum(JS*dT,'reverse');
dVdtint = rho0*Cp*cumsum(dVdt*dT,'reverse');
RES = dVdtint-JSint-M-I-F;

figure;
% $$$ subplot(2,2,1);
% $$$ pcolPlot(lon,lat,area);
% $$$ subplot(2,2,2);
% $$$ pcolPlot(lon,lat,ht);
% $$$ subplot(2,2,3);
subplot(2,4,1);
plot(A,z,'-k');
xlabel('A(z) (m$^2$)');
ylabel('z (m)');
ylim([-5500 0]);
subplot(2,4,2);
plot(Voz,z,'-k');
xlabel('V(z) (m$^3$)');
ylabel('z (m)');
ylim([-5500 0]);
subplot(2,4,3);
plot(V,T,'-k');
xlabel('$\mathcal{V}(\Theta)$ (m$^3$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,4,4);
plot(zT,T,'-k');
xlabel('$z_\Theta(\Theta)$ (m)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,4,5);
plot(dzTdT,Te,'-k');
xlabel('$\frac{\partial z_\Theta}{\partial\Theta}(\Theta)$ (m$^\circ$C$^{-1}$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,4,6);
plot(I/1e15,Te,'-b');
hold on;
plot(M/1e15,Te,'-r');
plot(F/1e15,Te,'-k');
plot(JSint/1e15,Te,'-','color',[0 0.5 0]);
plot(dVdtint/1e15,Te,'-m');
plot(RES/1e15,Te,'--k');
xlabel('Diathermal fluxes (PW)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
legend('$\mathcal{I}$','$\mathcal{M}$','$\mathcal{F}$',...
       '$\rho_0 C_p\int \mathcal{J}_S d\Theta$', ...
       '$\rho_0 C_p\int \frac{d\mathcal{V}}{dt} d\Theta$', ...
       'Residual');
subplot(2,4,7);
g = 9.81;
alpha = 0.2e-3;
fac = -g*alpha/Cp*dzTdT*dT;
plot(fac.*I/1e12,Te,'-b');
hold on;
plot(fac.*M/1e12,Te,'-r');
plot(fac.*F/1e12,Te,'-k');
plot(fac.*JSint/1e12,Te,'-','color',[0 0.5 0]);
plot(fac.*dVdtint/1e12,Te,'-m');
plot(fac.*RES/1e12,Te,'--k');
xlabel('Integrand (TW)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
box on;
It = nansum(-g*alpha/Cp*dT*dzTdT.*I)/1e12;
Mt = nansum(-g*alpha/Cp*dT*dzTdT.*M)/1e12;
Ft = nansum(-g*alpha/Cp*dT*dzTdT.*F)/1e12;
Jt = nansum(-g*alpha/Cp*dT*dzTdT.*JSint)/1e12;
Vt = nansum(-g*alpha/Cp*dT*dzTdT.*dVdtint)/1e12;
Rt = nansum(-g*alpha/Cp*dT*dzTdT.*RES)/1e12;
text(0.8,30,'Contributions to dBPE/dt');
text(0.6,25,sprintf('Total Numerical Mixing = %3.2f TW',It),'color','b');
text(0.6,20,sprintf('Total Vertical Mixing = %3.2f TW',Mt),'color','r');
text(0.6,15,sprintf('Total Surface Forcing = %3.2f TW',Ft),'color','k');
text(0.6,10,sprintf('Total Surface Volume = %3.2f TW',Jt),'color',[0 0.5 0]);
text(0.6,5,sprintf('Total Tendency = %3.2f TW',Vt),'color','m');
text(0.6,0,sprintf('Residual = %3.2f TW',Rt),'color','k');

% Effective diffusivity:
figure;
% $$$ subplot(2,3,1);
% $$$ plot(A,z,'-k');
% $$$ xlabel('A(z) (m$^2$)');
% $$$ ylabel('z (m)');
% $$$ ylim([-5500 0]);
subplot(2,3,1);
plot(AT,Te,'-k');
xlabel('A($\Theta$) (m$^2$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,3,2);
plot(V,T,'-k');
xlabel('$\mathcal{V}(\Theta)$ (m$^3$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,3,3);
plot(dzTdT,Te,'-k');
xlabel('$\frac{\partial z_\Theta}{\partial\Theta}(\Theta)$ (m$^\circ$C$^{-1}$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,3,4);
plot(dzTdT./AT,Te,'-k');
xlabel('$\frac{1}{A(\Theta)}\frac{\partial z_\Theta}{\partial\Theta}$ ($^\circ$C$^{-1}$ m$^{-1}$)');% (m$^{-1}$$^\circ$C$^{-1}$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
subplot(2,3,5);
plot(I/1e15,Te,'-b');
hold on;
plot(M/1e15,Te,'-r');
% $$$ plot(F/1e15,Te,'-k');
% $$$ plot(JSint/1e15,Te,'-','color',[0 0.5 0]);
% $$$ plot(dVdtint/1e15,Te,'-m');
% $$$ plot(RES/1e15,Te,'--k');
xlabel('Diathermal fluxes (PW)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
xlim([-1.5 0]);
legend('$\mathcal{I}$','$\mathcal{M}$','$\mathcal{F}$',...
       '$\rho_0 C_p\int \mathcal{J}_S d\Theta$', ...
       '$\rho_0 C_p\int \frac{d\mathcal{V}}{dt} d\Theta$', ...
       'Residual');
subplot(2,3,6);
plot(-I.*dzTdT/rho0/Cp./AT,Te,'-b');
hold on;
plot(-M.*dzTdT/rho0/Cp./AT,Te,'-r');
% $$$ plot(-F.*dzTdT/rho0/Cp./AT,Te,'-k');
% $$$ plot(-JSint.*dzTdT/rho0/Cp./AT,Te,'-','color',[0 0.5 0]);
% $$$ plot(-dVdtint.*dzTdT/rho0/Cp./AT,Te,'-m');
% $$$ plot(-RES.*dzTdT/rho0/Cp./AT,Te,'--k');
xlabel('$\kappa_{eff}$ (m$^2$s$^{-1}$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
box on;
set(gca,'xscale','log');
xlim([1e-6 1e-3]);
grid on;
set(gca,'xtick',[1e-6 1e-5 1e-4 1e-3]);

% Add in actual diffusivity average:
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
load([base model sprintf('_output%03d_diff_cbt_T.mat',121)]);%outputs(1))]);
diff_cbt_G = nanmonmean(diff_cbt_T_G,3,ndays);

KdiffR = squeeze(nanmean(diff_cbt_G,1));
[X,Y] = ndgrid(yt,T);

figure;


figure;
plot(-I.*dzTdT/rho0/Cp./AT,Te,'-b');
hold on;
plot(-M.*dzTdT/rho0/Cp./AT,Te,'-r');
plot(KdiffR,T,'--r');
legend('Numerical Mixing','Vertical Mixing','Model diffusivity');
xlabel('$\kappa_{eff}$ (m$^2$s$^{-1}$)');
ylabel('$\Theta$ ($^\circ$C)');
ylim([-2 34]);
box on;
set(gca,'xscale','log');
xlim([1e-6 1e-3]);
grid on;
set(gca,'xtick',[1e-6 1e-5 1e-4 1e-3]);



