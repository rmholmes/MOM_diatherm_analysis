
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

VG = monmean(VG,2,ndays); % volume in each temperature interval
                          % (centerd on T)
HG = monmean(HG,2,ndays);
V = cat(1,cumsum(VG,'reverse'),0);
H = cat(1,cumsum(HG,'reverse'),0);
HE = rho0*Cp*Te.*V;
HI = H-HE;

Vtot = V(1);
pV_ofT = V/Vtot;

pV = linspace(0,1,1000);
T_ofpV = interp1(pV_ofT+(1:length(pV_ofT))'/1e10,Te,pV,'linear');
T_ofpV(1) = Te(end);
H_ofpV = interp1(Te,H,T_ofpV,'linear');
HE_ofpV = interp1(Te,HE,T_ofpV,'linear');
HI_ofpV = interp1(Te,HI,T_ofpV,'linear');

HCh = cumsum(rho0*Cp*T_ofpV*Vtot*(pV(2)-pV(1)));

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
% $$$ subplot(2,2,1);
% $$$ plot(Te,V);
% $$$ ylabel('Ocean Volume $V(\Theta)$ ($m^3$)');
% $$$ xlabel('Temperature $\Theta$ ($^\circ$C)');
% $$$ grid on;
% $$$ xlim([-2 34]);
subplot(3,1,2);
plot(H,Te);
hold on;
plot(HE,Te,'-r');
plot(HI,Te,'-b');
ylim([-2 34]);
xlim([-1 3.5]*1e25);
grid on;
legend('$H(\Theta)$','$H_E(\Theta)$','$H_I(\Theta)$');
xlabel('Heat Content $H(\Theta)$ (J)');
ylabel('Temperature $\Theta$ ($^\circ$C)');

subplot(3,1,1);
plot(T_ofpV,pV);
ylabel('Volume Fraction $p=\mathcal{V}/\mathcal{V}_T$');
xlabel('Temperature $\Theta(p)$ ($^\circ$C)');
set(gca,'ydir','reverse');
grid on;
xlim([-2 34]);

subplot(3,1,3);
plot(H_ofpV/rho0/Cp/Vtot./pV,pV);
hold on;
plot(HE_ofpV,pV,'-r');
plot(HI_ofpV,pV,'-b');
% $$$ plot(HCh,pV,'--m');
set(gca,'ydir','reverse');
xlim([-1 3.5]*1e25);
grid on;
legend('$H(p)$','$H_E(p)$','$H_I(p)$');
ylabel('Volume Fraction $p$');
xlabel('Heat Content (J)');

% $$$ gname = [base 'ocean_grid.nc'];
% $$$ area = ncread(gname,'area_t');
% $$$ ht = ncread(gname,'ht');
% $$$ max(max(ht))
% $$$ dz = 10;
% $$$ z = -5500:dz:0;
% $$$ A = zeros(size(z));
% $$$ for zi=1:length(z)
% $$$     A(zi) = nansum(nansum(area(ht>-z(zi))));
% $$$ end
% $$$ Voz = cumsum(A*dz);%,'reverse');
% $$$ Voz = Voz(end)-Voz;
% $$$ 
% $$$ Voz(1) = Voz(1)+1e12;
% $$$ zT = zeros(size(V));
% $$$ for Ti=1:TL
% $$$     zT(Ti) = interp1(Voz,z,V(Ti),'linear');
% $$$ end
% $$$ 
% $$$ dzTdT = [0; diff(zT)/dT; 0];
% $$$ 
% $$$ AT = [A(1); -diff(V)/dT./dzTdT(2:end-1); A(end)];

% $$$ load('MOM025_kb3seg_output101-110_GlobalBudget_Gterms.mat');
% $$$ 
% $$$ JSint = rho0*Cp*cumsum(JS*dT,'reverse');
% $$$ dVdtint = rho0*Cp*cumsum(dVdt*dT,'reverse');
% $$$ RES = dVdtint-JSint-M-I-F;

