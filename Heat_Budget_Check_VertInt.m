% Check vertically-integrated budgets in different runs

addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));

% $$$ % MOM-SIS025:
% $$$ baseL = '/short/e14/rmh561/mom/archive/';
% $$$ model = 'MOM025_kb3seg';
% $$$ baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
% $$$ output = 96;
% $$$ post = ''; % For MOM-SIS.
% $$$ haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
% $$$ haveGM = 0; % 1 = GM is on, 0 = off;
% $$$ haveMDS = 0; % 1 = MDS is on, 0 = off;
% $$$ haveSIG = 0; % 1 = SIG is on, 0 = off;
% $$$ 
% $$$ % AOM2-025:
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ model = 'ACCESS-OM2_025deg_jra55_ryf8485_gmredi6';
% $$$ baseD = [baseL '025deg_jra55_ryf8485_gmredi6/']; %Data Directory.
% $$$ output = 148;
% $$$ post = 'ocean/'; % For ACCESS-OM2 output coulpled;
% $$$ haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
% $$$ haveGM = 1; % 1 = GM is on, 0 = off;
% $$$ haveMDS = 0; % 1 = MDS is on, 0 = off;
% $$$ haveSIG = 0; % 1 = SIG is on, 0 = off;
% $$$ 
% AOM2-1:
baseL = '/short/e14/rmh561/access-om2/archive/';
model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_july';
baseD = [baseL '1deg_jra55_ryf8485_kds50_july/']; %Data Directory.
output = 38;
post = 'ocean/'; % For ACCESS-OM2 output coulpled;
haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
haveGM = 1; % 1 = GM is on, 0 = off;
haveMDS = 1; % 1 = MDS is on, 0 = off;
haveSIG = 1; % 1 = SIG is on, 0 = off;

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

% term options:
haveSUB = 1; % 1 = submeso is on, 0 = off;

% scaling constant on the transports:
if (strcmp(model(1),'A')) %ACCESS-OM2, transport in kg/s
    tsc = 1;
else % MOM-SIS, transport in 1e9 kg/s
    tsc = 1e9;
end
tsc = 1e9;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
baser = [baseD sprintf('restart%03d/',output-1) post];
gname = [base 'ocean_grid.nc'];
wnameI = [base 'ocean_wmass_integrated.nc'];
wname = [base 'ocean_wmass.nc'];
fname = [base 'ocean.nc'];
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
xt = ncread(gname,'xt_ocean');xu = ncread(gname,'xu_ocean');
yt = ncread(gname,'yt_ocean');yu = ncread(gname,'yu_ocean');

% Time  -----------------------------------------
ndays = ncread(wname,'average_DT');
tL = length(ndays);

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% Mask ------------------------------------------
mask = ncread(fname,'temp',[1 1 1 1],[xL yL 1 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

%%%% Check numdif terms:
% $$$ [tmp ind1] = min(abs(xt+220));
% $$$ [tmp ind2] = min(abs(xt+140));
% $$$ [tmp ind3] = min(abs(yt+20));
% $$$ [tmp ind4] = min(abs(yt-20));
ind1 = 1; ind2=xL; ind3=1;ind4=yL;

xvec = ind1:ind2;
yvec = ind3:ind4;
xL = ind2-ind1+1;
yL = ind4-ind3+1;
doSGMviac = 0;

clear SS;

SS.JS = monmean(ncread(wnameI,'mass_pmepr',[ind1 ind3 1],[xL yL tL]),3,ndays)./area(ind1:ind2,ind3:ind4)/rho0;
SS.dVdt = monmean(double(ncread(wnameI,'dVdt',[ind1 ind3 1],[xL yL tL])),3,ndays)*1e9/rho0./area(ind1:ind2,ind3:ind4);
SS.dHdt = monmean(double(ncread(wnameI,'dHdt',[ind1 ind3 1],[xL yL tL])),3,ndays)./area(ind1:ind2,ind3:ind4);
SS.ten = monmean(double(ncread(wnameI,'temp_tendency',[ind1 ind3 1],[xL yL tL])),3,ndays);
SS.diff_cbt = monmean(ncread(wnameI,'temp_vdiffuse_diff_cbt',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.KPPNL = monmean(ncread(wnameI,'temp_nonlocal_KPP',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.sbc = monmean(ncread(wnameI,'temp_vdiffuse_sbc',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.swh = monmean(ncread(wnameI,'sw_heat',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.frz = monmean(ncread(wnameI,'frazil',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.ets = monmean(ncread(wnameI,'temp_eta_smooth',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.pme = monmean(ncread(wnameI,'sfc_hflux_pme',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.rmx = monmean(ncread(wnameI,'temp_rivermix',[ind1 ind3 1],[xL yL tL]),3,ndays);
SS.adv = monmean(double(ncread(wnameI,'temp_advection',[ind1 ind3 1],[xL yL tL])),3,ndays);
SS.sub = monmean(double(ncread(wnameI,'temp_submeso',[ind1 ind3 1],[xL yL tL])),3,ndays);
if (haveGM)
    SS.gm = monmean(ncread(wnameI,'neutral_gm',[ind1 ind3 1],[xL yL tL]),3,ndays);
else
    SS.gm = zeros(size(SS.dVdt));
end 
if (haveRedi)
    SS.k33 = monmean(ncread(wnameI,'temp_vdiffuse_k33',[ind1 ind3 1],[xL yL tL]),3,ndays);
    SS.ndi = monmean(ncread(wnameI,'neutral_diffusion',[ind1 ind3 1],[xL yL tL]),3,ndays);
else
    SS.k33 = zeros(size(SS.dVdt));
    SS.ndi = zeros(size(SS.dVdt));
end
if (haveMDS)
    SS.mds = monmean(ncread(wnameI,'mixdownslope_temp',[ind1 ind3 1],[xL yL tL]),3,ndays);
else
    SS.mds = zeros(size(SS.dVdt));
end    
if (haveSIG)
    SS.sig = monmean(ncread(wnameI,'temp_sigma_diff',[ind1 ind3 1],[xL yL tL]),3,ndays);
else
    SS.sig = zeros(size(SS.dVdt));
end    
SS.dift = SS.diff_cbt+SS.KPPNL+SS.sbc+SS.swh+SS.frz+SS.ets+SS.pme+SS.rmx+SS.k33+SS.ndi+SS.mds+SS.sig;
    
txtrans = monmean(ncread(wnameI,'tx_trans',[ind1 ind3 1],[xL yL tL]),3,ndays)*tsc/rho0;
tytrans = monmean(ncread(wnameI,'ty_trans',[ind1 ind3 1],[xL yL tL]),3,ndays)*tsc/rho0;
txtrans(isnan(txtrans)) = 0;
tytrans(isnan(tytrans)) = 0;

qxtrans = monmean(ncread(wnameI,'temp_xflux_adv',[ind1 ind3 1],[xL yL tL]),3,ndays);
qytrans = monmean(ncread(wnameI,'temp_yflux_adv',[ind1 ind3 1],[xL yL tL]),3,ndays);

% Submesoscale and GM: two options:
if (haveSUB)
    if (doSGMviac)
        SS.dift = SS.dift + SS.sub;
    else
        qxtrans = qxtrans + monmean(ncread(wnameI,'temp_xflux_submeso',[ind1 ind3 1],[xL yL tL]),3,ndays);
        qytrans = qytrans + monmean(ncread(wnameI,'temp_yflux_submeso',[ind1 ind3 1],[xL yL tL]),3,ndays);
    end
end
if (haveGM)
    if (doSGMviac)
        SS.dift = SS.dift + SS.gm;
    else
        qxtrans = qxtrans + monmean(ncread(wnameI,'temp_xflux_gm',[ind1 ind3 1],[xL yL tL]),3,ndays);
        qytrans = qytrans + monmean(ncread(wnameI,'temp_yflux_gm',[ind1 ind3 1],[xL yL tL]),3,ndays);
    end
end
SS.JI = zeros(size(SS.JS));    
SS.QI = zeros(size(SS.JS));    
SS.JI(2:end,2:end) = (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                      +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(ind1+1:ind2,ind3+1:ind4);
SS.JI(1,2:end) = (txtrans(end,2:end) - txtrans(1,2:end) ...
                  +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(ind1,ind3+1:ind4);        
SS.QI(2:end,2:end) = (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                      +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(ind1+1:ind2,ind3+1:ind4);
SS.QI(1,2:end) = (qxtrans(end,2:end) - qxtrans(1,2:end) ...
                  +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(ind1,ind3+1:ind4);        

% Calculate other terms and residuals:
SS.SFC = SS.sbc + SS.swh + SS.frz + SS.ets + SS.pme + SS.rmx;
SS.VMX = SS.diff_cbt + SS.KPPNL + SS.k33;
SS.LMX = SS.ndi + SS.mds + SS.sig;

SS.HBresTEN = SS.ten - SS.adv - SS.sub - SS.gm - SS.SFC - SS.VMX - SS.LMX; % Residual from full Eulerian tendencies
SS.HBlcon = SS.adv + SS.gm + SS.sub - SS.QI;
SS.tenres = SS.ten - SS.dHdt;
SS.HBres = SS.dHdt - SS.SFC - SS.LMX - SS.QI;

SS.dHIdt = SS.dHdt - (SS.dVdt - SS.JI - SS.JS)*rho0*Cp*34.0;
SS.ndif = SS.dHIdt - SS.dift - SS.QI;

% Volume budget:
SS.VBres = SS.dVdt - SS.JI - SS.JS;

%% Global values:
names = fieldnames(SS);
tarea = sum(area(mask));
for i=1:length(names)
    eval(['total = sum(SS.' names{i} '(mask).*area(mask));']);
    disp(sprintf(['% 6.2e Wm-2 =  % 6.2e PW ' names{i}],total/tarea,total/1e15));
end

vnames = {'JS','dVdt','JI','VBres'};
for i=1:length(vnames)
    eval(['total = sum(SS.' vnames{i} '(mask).*area(mask));']);
    disp(sprintf(['% 6.2e Sv ' vnames{i}],total/1e6));
end

%%%% Plotting vertically-integrated Volume budget:
% $$$ xvec = 800:1:xL;
% $$$ yvec = 800:1:yL;
xvec = 1:1:xL;
yvec = 1:1:yL;
figure;
subplot(2,2,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.dVdt(xvec,yvec));
title('dVdt (m3s-1 / m2)');
caxis([-1 1]*1e-9);
subplot(2,2,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.JI(xvec,yvec));
title('JI');
caxis([-1 1]*1e-7);
subplot(2,2,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.JS(xvec,yvec));
title('JS');
caxis([-1 1]*1e-7);
subplot(2,2,4);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.VBres(xvec,yvec));
title('Residual');
caxis([-1 1]*1e-7);
colormap(redblue);


%%%% Compare to eta_t budget in ACCESS-OM2-1:
% Compare to eta_t tendency budget:
ename = [base 'ocean_btr.nc'];
SS.eTEN = ncread(ename,'eta_t_tendency');
% $$$ SS.deta_dt = ncread(ename,'deta_dt'); ! This is the same as eta_t_tendency
SS.eSMO = ncread(ename,'eta_smoother');
SS.eCON = ncread(ename,'conv_rho_ud_t')/rho0;
SS.eTENres = SS.eTEN(xvec,yvec)-SS.eSMO(xvec,yvec)-SS.eCON(xvec,yvec)-SS.JS(xvec,yvec);

figure;
subplot(2,2,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.eTEN(xvec,yvec));
title('eta tendency (m3s-1 / m2)');
caxis([-1 1]*1e-9);
subplot(2,2,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.eSMO(xvec,yvec));
title('eta smoother (m3s-1 / m2)');
caxis([-1 1]*1e-7);
subplot(2,2,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.eCON(xvec,yvec));
title('eta convergence (m3s-1 / m2)');
caxis([-1 1]*1e-7);
subplot(2,2,4);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.eTENres);
title('Residual (eta tend-smoo-con-JS) (m3s-1 / m2)');
caxis([-1 1]*1e-13);
colormap(redblue);

% snapshots eta_t and dzt calculation:
sname = [base 'ocean_snap.nc'];
sbname = [baser 'ocean_barotropic.res.nc'];
srname = [baser 'ocean_thickness.res.nc'];
eta1 = ncread(sbname,'eta_t');
eta2 = ncread(sname,'eta_t');
dzt1 = ncread(srname,'rho_dzt')/rho0;
dzt2 = ncread(sname,'dzt');
dzt1(isnan(dzt2)) = NaN;
SS.detadt = (eta2-eta1)/(ndays*86400);
SS.ddztdt = (nansum(dzt2,3)-nansum(dzt1,3))/(ndays*86400);

% Difference tendency calculations:
figure;
subplot(2,2,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.detadt(xvec,yvec)-SS.dVdt(xvec,yvec));
title('detadt - dVdt');
caxis([-1 1]*1e-11);
subplot(2,2,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.ddztdt(xvec,yvec) - SS.dVdt(xvec,yvec));
title('ddztdt - dVdt');
caxis([-1 1]*1e-11);
subplot(2,2,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.detadt(xvec,yvec) - SS.eTEN(xvec,yvec));
title('detadt - eta tendency (m3s-1 / m2)');
caxis([-1 1]*1e-9);

% Final comparisons:
figure;
subplot(2,2,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.dVdt(xvec,yvec)- SS.eTEN(xvec,yvec));
title('dVdt - eta tendency (m3s-1 / m2)');
caxis([-1 1]*1e-9);
subplot(2,2,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.JI(xvec,yvec)-SS.eCON(xvec,yvec));
title('JI - eta convergence (m3s-1 / m2)');
caxis([-1 1]*1e-9);
subplot(2,2,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.dVdt(xvec,yvec) - SS.JI(xvec,yvec) - SS.eSMO(xvec,yvec)-SS.JS(xvec,yvec));
title('dVdt - JI - JS - eta smoother (m3s-1 / m2)');
caxis([-1 1]*1e-9);
colormap(redblue);

vnames = {'JS','dVdt','JI','VBres','eTEN','eSMO','eCON','eTENres','detadt','ddztdt'};
for i=1:length(vnames)
    eval(['total = sum(SS.' vnames{i} '(mask).*area(mask));']);
    disp(sprintf(['% 6.2e Sv ' vnames{i}],total/1e6));
end


%%%% Check vertically-integrated Heat budget:
xvec = 1:4:xL;
yvec = 1:4:yL;
% $$$ xvec = 1:1:xL;
% $$$ yvec = 1:1:yL;
figure;
subplot(3,3,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.SFC(xvec,yvec));
title('Surface Fluxes');
caxis([-100 100]);
subplot(3,3,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.QI(xvec,yvec)+SS.LMX(xvec,yvec));
title('Lateral fluxes (QI, MDS, SIG, REDI)');
caxis([-100 100]);
subplot(3,3,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.dHdt(xvec,yvec));
title('dHdt');
caxis([-100 100]);
subplot(3,3,4);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.HBres(xvec,yvec));
title('dHdt - Surface - Lateral fluxes');
caxis([-5 5]);
subplot(3,3,5);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.HBresTEN(xvec,yvec));
title('Residual Eul. Tendency (Wm-2)');
caxis([-2e-4 2e-4]);
subplot(3,3,6);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.VMX(xvec,yvec));
title('Vertical Mixing terms');
caxis([-1e-5 1e-5]);
subplot(3,3,7);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.HBlcon(xvec,yvec));
title('Eul. adv+gm+sub  - QI (Wm-2)');
caxis([-1e-1 1e-1]);
subplot(3,3,8);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.tenres(xvec,yvec));
title('Eul. Tendency - dHdt');
caxis([-5 5]);
subplot(3,3,9);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),SS.tenres(xvec,yvec)+SS.HBres(xvec,yvec));
title('(Eul. Tendency - dHdt) + (dHdt - Surface - Lateral fluxes)');
caxis([-0.1 0.1]);
colormap(redblue);

% $$$ AOM2SS = SS;
% $$$ MOMSS = SS;
% $$$ AOM21SS = SS;
% $$$ 
% $$$ SS = AOM21SS;
% $$$ 
% $$$ %%%%%% Load ZA:
% $$$ region = 'Global';
% $$$ load([outD model sprintf('_output%03d',output) '_' region '_ZAHBud.mat']);
% $$$ 
% $$$ % ZA_I calculation:
% $$$ Tef = flipud(Te); % Eureka!!! But why?
% $$$ % $$$ Tef = Te;
% $$$ % $$$ ZA.AE = rho0*Cp*repmat(Tef',[yL 1]).*ZA.PSI;
% $$$ % $$$ ZA.AIadv = ZA.AHD -ZA.AE;
% $$$ % $$$ ZA.AI = ZA.AIadv+ZA.AHDGM+ZA.AHDSUB+ZA.AHDR;
% $$$ % $$$ ZA.JSH = ZA.JS.*repmat(Tef',[yL 1])*rho0*Cp;
% $$$ % $$$ ZA.PI  = ZA.P - ZA.JSH;
% $$$ % $$$ ZA.N = ZA.dHdt - ZA.dVdt.*repmat(Tef',[yL 1])*rho0*Cp;
% $$$ % $$$ 
% $$$ % $$$ ZA.AEdiff = diff(cat(1,zeros(1,TL+1),ZA.AE),[],1);
% $$$ % $$$ % $$$ dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA.AI-ZA.AHDR),[],1);
% $$$ % $$$ dAI_mR_dphi = ZA.AHDdiff-ZA.AEdiff+ZA.AHDGMdiff+ZA.AHDSUBdiff;
% $$$ 
% $$$ 
% $$$ % Take convergence:
% $$$ ZA.JI = -diff(cat(1,zeros(1,TL+1),ZA.PSI),[],1);
% $$$ ZA.QIadv = -diff(cat(1,zeros(1,TL+1),ZA.AHD),[],1);
% $$$ ZA.QIgm  = -diff(cat(1,zeros(1,TL+1),ZA.AHDGM),[],1);
% $$$ ZA.QIsub = -diff(cat(1,zeros(1,TL+1),ZA.AHDSUB),[],1);
% $$$ ZA.QI = ZA.QIadv+ZA.QIgm+ZA.QIsub;
% $$$ % $$$ ZA.AHDRdiff = diff(cat(1,zeros(1,TL+1),ZA.AHDR),[],1);
% $$$ 
% $$$ % Volume imbalance terms (should be zero integrated in T):
% $$$ ZA.VIMB = -(ZA.dVdt- ZA.JI-ZA.JS).*repmat(Tef',[yL 1])*rho0*Cp;
% $$$ 
% $$$ ZA.NNET = ZA.M+ZA.K33; % These have totals less than 10^5
% $$$ 
% $$$ ZA.I = ZA.dHdt + ZA.VIMB + ...
% $$$        -ZA.QI + ...
% $$$        -ZA.NNET-ZA.F-ZA.P-ZA.RED-ZA.MDS-ZA.SIG;
% $$$ % To be compared with:
% $$$ % $$$ SS.ndif = SS.dHdt - (SS.dVdt - SS.JI - SS.JS)*rho0*Cp*Te(Ti) - SS.dift - SS.QI;

% Check with zonally-averaged budget:
plot(yt,nansum(SS.VIVB.*area,1),'-k','linewidth',3);
hold on;
latVIVB = ZA.dVdt- ZA.JI-ZA.JS;
plot(yt,latVIVB(:,end),'-r','linewidth',3);


% Check with zonally-averaged budget:
plot(yt,nansum(SS.VIHB.*area,1),'-k','linewidth',3);
hold on;
latVIHB = ZA.dHdt - ZA.F - ZA.P - ZA.RED - ZA.MDS-ZA.SIG - ZA.QI;
plot(yt,latVIHB(:,end),'-r','linewidth',3);

