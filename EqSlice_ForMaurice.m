%% This script extracts an equatorial slice in the Pacific of data
%% from MOM025 runs

model = 'MOM025';
baseD = '/short/e14/rmh561/mom/archive/MOM_HeatDiag_kb3seg/'; %Data Directory.

for output=75:79

    base = [baseD sprintf('output%03d/',output)];
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
    m3name = [base 'ocean.nc'];
else
    fname = [base 'ocean.nc'];
end
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
         
lon = ncread(hname,'geolon_t');lat = ncread(hname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonu = ncread(hname,'geolon_c');latu = ncread(hname,'geolat_c');

z = ncread(hname,'st_ocean');zL = length(z);

time = ncread(hname,'time');
tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

%% Get equatorial slices of variables:
latsl = 0;
[tmp ltind] = min(abs(lat(1,:)-latsl));
[tmp ln1] = min(abs(lon(:,ltind)+240));
[tmp ln2] = min(abs(lon(:,ltind)+70));


temp = squeeze(ncread(fname,'temp',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
u = squeeze(ncread(fname,'u',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
v = squeeze(ncread(fname,'v',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
kappa = squeeze(ncread(fname,'diff_cbt_t',[ln1 ltind 1 1],[ln2-ln1+1 1 zL tL]));
taux = squeeze(ncread(fname,'tau_x',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
tauy = squeeze(ncread(fname,'tau_y',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
mld = squeeze(ncread(fname,'mld',[ln1 ltind 1],[ln2-ln1+1 1 tL]));
vdif = squeeze(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
vnlc = squeeze(ncread(wname,'temp_nonlocal_KPP_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
pmer = squeeze(ncread(wname,'sfc_hflux_pme_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_rivermix_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
sufc = squeeze(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'frazil_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'temp_eta_smooth_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]) + ...
               ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));
swrd = squeeze(ncread(wname,'sw_heat_on_nrho',[ln1 ltind 1 1],[ln2-ln1+1 1 TL tL]));

[Xt,Zt] = ndgrid(lon(ln1:ln2,ltind),z);
[Xu,Zu] = ndgrid(lonu(ln1:ln2,ltind),z);

name = [outD 'mat_data/' model sprintf('_output%03d',output) '_varsat_Eq.mat']
save(name,'Xt','Zt','Xu','Zu','temp','u','v','kappa','taux','tauy','mld', ...
     'vdif','vnlc','pmer','sufc','swrd');

end


%%% Plot Equatorial Slices:

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
outputs = [75:79];

ndays = [31 28 31 30 31 30 31 31 30 31 30 31];

% Load Variable and calculate mean:
load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(1))]);
vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd'};
for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
end
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(i))]);
    for i=1:length(vars)
        eval([vars{i} 'a = ' vars{i} 'a + ' vars{i} ';']);
    end
end
for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a/length(outputs);']);
    eval(['clear ' vars{i} 'a;']);
end
[xL,zL,tL] = size(temp);
TL = length(T);

% Depth of isotherms:
Zi = zeros(xL,TL,tL);
for ti=1:tL
    for xi=1:xL
        tvec = squeeze(temp(xi,:,ti));
        zvec = -Zt(xi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(xi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(xi,:,ti)),1,'last');
        Zi(xi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Xi = repmat(Xt(:,1),[1 TL]);

% $$$ var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux
% $$$ clim = [-400 0];
% $$$ sp = 10;
% $$$ doWMT = 0;

var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
% $$$ var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
clim = [-1 1]*1e-5*86400; % FOR WMT
sp = 0.1*1e-5*86400;
doWMT = 1;

months = {[1:12],[3],[7],[11]};
monthsu01 = {[1:4],[1],[3],[4]};
labels = {'Annual','March','July','November'};

%Colormap:
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

if (doWMT)
    cmap = redblue(npts-3);
else
    cmap = parula(npts-3);
    cmap(end,:) = [1 1 1];
    cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
end

figure;
set(gcf,'Position',[1          36        1920         970]);
for i=1:length(months)
subplot(2,2,i);
contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
[c,h] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
clabel(c,h,[0:2:35]);
[c,h] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[23 23],'-k','linewidth',2);
if (strcmp(model,'MOM01'))
    mnu = monthsu01{i};
else
    mnu = months{i};
end
ucol = [0.8706    0.4902         0];
[c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
                'color',ucol);
[c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
                'color',ucol);
plot(Xu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
% $$$ clabel(c,h,'color','w');
ylim([-300 0]);
xlim([-220 -80]);
cb = colorbar;
if (doWMT)
    ylabel(cb,'m/day');
else
    ylabel(cb,'Wm$^{-2}$');
end
xlabel('Longitude ($^\circ$E)');
ylabel('Depth (m)');
caxis(clim);
text(-218,-288,labels{i},'Backgroundcolor','w','FontSize',20);

LabelAxes(gca,i,20,0.008,0.95);
end
colormap(cmap);
