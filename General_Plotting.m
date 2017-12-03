% This script analyses the heat budget in MOM025-CORENYF

%% file-names and grid properties:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/output007/';
baser = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/restart006/';
hname = [base 'ocean_heat.nc'];
fname = [base 'ocean.nc'];
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
rnameT = [baser 'ocean_temp_salt.res.nc'];
rnameZ = [baser 'ocean_thickness.res.nc'];
         
lon = ncread(hname,'geolon_t');
lat = ncread(hname,'geolat_t');
area = ncread(gname,'area_t');
[xL,yL] = size(lon);

z = ncread(hname,'st_ocean');
zL = length(z);

time = mod(ncread(hname,'time'),365);
time_snap = [0; mod(ncread(sname,'time'),365)];
if (time_snap(end) == 0) time_snap(end) = 365;end
tL = length(time);

TL = 161;
T = linspace(-3,35,TL);
dT = T(2)-T(1);
Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

% $$$ xvec = 1:3:xL;
% $$$ yvec = 1:3:yL;
% $$$ xvec2 = 1:20:xL;
% $$$ yvec2 = 1:20:yL;
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%%% Plot Mean Heat Flux and SST:
xvec = 1:1:xL;
yvec = 1:1:yL;

shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL 12]);
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 12]));

months = {[1:12], ...
          [3], ...
          [7], ...
          [11]};
monconv = datevec(mod(time,365));
monconv = monconv(:,2);
tmp = monconv;
for i=1:length(monconv)
    [pp tmp(i)] = min(abs(monconv-i));
end
monconv = tmp;

labels = {'Annual Average', ...
          'March', ...
          'July', ...
          'November'};

figure;
set(gcf,'Position',[1          36        1920         970]);
%set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
poss = [0.1300    0.4553    0.7402    0.4697; ...
        0.1300    0.1389    0.2343    0.2680; ...
        0.3951    0.1389    0.2343    0.2680; ...
        0.6681    0.1389    0.2343    0.2680];
for i=1:length(months)
    if (i == 1)
        subplot(5,3,[1 9]);
    else
        subplot(5,3,[10 13]+(i-2));
    end
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    Z = mean(shflux(:,:,monconv(months{i})),3);
    Z2 = mean(SST(:,:,monconv(months{i})),3);
    Z = Z(xvec,yvec);
    Z2 = Z2(xvec,yvec);
    contourf(X,Y,Z,[-1e10 -500:20:500 1e10],'linestyle','none');
    hold on;
% $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$     if (i==1)
        [c,h] = contour(X,Y,Z2,[-3:2:35],'-k');
% $$$     else
% $$$         [c,h] = contour(X,Y,Z2,[-3:4:35],'-k');
% $$$     end
    clabel(c,h);
    caxis([-200 200]);
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
    end
    ylim([-75 60]);
    if (i>1)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i<=2)
        ylabel('Latitude ($^\circ$N)');
    end
    if (i>1)
        text(-276,53,labels{i},'BackgroundColor','w');
    else
        text(-278,55,labels{i},'BackgroundColor','w');
    end        
    set(gca,'Position',[poss(i,:)]);
    set(gca,'color','k');
end 
colormap(redblue);

%%% Plot Seasonal cycle of winds and SST:
xvec2 = 1:20:xL;
yvec2 = 1:20:yL;
figure;
set(gcf,'Position',get(0,'ScreenSize'));
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',20);
for i=1:6
    month = 2*i-1;
    subplot(3,2,i);
    SST = ncread(fname,'temp',[1 1 1 month],[xL yL 1 1]);
    tau_x = ncread(fname,'tau_x',[1 1 month],[xL yL 1]);
    tau_y = ncread(fname,'tau_y',[1 1 month],[xL yL 1]);
    contourf(lon(xvec,yvec),lat(xvec,yvec),SST(xvec,yvec),[-5:1:40],'linestyle','none');
    hold on;
    quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$     contour(lon(xvec,yvec),lat(xvec,yvec),sqrt(tau_x(xvec,yvec).^2+tau_y(xvec,yvec).^2),[0.15 ...
% $$$                    0.2 0.25 0.3],'-k');
    cb = colorbar;
    caxis([0 32]);
    ylim([-75 65]);
    xlim([-280 80]);
    set(gca,'color','k');
    xlabel('Longitude ($^\circ$E)');
    ylabel('Latitude ($^\circ$N)');
    title(['MOM025 ' txtmonth{month} ' SST ($^\circ$C) and wind stress, ']);
end
colormap(redblue);


%% Plot Depth averaged T:
load('MOM025_WMHeatBudget.mat');

[X,Y] = ndgrid(0.5:1:11.5,-z);
figure;
contourf(X,Y,Temp'-repmat(mean(Temp',1),[12 1]),[-2:0.025:2],'linestyle','none');
ylim([-150 0]);
ylabel('Depth (m)');
xlabel('Month');
title('MOMO25 Global Temperature Anomaly ($^\circ$C)');
cb = colorbar;
caxis([-0.75 0.75]);

%% Plot zonal averaged T:
TZ = zeros(yL,zL);
[X,Y] = ndgrid(lat(1,:),-z);
for zi=1:zL
    zi
    TZ(:,zi) = squeeze(nanmean(nanmean(ncread(fname,'temp',[1 1 zi 1],[xL ...
                        yL 1 12]),4),1));
end

[c,h] = contourf(X,Y,TZ,[-2:1:35],'-k');
clabel(c,h);

xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
title('MOM025 Annual Average Temperature ($^\circ$C)');
cb = colorbar;
colormap(redblue);
ylim([-1000 0]);
caxis([5 28]);



