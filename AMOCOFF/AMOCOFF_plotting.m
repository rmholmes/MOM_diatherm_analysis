% Create anomalies and plot for AMOCOFF MOM025 simulation.
%
%

base = '/srv/ccrc/data03/z3500785/mom/MOM_HeatDiag_AMOCOFF/';

pre = 'v1p2p2_B1850_f19g16_conf9_BrR1000_';
mid = 'dtcount29.';
post = '.030101-040012';
ano = '1SvNAhosing_';

outpost = '_AMOCOFFanom.nc';

% XXX: Plot JFMA (seems to be North Atlantic deep mld seasons in
% MOM025) MLD and use a criteria (1000m?) to calculate an average
% dT, dQ and dSSS from Bryams runs in these regions. Then apply
% this over MOM025?

%%% Interpolated anomaly fields 
% MOM025 fields:
fname = [base 'kb3out120_SSTSSS.nc'];
lon = ncread(fname,'geolon_t');
lat = ncread(fname,'geolat_t');
SST = ncread(fname,'temp'); % Annual mean kb3seg SST 
SSS = ncread(fname,'salt'); % Annual mean kb3seg SSS
% ncrcat -v temp,salt -d st_ocean,0,0 ../../archive/MOM_HeatDiag_kb3seg/output120/ocean.nc kb3out120_SSTSSS.nc
sname = [base 'salt_sfc_restore_avg.nc']; % Annual mean SSS
                                          % restoring field
% ncra -v SALT salt_sfc_restore.nc salt_sfc_restore_avg.nc
snameout = [base 'salt_sfc_restore' outpost]; % out file for SSS anomaly.
copyfile(sname,snameout);
SSSrst = ncread(sname,'SALT'); % SSS restoring field

tname = [base 't_10_avg.nc']; % Tair
lonA = ncread(tname,'LON');
latA = ncread(tname,'LAT');
[LONA,LATA] = ndgrid(lonA,latA);
% ncra -v T_10_MOD t_10.nc t_10_avg.nc
tnameout = [base 't_10' outpost]; % out file for SSS anomaly.
copyfile(tname,tnameout);
T10 = ncread(tname,'T_10_MOD'); % Tair base

qname = [base 'q_10_avg.nc']; % Qair
% ncra -v Q_10_MOD q_10.nc q_10_avg.nc
qnameout = [base 'q_10' outpost]; % out file for SSS anomaly.
copyfile(qname,qnameout);
Q10 = ncread(qname,'Q_10_MOD'); % Tair base

% Atmospheric anomaly fields:

% Tair and Qair:
fname = [base pre mid 'TREFHT' post '.ncra50.nc'];
lonAAA = ncread(fname,'lon');
latAAA = ncread(fname,'lat');
[LONAAA,LATAAA] = ndgrid(lonAAA,latAAA);
TAIRano = ncread([base pre ano mid 'TREFHT' post '.ncra50.diff.nc'],'TREFHT');
QAIRano = ncread([base pre ano mid 'QREFHT' post '.diff.ncra50.nc'],'QREFHT');

% Interpolate:
X = cat(1,lonAAA(end)-360,lonAAA,lonAAA(1)+360);
Z = cat(1,TAIRano(1,:),TAIRano,TAIRano(end,:));
TAIR = interp2(X,latAAA,Z',LONA,LATA,'linear');
Z = cat(1,QAIRano(1,:),QAIRano,QAIRano(end,:));
QAIR = interp2(X,latAAA,Z',LONA,LATA,'linear');

% MASK:
mask = (LONA>250 & LONA<360 | LONA < 30) & LATA>40 & LATA<70;
buf = 5;
mask = filter_field(mask,buf,'-s');
% fix end points:
mask(1:(buf-1)/2,:) = repmat(mask((buf-1)/2+1,:),[(buf-1)/2 1]);
mask(end-(buf-1)/2+1:end,:) = repmat(mask(end-(buf-1)/2,:),[(buf-1)/2 1]);

ncid = netcdf.open(tnameout,'NC_WRITE');
netcdf.putVar(ncid,netcdf.inqVarID(ncid,'T_10_MOD'),TAIR.*mask);
netcdf.close(ncid);
ncid = netcdf.open(qnameout,'NC_WRITE');
netcdf.putVar(ncid,netcdf.inqVarID(ncid,'Q_10_MOD'),QAIR.*mask);
netcdf.close(ncid);

% SSS:
fname = [base pre mid 'sst' post '.ncra50.nc'];
lonSSS = ncread(fname,'TLON');
latSSS = ncread(fname,'TLAT');
SSSano = ncread([base pre ano mid 'SALT' post '.ncra50.diff.nc'],'SALT');

% Interpolate:
xvec = [lonSSS(:)-360; lonSSS(:); lonSSS(:)+360];yvec = [latSSS(:); ...
                    latSSS(:); latSSS(:)];zvec = [SSSano(:); SSSano(:); SSSano(:)];
NaNs = isnan(xvec) | isnan(zvec);
xvec = xvec(~NaNs);yvec = yvec(~NaNs);zvec = zvec(~NaNs);
xvecO = lon(:);
xvecO(xvecO<0) = xvecO(xvecO<0)+360;
yvecO = lat(:);
vq = griddata(xvec(:),yvec(:),zvec(:),xvecO(:),yvecO(:));
SSS = reshape(vq,[length(lon(:,1)) length(lon(1,:))]);
SSS(isnan(SSS)) = 0;

% Interpolate mask:
xvec = [LONA(:)-360; LONA(:); LONA(:)+360];yvec = [LATA(:); LATA(:); ...
                    LATA(:)];zvec = [mask(:); mask(:); mask(:)];
vq = griddata(xvec(:),yvec(:),zvec(:),xvecO(:),yvecO(:));
maskHR = reshape(vq,[length(lon(:,1)) length(lon(1,:))]);
maskHR(isnan(maskHR)) = 0;

ncid = netcdf.open(snameout,'NC_WRITE');
netcdf.putVar(ncid,netcdf.inqVarID(ncid,'SALT'),SSS.*maskHR);
netcdf.close(ncid);

% Plotting:
% $$$ figure;
clf;

load('coasts.mat');
lonc = lon_coast;
latc = lat_coast;
lonc(lonc<(80-360)) = lonc(lonc<(80-360))+360;
lonc(abs(lonc)<0.1) = NaN;

subplot(2,2,1);
X = LONA;Y = LATA;Z = TAIR.*mask;
[tmp ind] = min(abs(LONA(:,1)-80));
X = cat(1,X((ind+1):end,:)-360,X(1:ind,:));
Y = cat(1,Y((ind+1):end,:),Y(1:ind,:));
Z = cat(1,Z((ind+1):end,:),Z(1:ind,:));
pcolPlot(X,Y,Z);
caxis([-12 12]);
title('Tair Anomaly ($^\circ$C)');
hold on;
plot(lonc,latc,'.k','markersize',1);
xlim([-120 100]);
ylim([20 80]);

subplot(2,2,2);
X = LONA;Y = LATA;Z = QAIR.*mask;
[tmp ind] = min(abs(LONA(:,1)-80));
X = cat(1,X((ind+1):end,:)-360,X(1:ind,:));
Y = cat(1,Y((ind+1):end,:),Y(1:ind,:));
Z = cat(1,Z((ind+1):end,:),Z(1:ind,:));
pcolPlot(X,Y,Z);
caxis([-0.004 0.004]);
hold on;
plot(lonc,latc,'.k','markersize',1);
title('Qair Anomaly (g/kg)');
xlim([-120 100]);
ylim([20 80]);

subplot(2,2,3);
X = lon;Y = lat;Z = SSS.*maskHR;
% $$$ [tmp ind] = min(abs(lon(:,1)-100));
% $$$ X = cat(1,X((ind+1):end,:)-360,X(1:ind,:));
% $$$ Y = cat(1,Y((ind+1):end,:),Y(1:ind,:));
% $$$ Z = cat(1,Z((ind+1):end,:),Z(1:ind,:));
pcolPlot(X,Y,Z);
caxis([-15 15]);
colormap(redblue);
hold on;
plot(lonc,latc,'.k','markersize',1);
title('SSS Anomaly (psu)');
xlim([-120 100]);
ylim([20 80]);






























%%%%%%%%%%%%%%%%%% PLOTTING 

% Plot spatial pattern of last 50 year mean and anomalies:

% MOM025 KB3seg fields:
fname = [base 'kb3out120_SSTSSS.nc'];
lon = ncread(fname,'geolon_t');
lat = ncread(fname,'geolat_t');
SST = ncread(fname,'temp');
SSS = ncread(fname,'salt');
SSSrst = ncread([base 'salt_sfc_restore_avg.nc'],'SALT');

xlims = [min(min(lon)) max(max(lon))];
ylims = [-90 90];
xlims = [250 360]-360;
ylims = [0 80];

figure;
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,3,1);
pcolPlot(lon,lat,SST);
caxis([-2 25]);
colorbar;xlim(xlims);ylim(ylims);
title('MOM025 SST ($^\circ$C)');
set(gca,'color','k');
subplot(2,3,2);
pcolPlot(lon,lat,SSS);
caxis([28 38]);
colorbar;xlim(xlims);ylim(ylims);
title('MOM025 SSS (psu)');
set(gca,'color','k');
subplot(2,3,3);
pcolPlot(lon,lat,SSSrst);
caxis([28 38]);
colorbar;xlim(xlims);ylim(ylims);
title('MOM025 restoring SSS (psu)');
set(gca,'color','k');
subplot(2,3,4);
pcolPlot(lon,lat,SSS-SSSrst);
caxis([-2 2]);
colorbar;xlim(xlims);ylim(ylims);
title('MOM025 SSS - restoring SSS (psu)');
set(gca,'color','k');
colormap(redblue);


% SST and SSS from Bryam's runs:
fname = [base pre mid 'sst' post '.ncra50.nc'];
lonO = ncread(fname,'TLON');
latO = ncread(fname,'TLAT');
SSTcon = ncread(fname,'sst');
SSTper = ncread([base pre ano mid 'sst' post '.ncra50.nc'],'sst');
SSTano = ncread([base pre ano mid 'sst' post '.ncra50.diff.nc'],'sst');
SSScon = ncread([base pre mid 'SALT' post '.ncra50.nc'],'SALT');
SSSper = ncread([base pre ano mid 'SALT' post '.ncra50.nc'],'SALT');
SSSano = ncread([base pre ano mid 'SALT' post '.ncra50.diff.nc'],'SALT');
[xL,yL] = size(lonO);

NaNs = latO>50 & lonO<5;
lonO(NaNs) = NaN;
latO(NaNs) = NaN;

% Shift longitude axis:
ind = 36;
vars = {'lonO','latO','SSTcon','SSTano','SSScon','SSSano','SSTper','SSSper'};
for fi=1:length(vars)
    eval([vars{fi} ' = cat(1,' vars{fi} '(ind+1:end,:),' vars{fi} ...
          '(1:ind,:));']);
end

xlims = [0 360];
ylims = [-90 90];
xlims = [250 360];
ylims = [0 80];

figure;
set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'defaulttextfontSize',10);
% $$$ set(gcf,'defaultaxesfontSize',10);
subplot(2,3,1);
pcolPlot(lonO,latO,SSTcon);
caxis([-2 25]);
colorbar;xlim(xlims);ylim(ylims);
title('SST Control ($^\circ$C)');
set(gca,'color','k');
subplot(2,3,2);
pcolPlot(lonO,latO,SSTper);
caxis([-2 25]);
colorbar;xlim(xlims);ylim(ylims);
title('SST Perturbed ($^\circ$C)');
set(gca,'color','k');
subplot(2,3,3);
pcolPlot(lonO,latO,SSTano);
caxis([-10 10]);
colorbar;xlim(xlims);ylim(ylims);
title('SST Anomaly ($^\circ$C)');
set(gca,'color','k');
subplot(2,3,4);
pcolPlot(lonO,latO,SSScon);
caxis([28 38]);
colorbar;xlim(xlims);ylim(ylims);
title('SSS Control (psu)');
set(gca,'color','k');
subplot(2,3,5);
pcolPlot(lonO,latO,SSSper);
caxis([28 38]);
colorbar;xlim(xlims);ylim(ylims);
title('SSS Perturbed (psu)');
set(gca,'color','k');
subplot(2,3,6);
pcolPlot(lonO,latO,SSSano);
caxis([-10 10]);
colorbar;xlim(xlims);ylim(ylims);
title('SSS Anomaly (psu)');
set(gca,'color','k');
colormap(redblue);
plot([300 359 359 300 300],[40 40 65 65 40],'-m','linewidth',2);

% Tair and Qair:
fname = [base pre mid 'TREFHT' post '.ncra50.nc'];
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
[lonO,latO] = ndgrid(lon,lat);
TAIRcon = ncread(fname,'TREFHT')-273.15;
TAIRper = ncread([base pre ano mid 'TREFHT' post '.ncra50.nc'],'TREFHT')-273.15;
TAIRano = ncread([base pre ano mid 'TREFHT' post '.ncra50.diff.nc'],'TREFHT');
QAIRcon = ncread([base pre mid 'QREFHT' post '.ncra50.nc'],'QREFHT')*1000;
QAIRper = ncread([base pre ano mid 'QREFHT' post '.ncra50.nc'],'QREFHT')*1000;
QAIRano = ncread([base pre ano mid 'QREFHT' post '.diff.ncra50.nc'],'QREFHT')*1000;
[xL,yL] = size(lonO);

xlims = [0 360];
ylims = [-90 90];
xlims = [250 360];
ylims = [0 80];
load('coasts.mat');


figure;
set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'defaulttextfontSize',10);
set(gcf,'defaultaxesfontSize',10);
subplot(2,3,1);
pcolPlot(lonO,latO,TAIRcon);
caxis([-25 25]);
colorbar;xlim(xlims);ylim(ylims);
title('TAIR Control ($^\circ$C)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
subplot(2,3,2);
pcolPlot(lonO,latO,TAIRper);
caxis([-25 25]);
colorbar;xlim(xlims);ylim(ylims);
title('TAIR Perturbed ($^\circ$C)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
subplot(2,3,3);
pcolPlot(lonO,latO,TAIRano);
caxis([-10 10]);
colorbar;xlim(xlims);ylim(ylims);
title('TAIR Anomaly ($^\circ$C)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
subplot(2,3,4);
pcolPlot(lonO,latO,QAIRcon);
caxis([0 20]);
colorbar;xlim(xlims);ylim(ylims);
title('QAIR Control (g/kg)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
subplot(2,3,5);
pcolPlot(lonO,latO,QAIRper);
caxis([0 20]);
colorbar;xlim(xlims);ylim(ylims);
title('QAIR Perturbed (g/kg)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
subplot(2,3,6);
pcolPlot(lonO,latO,QAIRano);
caxis([-5 5]);
colorbar;xlim(xlims);ylim(ylims);
title('QAIR Anomaly (g/kg)');
set(gca,'color','k');hold on;plot(lon_coast(:)+360,lat_coast(:),'-k','linewidth',0.5);
colormap(redblue);


% Time series:
fname = [base pre ano mid 'sst' post '.diff.nc'];
lonO = ncread(fname,'TLON');
latO = ncread(fname,'TLAT');
time = ncread(fname,'time');
SSTano = ncread(fname,'sst');
SSSano = ncread([base pre ano mid 'SALT' post '.diff.nc'],'SALT');
fname = [base pre ano mid 'TREFHT' post '.diff.nc'];
lon = ncread(fname,'lon');
lat = ncread(fname,'lat');
[lonA,latA] = ndgrid(lon,lat);
TAIRano = ncread(fname,'TREFHT');
[xL,yL,tL] = size(SSTano);

reg = lonO > 300 & lonO < 359 & latO > 40 & latO < 65;
regA = lonA > 300 & lonA < 359 & latA > 40 & latA < 65;
anomSST = zeros(tL,1);
anomSSS = zeros(tL,1);
anomTAIR = zeros(tL,1);
for i=1:tL
    i
    tmp = SSTano(:,:,i);
    anomSST(i) = nanmean(nanmean(tmp(reg),1),2);
    tmp = SSSano(:,:,i);
    anomSSS(i) = nanmean(nanmean(tmp(reg),1),2);
    tmp = TAIRano(:,:,i);
    anomTAIR(i) = nanmean(nanmean(tmp(regA),1),2);
end

subplot(2,1,1);
plot(time/365,anomSST,'-k');
hold on;
plot(time/365,anomTAIR,'-r');
plot(time/365,filter_field(anomSST,13,'-t'),'-k','linewidth',2);
plot(time/365,filter_field(anomTAIR,13,'-t'),'-r','linewidth',2);
legend('SST','Tair');
title('North Atlantic SST anomaly');
xlabel('Year');
ylabel('$^\circ$C');
grid on;
ylim([-15 2]);
xlim([301 401]);
subplot(2,1,2);
plot(time/365,anomSST,'-k');
hold on;
plot(time/365,filter_field(anomSST,13,'-t'),'-k','linewidth',2);
title('North Atlantic SSS anomaly');
xlabel('Year');
ylabel('psu');
grid on;
ylim([-8 0]);
xlim([301 401]);
