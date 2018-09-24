
%This script compares the TPBAS4 seasonal cycle bulk simulations to
%observations. 

addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));
startup;

% OSU details:
OSUbase = '/g/data/e14/rmh561/OSU/';

% MOM025 details:
MOMname = '/short/e14/rmh561/mom/archive/MOM_HeatDiag_kb3seg/output080/ocean.nc';
MOMnameH = '/short/e14/rmh561/mom/archive/MOM_HeatDiag_kb3seg/output080/ocean_heat.nc';
MOMlon = ncread(MOMname,'xt_ocean');
MOMlat = ncread(MOMname,'yt_ocean');
MOMdep = ncread(MOMname,'st_ocean');

%% 140W, 0N TAO comparison:

% Get MOM data:
[tmp lnind] = min(abs(MOMlon+140));
[tmp ltind] = min(abs(MOMlat));
MOMU = squeeze(ncread(MOMname,'u',[lnind ltind 1 1],[1 1 50 12]));
MOMT = squeeze(ncread(MOMname,'temp',[lnind ltind 1 1],[1 1 50 12]));
MOMF = squeeze(ncread(MOMnameH,'temp_vdiffuse_diff_cbt',[lnind ltind 1 1],[1 1 50 12]));
MOMFkppiw = squeeze(ncread(MOMnameH,'temp_vdiffuse_diff_cbt_kppiw',[lnind ltind 1 1],[1 1 50 12]));
MOMFkppish = squeeze(ncread(MOMnameH,'temp_vdiffuse_diff_cbt_kppish',[lnind ltind 1 1],[1 1 50 12]));
MOMFkppbl = squeeze(ncread(MOMnameH,'temp_vdiffuse_diff_cbt_kppbl',[lnind ltind 1 1],[1 1 50 12]));
MOMFkppNL = squeeze(ncread(MOMnameH,'temp_nonlocal_KPP',[lnind ltind 1 1],[1 1 50 12]));
MOMFswh = squeeze(ncread(MOMnameH,'sw_heat',[lnind ltind 1 1],[1 1 50 12]));
MOMFsbc = squeeze(ncread(MOMnameH,'temp_vdiffuse_sbc',[lnind ltind 1 1],[1 1 50 12]));
intvars = {'MOMF','MOMFkppiw','MOMFkppish','MOMFkppbl','MOMFswh','MOMFkppNL', ...
           'MOMFsbc'};
for i = 1:length(intvars)
    eval([intvars{i} '(isnan(' intvars{i} ')) = 0;']);
    eval([intvars{i} ' = -cumsum(' intvars{i} ',1,''reverse'');']);
end

% Load OSU data:
OSUname = [OSUbase 'vel_0_140W_daily.mat'];
load(OSUname);
OSUname = [OSUbase 'T_0_140W_daily.mat'];
load(OSUname);
OSUname = [OSUbase 'chipods_0_140W_alldepths_daily.mat'];
load(OSUname);

Utao = veldy.u;
dUdztao = veldy.dudz;
Ttao = tdy.T;
Ztao = veldy.depth;
TAOnum = veldy.time;
TAOvec = datevec(TAOnum);
CHIJq = cdy.Jq;
CHIJqZ = cdy.depth;
CHIvec = datevec(cdy.time);

% Set time period:
TAOinds = find(TAOvec(:,2) >= 9);
MOMinds = (1:12) >=9;
CHIinds = find(CHIvec(:,2) >= 9);

figure;
subplot(1,3,1);
plot(nanmean(Utao(:,TAOinds),2),-Ztao,'-k','linewidth',2);
hold on;
plot(nanmean(MOMU(:,MOMinds),2),-MOMdep,'--k','linewidth',2);
ylim([-250 0]);
xlabel('U (ms$^{-1}$)');ylabel('Depth (m)');legend('TAO','MOM025 kb3seg');
xlim([-0.5 1.5]);grid on;
subplot(1,3,2);
plot(nanmean(Ttao(:,TAOinds),2),-Ztao,'-k','linewidth',2);
hold on;
plot(nanmean(MOMT(:,MOMinds),2),-MOMdep,'--k','linewidth',2);
ylim([-250 0]);
xlabel('T ($^\circ$C)');ylabel('Depth (m)');title('Annual');
xlim([10 27]);grid on;
subplot(1,3,3);
plot(nanmean(CHIJq(:,CHIinds),2),-CHIJqZ,'-k','linewidth',2);
hold on;
plot(nanmean(MOMF(:,MOMinds),2),-MOMdep,'--k','linewidth',2);
plot(nanmean(MOMFkppiw(:,MOMinds),2),-MOMdep,'--r','linewidth',2);
plot(nanmean(MOMFkppish(:,MOMinds),2),-MOMdep,'--b','linewidth',2);
plot(nanmean(MOMF(:,MOMinds)-MOMFkppish(:,MOMinds)-MOMFkppiw(:,MOMinds),2),-MOMdep,'--','color',[0 0.5 0],'linewidth',2);
% $$$ plot(nanmean(MOMFkppbl(:,MOMinds),2),-MOMdep,'--','color',[0 0.5 0],'linewidth',2);
plot(nanmean(MOMFkppNL(:,MOMinds),2),-MOMdep,':','color',[0 0.5 0],'linewidth',2);
plot(nanmean(MOMFswh(:,MOMinds),2)+nanmean(MOMFsbc(:,MOMinds),2),-MOMdep,'--m','linewidth',2);
legend('Chipods','MOM Total','MOM background','MOM Shear','MOM BL','MOM Non-local','MOM Surface Forcing');
ylim([-250 0]);
xlabel('$J_q$ (Wm$^{-2}$)');ylabel('Depth (m)');
xlim([-200 0]);grid on;
