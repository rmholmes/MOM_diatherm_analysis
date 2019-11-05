% Quick calculation of WWV diabatic terms

clear all;

%base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
base = '/g/data/e14/rmh561/access-om2/archive/025deg_jra55_iaf_Maurice/mat_data/';
model = 'ACCESS-OM2_025deg_jra55_iaf';
outputs = [39];

% mask:
load([base 'wwv_mask.mat']);

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
    ndays = ndays(1:12);
end
if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
region = 'Global';
nyrs = tL/12;szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
yrs = 1:nyrs;
months = 1:12;
ycur = 1;

% Get JS, JA and dWWVdt:
wname = '/g/data/e14/rmh561/access-om2/archive/025deg_jra55_iaf_Maurice/output039/ocean/ocean_wmass.nc';

Tl = 20;
[tmp ind] = min(abs(T-Tl));
dWWVdt = squeeze(nansum(ncread(wname,'dVdt',[1 1 ind 1],[xL yL TL-ind+1 tL])/rho0*1e9,3));
dWWVdt = double(squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo,[1 1 tL]).*dWWVdt,1),2)));
JS = squeeze(nansum(ncread(wname,'mass_pmepr_on_nrho',[1 1 ind 1],[xL yL TL-ind+1 tL])/rho0,3));
JS = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo,[1 1 tL]).*JS,1),2));
tytrans = squeeze(nansum(ncread(wname,'ty_trans_nrho',[1 1 ind 1],[xL yL TL-ind+1 tL])/rho0,3));
txtrans = squeeze(nansum(ncread(wname,'tx_trans_nrho',[1 1 ind 1],[xL yL TL-ind+1 tL])/rho0,3));

mask_S = zeros(size(wwv_mask_2S_borneo));
mask_N = mask_S;mask_W = mask_S; mask_E = mask_S;
for i=1:xL
    ind = find(wwv_mask_2S_borneo(i,:),1,'first');
    if (length(ind)>0)
        mask_S(i,ind-1) = 1;
    end
    ind = find(wwv_mask_2S_borneo(i,:),1,'last');
    if (length(ind)>0)
        mask_N(i,ind) = 1;
    end
end
for i=1:yL
    ind = find(wwv_mask_2S_borneo(:,i),1,'first');
    if (length(ind)>0)
        mask_W(ind-1,i) = 1;
    end
    ind = find(wwv_mask_2S_borneo(:,i),1,'last');
    if (length(ind)>0)
        mask_E(ind-1,i) = 1;
    end
end

JA_S = squeeze(nansum(nansum(repmat(mask_S,[1 1 tL]).*tytrans,1),2));
JA_N = squeeze(-nansum(nansum(repmat(mask_N,[1 1 tL]).*tytrans,1),2));
JA_W = squeeze(nansum(nansum(repmat(mask_W,[1 1 tL]).*txtrans,1),2));
JA_E = squeeze(-nansum(nansum(repmat(mask_E,[1 1 tL]).*txtrans,1),2));

% Get WMT terms:
Tl = 19.75;
TYPE = 'WMT';
load([base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']);

G_M = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*WMTM,1),2));
G_I = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*WMTI,1),2));
G_R = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*(WMTK+WMTR),1),2));
G_F = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*WMTF,1),2));
G_P = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*WMTP,1),2));
G_SP = squeeze(nansum(nansum(repmat(wwv_mask_2S_borneo.*area,[1 1 tL]).*WMTSP,1),2));

save([base model '_Ryans_WWV_WMTterms_output039.mat'],'G_M','G_I','G_R','G_F','G_P','G_SP','JS','dWWVdt','JA_W','JA_E','JA_S','JA_N');

load([base model '_Ryans_WWV_WMTterms_output039.mat'],'G_M','G_I','G_R','G_F','G_P','G_SP','JS','dWWVdt','JA_W','JA_E','JA_S','JA_N');
JA = JA_S+JA_N+JA_W+JA_E;
I_RES = dWWVdt-JA-G_M-G_R-G_F;


figure;
subplot(2,1,1);
plot(months,G_M/1e6,'-r','linewidth',2);
hold on;
plot(months,G_I/1e6,'-b','linewidth',2);
plot(months,G_F/1e6,'-k','linewidth',2);
plot(months,G_R/1e6,'-','color',[0 0.5 0],'linewidth',2);
plot(months,G_P/1e6,'-m');
plot(months,G_SP/1e6,'-c');
plot(months,I_RES/1e6,'--b','linewidth',2);
xlabel('Month of 1998');
ylabel('WMT (Sv)');
title('Diabatic WMT terms in WWV region $19.75^\circ$C, 025deg\_jra55\_iaf output039');
legend('Vertical Mixing','Numerical Mixing','Surface Forcing', ...
       'Neutral Diffusion','PME Surface Forcing','SW Penetration', ...
       'Numerical Mixing by-Residual');
xlim([0.5 12.5]);
ylim([-50 30]);
grid on;
set(gca,'xtick',1:12);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul', ...
                    'Aug','Sep','Oct','Nov','Dec'});

subplot(2,1,2);
plot(months,JA/1e6,'-k','linewidth',2);
hold on;
plot(months,dWWVdt/1e6,'-b','linewidth',2);
plot(months,JS/1e6,'-r','linewidth',2);
xlabel('Month of 1998');
ylabel('Volume flux (Sv)');
title('Adiabatic WMT terms in WWV region $19.75^\circ$C, 025deg\_jra55\_iaf output039');
legend('Lateral Transport','dWWV/dt','P-E+R');
xlim([0.5 12.5]);
ylim([-50 30]);
grid on;
set(gca,'xtick',1:12);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul', ...
                    'Aug','Sep','Oct','Nov','Dec'});




