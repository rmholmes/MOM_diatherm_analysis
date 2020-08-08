% This script makes plots depth-longitude or depth-latitude slices 

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[4567]}, ...
% $$$ % MOM025-SIS:
% $$$     {'MOM025_kb3seg',[101120],'(a) MOM025 Control'}, ...
% $$$     {'MOM025',[16:19],'(b) MOM025-kb0'}, ...
% $$$     {'MOM025_kb1em5',[95:99],'(c) MOM025-kb5'}, ...
% $$$     {'MOM025_kb1em6',[30]}, ...
% $$$ % ACCESS-OM2 Gadi runs:
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[7680],'(a) ACCESS-OM2-025-NG'}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kbvar',[7781],'(b) ACCESS-OM2-025-NG-kbv'}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5',[7781],'(c) ACCESS-OM2-025-NG-kb5'}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[81]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_kds75',[7680]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[80]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[300]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
% $$$          {'ACCESS-OM2_01deg_jra55_ryf',[636639]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf',[31],'(d) ACCESS-OM2-1-KDS50'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_sgl',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31],'(e) ACCESS-OM2-1-GFDL50'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135],'ACCESS-OM2-1-KDS75'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135],'ACCESS-OM2-1-KDS100'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135],'(f) ACCESS-OM2-1-KDS135'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf',[51]}, ...
       };
cols = {'b',[0.3020    0.7451    0.9333],'b','b','b'};
typs = {'-','-','--',':','-.'};

figure;
set(gcf,'Position',[1923           5        1366         998]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

rr = 1;
% $$$ for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

% $$$     clearvars -except base RUNS rr outputs model;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
        ndays = ndays(1:12);
    end
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
    nyrs = tL/12;
    if (round(nyrs)~=nyrs)
        anavg = 1;
        nyrs = tL;
        months = {[1:1]};
    else
        anavg = 0;
        months = {[1:12]};
    end
    yrs = 1:nyrs;
    ycur = 1;

% Load Variable and calculate mean:
reg = 'EqPM2';
load([base model sprintf(['_output%03d_varsat_' reg '_monthVar.mat'],outputs(1))]);

ndays = [31 28 31 30 31 30 31 31 30 31 30 31];
ndays = cat(2,ndays,ndays)

mld = monmean(mld,2,ndays);
temp = monmean(temp,3,ndays);
Tdh = sqrt(monmean(0.5*(Tdxsq+Tdysq),3,ndays));
Tdv = sqrt(monmean(Tdzsq,3,ndays));
udhd = sqrt(monmean(0.5*(udxsq+vdysq),3,ndays));

[xL,zL,tL] = size(temp);
TL = length(T);

if (max(max(max(temp)))>100)
    temp = temp-273.15;
end

VARS = {'udhd','Tdh','Tdv'};
% $$$ names = {'(b) $\sqrt{\frac{1}{4}\left(|\Delta_x u|^2+|\Delta_y u|^2+|\Delta_x v|^2+|\Delta_y v|^2\right)}$', ...
names = {'(c) $\sqrt{\frac{1}{2}\left(|\Delta_x u|^2+|\Delta_y v|^2\right)}$', ...
         '(d) $\sqrt{\frac{1}{2}\left(|\Delta_x T|^2+|\Delta_y T|^2\right)}$', ...
         '(e) $\sqrt{|\Delta_z T|^2}$'};
units = {'$ms^{-1}$','$^\circ C$','$^\circ C$'};
clims = {[0 0.05],[0 0.5],[0 5]};

nlv = 50;
  
cmap = parula(nlv-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
cmap = flipud(cmap);
cpts = {};
    for i=1:length(VARS)
        sp = (clims{i}(2)-clims{i}(1))/(nlv-3);
        cpts{i} = [-1e10 clims{i}(1):sp:clims{i}(2) 1e10];
    end

poss = [0.1300    0.4800    0.4154    0.3355; ...
        0.1300    0.1100    0.4154    0.3355; ...    
        0.1300    0.1100    0.4154    0.3355;];    

set(gcf,'Position',[1923           5        1366         998]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
for i=1:length(VARS)
    subplot(2,2,i);
    eval(['Z = ' VARS{i} ';']);
    Z(Z<clims{i}(1)) = clims{i}(1);
    contourf(Xt,-Zt,Z,cpts{i},'linestyle','none');
    hold on;
    [c,h] = contour(Xt,-Zt,temp,[0:1:35],'-k');
    plot(Xt(1:20:end,:),-Zt(1:20:end,:),'.','color',[0.5 0.5 0.5]);
    clabel(c,h,[0:2:35]);
    plot(Xt(:,1),mld,'--','color',[0 0.5 0],'linewidth',3);
    ylim([-200 0]);
    xlim([-200 -80]);
    cb = colorbar;
    ylabel(cb,units{i});
        xlabel('Longitude ($^\circ$E)');
        ylabel('Depth (m)');
    caxis(clims{i});
    text(-199,-10,names{i},'Backgroundcolor','w','FontSize',15,'margin',0.5);
    set(gca,'Position',poss(i,:));
colormap(gca, cmap);
end

end
