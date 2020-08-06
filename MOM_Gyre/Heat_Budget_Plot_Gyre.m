% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

outs = [1:5];
RUNS = { ...
% $$$          {'MOM_Gyre',outs}, ...
% $$$          {'MOM_Gyre_Run002',outs}, ...
% $$$          {'MOM_Gyre_Run003',[1:0]}, ...
         {'MOM_Gyre_Run004',outs,'$k_{smag}=2$, $\kappa_B=1\times10^{-6}$'}, ... %, $\kappa_L=0$'}, ...
         {'MOM_Gyre_Run006',outs,'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$'}, ... %$\kappa_L=0$'}, ...
% $$$          {'MOM_Gyre_Run007',outs}, ...
         {'MOM_Gyre_Run009',outs,'$k_{smag}=20$, $\kappa_B=5\times10^{-5}$'}, ... %, $\kappa_L=0$'}, ...
% $$$          {'MOM_Gyre_Run009a',outs}, ...
% $$$          {'MOM_Gyre_Run010',outs}, ...
% $$$          {'MOM_Gyre_Run011',outs,'$k_{smag}=20$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run012',outs,'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run013',[0],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
% $$$          {'MOM_Gyre_Run014',outs,'Control'}, ...
% $$$ % $$$          {'MOM_Gyre_Run015',outs,'dt=1800s'}, ...
% $$$          {'MOM_Gyre_Run016',outs,'$\kappa_R=300$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run017',outs,'$\kappa_R=300$m$^2$s$^{-1}$, $\kappa_v = 10^{-4}$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run019',outs,'$\kappa_v = 10^{-6}$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run020',outs,'Salt H IC'}, ...
% $$$          {'MOM_Gyre_Run018',[1:12],'Double Res.'}, ...
% $$$          {'MOM_Gyre_Run021',[1:12],'Double Res. SSH smooth'}, ...
% $$$          {'MOM_Gyre_Run013',[2],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
% $$$          {'MOM_Gyre_Run013',[3],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
% $$$          {'MOM_Gyre_Run013',[4],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
% $$$          {'MOM_Gyre_Run013',[5],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
%         {'MOM_Gyre_Run013',outs,'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
       };
ltype = {'-','--',':','-.','-','--',':','-.'};
lthic = [2,2,2,2,1];%1,1,1,1,2,2,2];
          
figure;
set(gcf,'Position',[207          97        1609         815]);

rr = 1;
for rr=1:length(RUNS)
clearvars -except base RUNS ltype lthic rr;
outputs = RUNS{rr}{2};
model = RUNS{rr}{1};
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
region = '';


%% Global Calculations:
for i=1:length(outputs)
    
    load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
    
% Fluxes:
M(:,:,i) = GWB.VDF; % Vertical mixing flux (W)
if (isfield(GWB,'RED')) % Redi Diffusion
    R(:,:,i) = GWB.K33+GWB.RED; 
else
    R(:,:,i) = zeros(size(M(:,:,i)));
end
if (isfield(GWB,'VDS')) % Surface Restoring
    F(:,:,i) = GWB.VDS; 
else
    F(:,:,i) = zeros(size(M(:,:,i)));
end
if (isfield(GWB,'LTD')) % Lateral Diffusion
    L(:,:,i) = GWB.LTD; 
else
    L(:,:,i) = zeros(size(M(:,:,i)));
end
D(:,:,i) = GWB.TEN-GWB.ADV; % Material derivative of T (W)
NUM(:,:,i) = GWB.NUM; % Numerical mixing from heat bug

dVdt(:,:,i) = GWB.dVdt; % V Change (m3s-1)
dHdt(:,:,i) = GWB.dHdt; % H Change (W)

% Water-mass transformation:
G(:,:,i) = dVdt(:,:,i); %Water-mass transformation (m3s-1)

% Across-isotherm advective heat flux:
CIA(:,:,i) = G(:,:,i).*repmat(Te,[1 tL])*rho0*Cp;

% External HC Tendency:
EHC(:,:,i) = dVdt(:,:,i).*repmat(Te,[1 tL])*rho0*Cp;

% Internal HC Tendency:
N(:,:,i) = dHdt(:,:,i) - EHC(:,:,i);

% Implicit mixing:
I(:,:,i) = N(:,:,i) - M(:,:,i) - R(:,:,i) - L(:,:,i)-F(:,:,i);

% Non-advective flux into volume:
B(:,:,i) = M(:,:,i)+I(:,:,i)+R(:,:,i)+L(:,:,i)+F(:,:,i);

% WMT from B:
WMTM(:,:,i) = -diff(M(:,:,i),[],1)/dT/rho0/Cp;
WMTI(:,:,i) = -diff(I(:,:,i),[],1)/dT/rho0/Cp;
WMTR(:,:,i) = -diff(R(:,:,i),[],1)/dT/rho0/Cp;
WMTL(:,:,i) = -diff(L(:,:,i),[],1)/dT/rho0/Cp;
WMTF(:,:,i) = -diff(F(:,:,i),[],1)/dT/rho0/Cp;
WMT(:,:,i) = WMTM(:,:,i)+WMTI(:,:,i)+WMTR(:,:,i)+WMTL(:,:,i)+WMTF(:,:,i);

% WMT HB from B:
HWMTM(:,:,i) = rho0*Cp*WMTM(:,:,i).*repmat(T,[1 tL]);
HWMTI(:,:,i) = rho0*Cp*WMTI(:,:,i).*repmat(T,[1 tL]);
HWMTR(:,:,i) = rho0*Cp*WMTR(:,:,i).*repmat(T,[1 tL]);
HWMTL(:,:,i) = rho0*Cp*WMTL(:,:,i).*repmat(T,[1 tL]);
HWMTF(:,:,i) = rho0*Cp*WMTF(:,:,i).*repmat(T,[1 tL]);
HWMT(:,:,i) = HWMTM(:,:,i)+HWMTI(:,:,i)+HWMTR(:,:,i)+HWMTL(:,:,i)+HWMTF(:,:,i);

end

months = 1:length(M(1,:,1));
% $$$ months = 80:120;
% $$$ months = [12:24];

%%%%Heat Flux:
% Production fields:
fields = { ...
          {N(:,months,:), 'Internal HC Tendency $\mathcal{N}$','m',lthic(rr),ltype{rr}}, ...
% $$$           {F(:,months,:), 'Forcing $\mathcal{F}$','k',lthic(rr),ltype{rr}}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',lthic(rr),ltype{rr}}, ...
          {I(:,months,:), 'Numerical Mixing $\mathcal{I}$','b',lthic(rr),ltype{rr}}, ...
% $$$           {L(:,months,:), 'Lateral Mixing $\mathcal{L}$',[0 0.5 0],lthic(rr),ltype{rr}}, ...
% $$$           {NUM(:,months,:), 'Numerical Mixing Direct $\mathcal{I}$','c',lthic(rr),ltype{rr}}, ...
% $$$           {R(:,months,:), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],lthic(rr),ltype{rr}}, ...
% $$$           {dHdt(:,months,:), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',lthic(rr),'--'}, ...
          };

Fscale = 1/1e12;

%Fluxes only:
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
    if (length(fields{i}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
% $$$     for j=1:length(P(1,1,:));
% $$$         h = plot(Te,monmean(fields{i}{1}(:,:,j),2,ndays(months))*Fscale,fields{i}{5}, 'color',fields{i}{3} ...
% $$$              ,'linewidth',0.5);
% $$$     end
    legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    leg{i} = fields{i}{2};
% $$$     leg{i} = strrep(RUNS{rr}{1},'_',' ');
end
ylim([-15 2]);
xlim([0 24]);
box on; 
grid on;
ylabel('Diathermal heat flux (TW, $0.45$Wm$^{-2}$)');
xlabel('Temperature $\Theta$ ($^\circ$C)');
if (rr == 1)
    lg = legend(legh,leg);
end
% $$$ set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

end

% $$$ %%% Distributions change in time:
% $$$ figure;
% $$$ set(gcf,'Position',[3    40   956   963]);%207          97        1609         815]);
% $$$ % $$$ axs1 = subplot(2,1,1);
% $$$ % $$$ hold on;
% $$$ axs2 = gca;%subplot(2,1,2);
% $$$ hold on;
% $$$ for rr=1:length(RUNS)
% $$$ outputs = RUNS{rr}{2};
% $$$ model = RUNS{rr}{1};
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ 
% $$$ for i=1:length(outputs)
% $$$ load([base model sprintf('_output%03d',outputs(i)) ...
% $$$       '_VHza.mat']);
% $$$ Vs(:,:,i) = squeeze(nansum(V,1));
% $$$ Hs(:,:,i) = squeeze(nansum(H,1));
% $$$ end
% $$$ 
% $$$ V = cat(1,cumsum(Vs,1,'reverse'),zeros(1,length(Vs(1,:,1)),length(Vs(1,1,:))));
% $$$ H = cat(1,cumsum(Hs,1,'reverse'),zeros(1,length(Vs(1,:,1)),length(Vs(1,1,:))));
% $$$ HE = rho0*Cp*V.*repmat(Te,[1 length(Vs(1,:,1)) length(Vs(1,1,:))]);
% $$$ HI = H - HE;    
% $$$ 
% $$$ % $$$ if (rr==1)
% $$$ % $$$     plot(Te,V(:,1,1),'-k','Parent',axs1);
% $$$ % $$$ end
% $$$ % $$$ plot(Te,V(:,end,end),ltype{rr},'color','r','Parent',axs1);
% $$$ if (rr==1)
% $$$     plot(Te,H(:,1,1),'-k','Parent',axs2,'linewidth',2);
% $$$     plot(Te,HI(:,1,1),'-r','Parent',axs2,'linewidth',2);
% $$$     plot(Te,HE(:,1,1),'-b','Parent',axs2,'linewidth',2);
% $$$ end
% $$$ plot(Te,H(:,end,end),ltype{rr},'color','k','Parent',axs2);
% $$$ plot(Te,HI(:,end,end),ltype{rr},'color','r','Parent',axs2);
% $$$ plot(Te,HE(:,end,end),ltype{rr},'color','b','Parent',axs2);
% $$$ end

%%% Spatial Structure:
rr = 3;
model = RUNS{rr}{1};%'MOM_Gyre_Run011';
outputs = RUNS{rr}{2};
label = RUNS{rr}{3};
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
region = '';
VARS = {'FlI'};
TYPE = 'VertInt';
Tls = [20 14];
months = 1:tL;%length(M(1,:,1));
% $$$ months = 80:tL;%length(M(1,:,1));

%Mean of all months:
figure;
set(gcf,'Position',[3          59        1356         914]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
for Ti = 1:length(Tls)
    Tl = Tls(Ti);
for vi = 1:length(VARS)
    VAR = VARS{vi};
name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
eval(['load(name,''' VAR ''');']);
eval([VAR '(isnan(' VAR ')) = 0.0;']);
eval([VAR 'a = ' VAR ';']);
for i=2:length(outputs)
    name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
    eval(['load(name,''' VAR ''');']);
    eval([VAR '(isnan(' VAR ')) = 0.0;']);
    eval([VAR 'a = ' VAR 'a + ' VAR ';']);
end
eval([VAR ' = ' VAR 'a/length(outputs);']);
eval([VAR '(' VAR '==0) = NaN;']);
eval(['FlM = ' VAR ';']);

LAND = zeros(size(FlM(:,:,1)));

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;

clim = [-10 10];
sp = 0.1;
clim = [-50 50];
sp = 0.1;

cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

cmap = flipud(lbmap(npts-3,'RedBlue'));
cmap = redblue(npts-3);
climn = clim;
    
% $$$ subplot(length(Tls),2,2*(Ti-1)+vi);
subplot(1,4,Ti+2);
X = lon(xvec,yvec);
Y = lat(xvec,yvec);
tmp = FlM;
tmp(isnan(tmp)) = 0.0;
Z = monmean(tmp(:,:,months),3,ndays(months));
Z(Z == 0) = NaN;
Z = Z(xvec,yvec);
    
Z(Z<clim(1)) = clim(1);
contourf(X,Y,Z,cpts,'linestyle','none');
hold on;    
% $$$ plot(X(:),Y(:),'.k','MarkerSize',5);
caxis(climn);
cb = colorbar('Location','NorthOutside');
ylabel(cb,'$\mathcal{I}$ (Wm$^{-2}$)');
colormap(cmap);
xlabel('Longitude ($^\circ$E)');
if (Ti==1)
    ylabel('Latitude ($^\circ$N)');
else
    set(gca,'yticklabel',[]);
end
% $$$ if (Ti==1)
% $$$     if (vi == 1)
% $$$         title('Vertical Mixing');
% $$$     else
        title('Numerical Mixing');
% $$$     end
% $$$ end
text(0.25,min(min(lat))+0.2,[sprintf('%02d',Tl)  '$^\circ$C'],'BackgroundColor','w');
% $$$ if (Ti==1)
% $$$     text(0.25,max(max(lat))-0.25,label);
% $$$ end
title(label);
end
end

%%%% Solution general plot:
clear all;
close all;
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
basenc = '/srv/ccrc/data03/z3500785/mom/MOM_Gyre/Run009/';
model = 'MOM_Gyre_Run009';
label = '$k_{smag}=2$, $\kappa_B=1\times10^{-6}$, $\kappa_L=0$';
outputs = [1:4];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
Tl = 12;
[tmp ind] = min(abs(Te-Tl));
tsc = 1e9;

taux = squeeze(ncread([basenc 'tau.nc'],'taux',[1 1 1],[1 yL 1]));
txtrans = zeros(xL,yL);
u = txtrans;
u_rms = txtrans;
v = txtrans;
v_rms = txtrans;
tempSL = zeros(yL,zL);
tempSLrms = zeros(yL,zL);

for ii=1:length(outputs)
    output = outputs(ii)
    base = [basenc sprintf('output%03d/',output)];
    fname = [base 'ocean.nc'];
    
    txtrans = txtrans + mean(sum(ncread(fname,'tx_trans'),3),4);
    u = u + mean(ncread(fname,'u',[1 1 1 1],[xL yL 1 tL]),4);
    u_rms = u_rms + mean(ncread(fname,'u_rms',[1 1 1 1],[xL yL 1 tL]),4);
    v = v + mean(ncread(fname,'v',[1 1 1 1],[xL yL 1 tL]),4);
    v_rms = v_rms + mean(ncread(fname,'v_rms',[1 1 1 1],[xL yL 1 tL]),4);
    tempSL = tempSL + squeeze(mean(ncread(fname,'temp',[15 1 1 1],[1 yL zL tL]),4));
    tempSLrms = tempSLrms + squeeze(mean(ncread(fname,'temp_rms',[15 1 1 1],[1 yL zL tL]),4));
end
vars = {'txtrans','u','u_rms','v','v_rms','tempSL','tempSLrms'};
for ii=1:length(vars)
    eval([vars{ii} ' = ' vars{ii} '/length(outputs);']);
end
EKE = (u_rms.^2-u.^2)+(v_rms.^2-v.^2);
BT = cumsum(txtrans,2);

SSTsl = ncread(fname,'temp',[1 1 1 24],[xL yL 1 1]);
u = ncread(fname,'u',[1 1 1 24],[xL yL 1 1]);


figure;
set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
subplot(2,6,[1 7]);
plot(taux,lat(1,:),'-k');
ylim([min(lat(1,:)) max(lat(1,:))]);
ylabel('y ($^\circ$)');
xlabel('x ($^\circ$)');
title('Wind Stress (Nm$^{-2}$)');
subplot(2,6,[2 3 8 9]);
contourf(lon,lat,EKE,[0:0.002:0.12],'linestyle','none');
hold on;
[c,h] = contour(lon,lat,BT,[2:2:30],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,BT,[-30:2:-2],'--k');
clabel(c,h);
caxis([0 0.12]);
cmap = redblue(100);
cmap = cmap(1:50,:)
colormap(gca,flipud(cmap));
xlabel('x ($^\circ$)');
cb = colorbar('Location','NorthOutside');
ylabel(cb,'EKE (m$^2$s$^{-2}$) and BT $\Psi$ (Sv)','FontSize',15);
set(gca,'yticklabel',[]);
plot(lon(:),lat(:),'.k');
ylim([min(lat(1,:)) max(lat(1,:))]);
% $$$ subplot(2,3,[2 5]);
% $$$ pcolPlot(lon,lat,u);
% $$$ hold on;
% $$$ plot(lon(:),lat(:),'.k');
% $$$ [c,h] = contour(lon,lat,SSTsl,[15:0.1:25],'-k');
% $$$ clabel(c,h);
% $$$ caxis([-1 1]);
% $$$ colormap(gca,'redblue');
% $$$ title('Zonal velocity (ms$^{-1}$) and SST ($^\circ$C) Single Month');
% $$$ xlabel('x ($^\circ$)');
% $$$ ylabel('y ($^\circ$)');

Tvar = tempSLrms.^2-tempSL.^2;
subplot(2,6,[4 5 6]);
[X,Y] = ndgrid(yt,z);
contourf(X,Y,Tvar,[0:0.001:0.05],'linestyle','none');
hold on;
[c,h] = contour(X,Y,tempSL,[5:1:25],'-k');
clabel(c,h);
plot(X(:),Y(:),'.k');
colormap(gca,flipud(cmap));%'parula');
caxis([0 0.05]);
ylim([0 600]);
set(gca,'ydir','reverse');
ylabel('Depth (m)');
xlabel('y ($^\circ$)');
cb = colorbar('Location','NorthOutside');
ylabel(cb,'$\Theta$ Variance ($^\circ$C$^2$) and $\Theta$ ($^\circ$C), x=$5^\circ$');
caxis([0 0.05]);

%% Plot vertical structure:


%%%% Check the spatial calculation:
clear all;
close all;
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
basenc = '/srv/ccrc/data03/z3500785/mom/MOM_Gyre/Run004/';
model = 'MOM_Gyre_Run004';
haveKPP = 0;
outputs = [1:5];
months = [1:24];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
tsc = 1e9;

JIs = zeros(xL,yL,TL+1);
QIs = zeros(xL,yL,TL+1);
dVdts = zeros(xL,yL,TL+1);
dHdts = zeros(xL,yL,TL+1);
difts = zeros(xL,yL,TL+1);
ndifs = zeros(xL,yL,TL+1);

% $$$ TENs = zeros(xL,yL,TL+1);
% $$$ ADVs = zeros(xL,yL,TL+1);
% $$$ DIFs = zeros(xL,yL,TL+1);

for ii = 1:length(outputs)
    output = outputs(ii);
    base = [basenc sprintf('output%03d/',output)];
    wname = [base 'ocean_wmass.nc'];
    hname = [base 'ocean_heat.nc'];
    
    for ti=months
    
    JI = zeros(xL,yL,TL+1);
    QI = zeros(xL,yL,TL+1);
    dVdt = zeros(xL,yL,TL+1);
    dHdt = zeros(xL,yL,TL+1);
    dift = zeros(xL,yL,TL+1);
    ndif = zeros(xL,yL,TL+1);

    % warm-to-cold:
    for Ti=TL:-1:1
        Tm = Ti+1;
        Tf = Ti;
% $$$     % cold-to-warm:
% $$$     for Ti = 2:(TL+1)
% $$$         Tm = Ti-1;
% $$$         Tf = Ti-1;
        sprintf('Calculating numdif time %03d of %03d, Temp %03d of %03d, out %03d of %03d',ti-months(1)+1,length(months),Ti,TL,output,length(outputs))
        dVdt(:,:,Ti) = dVdt(:,:,Tm) + double(ncread(wname,'dVdt',[1 1 Tf ti],[xL yL 1 1]))*1e9/rho0./area;
        dHdt(:,:,Ti) = dHdt(:,:,Tm) + double(ncread(wname,'dHdt',[1 1 Tf ti],[xL yL 1 1]))./area;
        dift(:,:,Ti) = dift(:,:,Tm) + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Tf ti],[xL yL 1 1]) + ...
                      ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Tf ti],[xL yL 1 1]);
        if (haveKPP)
            dift(:,:,Ti) = dift(:,:,Ti) + ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Tf ti],[xL yL 1 1]);
        end
        
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Tf ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Tf ti],[xL yL 1 1])*tsc/rho0;
        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Tf ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Tf ti],[xL yL 1 1]);

        JI(2:end,2:end,Ti) = JI(2:end,2:end,Tm) + (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                             +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end,Ti) = JI(1,2:end,Tm) + (- txtrans(1,2:end) ...
                                     +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
        JI(2:end,1,Ti) = JI(2:end,1,Tm) + (txtrans(1:(end-1),1) - txtrans(2:end,1) ...
                       - tytrans(2:end,1))./area(2:end,1);        
        JI(1,1,Ti) = JI(1,1,Tm) + (-txtrans(1,1)-tytrans(1,1))./area(1,1);
        
        QI(2:end,2:end,Ti) = QI(2:end,2:end,Tm) + (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                             +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end,Ti) = QI(1,2:end,Tm) + (- qxtrans(1,2:end) ...
                                     +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        
        QI(2:end,1,Ti) = QI(2:end,1,Tm) + (qxtrans(1:(end-1),1) - qxtrans(2:end,1) ...
                       - qytrans(2:end,1))./area(2:end,1);        
        QI(1,1,Ti) = QI(1,1,Tm) + (-qxtrans(1,1)-qytrans(1,1))./area(1,1);

        ndif(:,:,Ti) = dHdt(:,:,Ti) - (dVdt(:,:,Ti) - JI(:,:,Ti))*rho0*Cp*Te(Ti) - dift(:,:,Ti) - QI(:,:,Ti);
    end
    
    JIs = JIs+JI*ndays(ti)/sum(ndays(months));
    QIs = QIs+QI*ndays(ti)/sum(ndays(months));
    dVdts = dVdts+dVdt*ndays(ti)/sum(ndays(months));
    dHdts = dHdts+dHdt*ndays(ti)/sum(ndays(months));
    difts = difts+dift*ndays(ti)/sum(ndays(months));
    ndifs = ndifs+ndif*ndays(ti)/sum(ndays(months));
    
% $$$     % Eulerian terms:
% $$$     TENs = TENs + sum(ncread(hname,'temp_tendency',[1 1 1 ti],[xL yL 50 1]),3)*ndays(ti)/sum(ndays(months));
% $$$     ADVs = ADVs + sum(ncread(hname,'temp_advection',[1 1 1 ti],[xL yL 50 1]),3)*ndays(ti)/sum(ndays(months));
% $$$     DIFs = DIFs + sum(ncread(hname,'temp_vdiffuse_diff_cbt',[1 1 1 ti],[xL yL 50 1]),3)*ndays(ti)/sum(ndays(months));
    end
end
JI = JIs/length(outputs);
QI = QIs/length(outputs);%+QI*ndays(ti)/sum(ndays);
dVdt = dVdts/length(outputs);%+dVdt*ndays(ti)/sum(ndays);
dHdt = dHdts/length(outputs);%+dHdt*ndays(ti)/sum(ndays);
dift = difts/length(outputs);%+dift*ndays(ti)/sum(ndays);
ndif = ndifs/length(outputs);%+ndif*ndays(ti)/sum(ndays);

% $$$ TENs = TENs/length(outputs);%+ndif*ndays(ti)/sum(ndays);
% $$$ ADVs = ADVs/length(outputs);%+ndif*ndays(ti)/sum(ndays);
% $$$ DIFs = DIFs/length(outputs);%+ndif*ndays(ti)/sum(ndays);

Tl = 14.5;
[tmp ind1] = min(abs(Te-Tl));
Tl = 14.5;
[tmp ind2] = min(abs(Te-Tl));
% $$$ ctow = 
% $$$ if 
JIs = mean(JI(:,:,ind1:ind2)*rho0*Cp.*repmat(permute(Te(ind1:ind2),[3 2 1]),[xL yL 1]),3);
QIs = mean(QI(:,:,ind1:ind2),3);
dVdts = mean(-dVdt(:,:,ind1:ind2)*rho0*Cp.*repmat(permute(Te(ind1:ind2),[3 2 1]),[xL yL 1]),3);
dHdts = mean(dHdt(:,:,ind1:ind2),3);
difts = mean(dift(:,:,ind1:ind2),3);
ndifs = mean(ndif(:,:,ind1:ind2),3);

figure;
subplot(3,4,1);
pcolPlot(lon,lat,dHdts);
caxis([-100 100]);
title('dHdt');
subplot(3,4,2);
pcolPlot(lon,lat,-QIs);
caxis([-2 2]*1e3);
title('-QI');
subplot(3,4,3);
pcolPlot(lon,lat,dHdts-QIs);
caxis([-2 2]*1e3);
title('dHdt-QI');
subplot(3,4,5);
pcolPlot(lon,lat,dVdts);
caxis([-20 20]);
title('-dVdt*rho0*Cp*T');
subplot(3,4,6);
pcolPlot(lon,lat,JIs);
caxis([-2 2]*1e3);
title('JI*rho0*Cp*T');
subplot(3,4,7);
pcolPlot(lon,lat,JIs+dVdts);
caxis([-2 2]*1e3);
title('(JI-dVdts)*rho0*Cp*T');
subplot(3,4,4);
pcolPlot(lon,lat,dHdts+dVdts);
caxis([-20 20]);
title('dHdt-dVdt*rho0*Cp*T');
subplot(3,4,8);
pcolPlot(lon,lat,-QIs+JIs);
caxis([-20 20]);
title('-QI+JI*rho0*Cp*T');
subplot(3,4,9);
pcolPlot(lon,lat,-difts);
caxis([-2 2]);
title('-F-M');
subplot(3,4,11);
pcolPlot(lon,lat,ndifs);
caxis([-20 20]);
title('I mean resid');
colormap(redblue);

% $$$ % Vertically integrated heat budget:
% $$$ subplot(3,4,1);
% $$$ pcolPlot(lon,lat,dHdts);
% $$$ %caxis([-200 200]);
% $$$ title('dHdt');
% $$$ subplot(3,4,2);
% $$$ pcolPlot(lon,lat,-QIs);
% $$$ %caxis([-2 2]*1e3);
% $$$ title('-QI');
% $$$ subplot(3,4,3);
% $$$ pcolPlot(lon,lat,dHdts-QIs);
% $$$ %caxis([-2 2]*1e3);
% $$$ title('dHdt-QI');
% $$$ subplot(3,4,4);
% $$$ pcolPlot(lon,lat,-difts);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('-F-M');
% $$$ subplot(3,4,8);
% $$$ pcolPlot(lon,lat,-DIFs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('DIF');
% $$$ subplot(3,4,5);
% $$$ pcolPlot(lon,lat,TENs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('TEN');
% $$$ subplot(3,4,6);
% $$$ pcolPlot(lon,lat,-ADVs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('-ADV');
% $$$ subplot(3,4,7);
% $$$ pcolPlot(lon,lat,TENs-ADVs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('TEN-ADV');
% $$$ subplot(3,4,9);
% $$$ pcolPlot(lon,lat,TENs-dHdts);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('TEN-dHdt');
% $$$ subplot(3,4,10);
% $$$ pcolPlot(lon,lat,-ADVs+QIs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('-ADV+QI');
% $$$ 
% $$$ subplot(3,4,10);
% $$$ pcolPlot(lon,lat,dHdts - (dVdts - JIs)*rho0*Cp*Te(ind) - difts - QIs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('I resid mean');
% $$$ subplot(3,4,11);
% $$$ pcolPlot(lon,lat,ndifs);
% $$$ % $$$ caxis([-20 20]);
% $$$ title('I mean resid');
% $$$ colormap(redblue);
% $$$ 

% $$$ 
