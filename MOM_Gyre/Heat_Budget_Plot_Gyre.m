% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

outs = [1:4];
RUNS = { ...
% $$$          {'MOM_Gyre',outs}, ...
% $$$          {'MOM_Gyre_Run002',outs}, ...
% $$$          {'MOM_Gyre_Run003',[1:0]}, ...
% $$$          {'MOM_Gyre_Run006',outs,'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=0$'}, ...
% $$$          {'MOM_Gyre_Run004',outs,'$k_{smag}=2$, $\kappa_B=1\times10^{-6}$, $\kappa_L=0$'}, ...
% $$$          {'MOM_Gyre_Run007',outs}, ...
% $$$          {'MOM_Gyre_Run009',outs,'$k_{smag}=20$, $\kappa_B=5\times10^{-5}$, $\kappa_L=0$'}, ...
% $$$          {'MOM_Gyre_Run009a',outs}, ...
% $$$          {'MOM_Gyre_Run010',outs}, ...
% $$$          {'MOM_Gyre_Run011',outs,'$k_{smag}=20$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run012',outs,'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run013',[0],'$k_{smag}=2$, $\kappa_B=5\times10^{-5}$, $\kappa_L=500$m$^2$s$^{-1}$, SRSTR'}, ...
         {'MOM_Gyre_Run014',outs,'Control'}, ...
% $$$          {'MOM_Gyre_Run015',outs,'dt=1800s'}, ...
% $$$          {'MOM_Gyre_Run016',outs,'$\kappa_R=300$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run017',outs,'$\kappa_R=300$m$^2$s$^{-1}$, $\kappa_v = 10^{-4}$m$^2$s$^{-1}$'}, ...
% $$$          {'MOM_Gyre_Run019',outs,'$\kappa_v = 10^{-6}$m$^2$s$^{-1}$'}, ...
         {'MOM_Gyre_Run018',[5:24],'Double Res.'}, ...
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
          {F(:,months,:), 'Forcing $\mathcal{F}$','k',lthic(rr),ltype{rr}}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',lthic(rr),ltype{rr}}, ...
          {I(:,months,:), 'Numerical Mixing $\mathcal{I}$','b',lthic(rr),ltype{rr}}, ...
% $$$           {L(:,months,:), 'Lateral Mixing $\mathcal{L}$',[0 0.5 0],lthic(rr),ltype{rr}}, ...
% $$$           {NUM(:,months,:), 'Numerical Mixing Direct $\mathcal{I}$','c',lthic(rr),ltype{rr}}, ...
          {R(:,months,:), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],lthic(rr),ltype{rr}}, ...
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
ylabel('Heat flux into fluid warmer than $\Theta$ (TW)');
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
rr = 1;
model = RUNS{rr}{1};%'MOM_Gyre_Run011';
outputs = RUNS{rr}{2};
label = RUNS{rr}{3};
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
region = '';
VARS = {'FlM','FlI'};
TYPE = 'VertInt';
Tls = [20 18 12];
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
    

subplot(length(Tls),2,2*(Ti-1)+vi);
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
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
colormap(cmap);
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
if (Ti==1)
    if (vi == 1)
        title('Vertical Mixing');
    else
        title('Numerical Mixing');
    end
end
text(0.25,min(min(lat))+1,[sprintf('%02d',Tl)  '$^\circ$C']);
if (Ti==1)
    text(0.25,max(max(lat))-2,label);
end
end
end

%%%% Solution general plot:
clear all;
close all;
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
basenc = '/srv/ccrc/data03/z3500785/mom/MOM_Gyre/Run009/';
model = 'MOM_Gyre_Run009';
label = '$k_{smag}=2$, $\kappa_B=1\times10^{-6}$, $\kappa_L=0$';
outputs = [1:5];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
Tl = 12;
[tmp ind] = min(abs(Te-Tl));
tsc = 1e9;

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
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
subplot(2,3,[1 4]);
pcolPlot(lon,lat,EKE);
hold on;
[c,h] = contour(lon,lat,BT,[2:2:30],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,BT,[-30:2:-2],'--k');
clabel(c,h);
caxis([0 0.15]);
colormap(gca,'parula');
title('EKE (m$^2$s$^{-2}$) and BT streamfunction (Sv) Average');
xlabel('x ($^\circ$)');
ylabel('y ($^\circ$)');
subplot(2,3,[2 5]);
pcolPlot(lon,lat,u);
hold on;
[c,h] = contour(lon,lat,SSTsl,[15:0.1:25],'-k');
clabel(c,h);
caxis([-1 1]);
colormap(gca,'redblue');
title('Zonal velocity (ms$^{-1}$) and SST ($^\circ$C) Single Month');
xlabel('x ($^\circ$)');
ylabel('y ($^\circ$)');

Tvar = tempSLrms.^2-tempSL.^2;
subplot(2,3,3);
[X,Y] = ndgrid(yt,z);
pcolPlot(X,Y,Tvar);
hold on;
[c,h] = contour(X,Y,tempSL,[5:1:25],'-k');
clabel(c,h);
colormap(gca,'parula');
ylim([0 1000]);
set(gca,'ydir','reverse');
ylabel('Depth (m)');
xlabel('y ($^\circ$)');
title('Temperature Variance ($^\circ$C$^2$) and temperature ($^\circ$C) x=L/2');
caxis([0 0.05]);

%%%% Check the spatial calculation:
clear all;
close all;
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
basenc = '/srv/ccrc/data03/z3500785/mom/MOM_Gyre/Run004/';
model = 'MOM_Gyre_Run004';
outputs = [1:5];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
Tl = 12;
[tmp ind] = min(abs(Te-Tl));
tsc = 1e9;

JIs = zeros(xL,yL);
QIs = zeros(xL,yL);
dVdts = zeros(xL,yL);
dHdts = zeros(xL,yL);
difts = zeros(xL,yL);
ndifs = zeros(xL,yL);

for ii = 1:length(outputs)
    output = outputs(ii);
    base = [basenc sprintf('output%03d/',output)];
    wname = [base 'ocean_wmass.nc'];
    
    for ti=1:tL
    
    JI = zeros(xL,yL);
    QI = zeros(xL,yL);
    dVdt = zeros(xL,yL);
    dHdt = zeros(xL,yL);
    dift = zeros(xL,yL);

% $$$     for Ti=TL:-1:ind
    for Ti=1:ind
        sprintf('Calculating numdif time %03d of %03d, Temp %03d of %03d',ti,tL,Ti,TL)
        dVdt = dVdt + double(ncread(wname,'dVdt',[1 1 Ti ti],[xL yL 1 1]))*1e9/rho0./area;
        dHdt = dHdt + double(ncread(wname,'dHdt',[1 1 Ti ti],[xL yL 1 1]))./area;
        dift = dift + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);

        JI(2:end,2:end) = JI(2:end,2:end) + (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                             +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end) = JI(1,2:end) + (- txtrans(1,2:end) ...
                                     +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
        JI(2:end,1) = JI(2:end,1) + (txtrans(1:(end-1),1) - txtrans(2:end,1) ...
                       - tytrans(2:end,1))./area(2:end,1);        
        JI(1,1) = JI(1,1) + (-txtrans(1,1)-tytrans(1,1))./area(1,1);
        
        QI(2:end,2:end) = QI(2:end,2:end) + (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                             +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end) = QI(1,2:end) + (- qxtrans(1,2:end) ...
                                     +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        
        QI(2:end,1) = QI(2:end,1) + (qxtrans(1:(end-1),1) - qxtrans(2:end,1) ...
                       - qytrans(2:end,1))./area(2:end,1);        
        QI(1,1) = QI(1,1) + (-qxtrans(1,1)-qytrans(1,1))./area(1,1);

    end
    ndif = dHdt - (dVdt - JI)*rho0*Cp*Te(Ti) - dift - QI;
    
    JIs = JIs+JI*ndays(ti)/sum(ndays);
    QIs = QIs+QI*ndays(ti)/sum(ndays);
    dVdts = dVdts+dVdt*ndays(ti)/sum(ndays);
    dHdts = dHdts+dHdt*ndays(ti)/sum(ndays);
    difts = difts+dift*ndays(ti)/sum(ndays);
    ndifs = ndifs+ndif*ndays(ti)/sum(ndays);
end
end
JIs = JIs/length(outputs);
QIs = QIs/length(outputs);%+QI*ndays(ti)/sum(ndays);
dVdts = dVdts/length(outputs);%+dVdt*ndays(ti)/sum(ndays);
dHdts = dHdts/length(outputs);%+dHdt*ndays(ti)/sum(ndays);
difts = difts/length(outputs);%+dift*ndays(ti)/sum(ndays);
ndifs = ndifs/length(outputs);%+ndif*ndays(ti)/sum(ndays);


figure;
subplot(3,3,1);
pcolPlot(lon,lat,dHdts);
caxis([-200 200]);
title('dHdt');
subplot(3,3,2);
pcolPlot(lon,lat,-dVdts*rho0*Cp*Te(ind));
caxis([-200 200]);
title('-dVdt*rho0*Cp*T');
subplot(3,3,3);
pcolPlot(lon,lat,dHdts-dVdts*rho0*Cp*Te(ind));
caxis([-20 20]);
title('dHdt-dVdt*rho0*Cp*T');
subplot(3,3,4);
pcolPlot(lon,lat,JIs*rho0*Cp*Te(ind));
caxis([-2 2]*1e3);
title('JI*rho0*Cp*T');
subplot(3,3,5);
pcolPlot(lon,lat,-QIs);
caxis([-2 2]*1e3);
title('-QI');
subplot(3,3,6);
pcolPlot(lon,lat,-QIs+JIs*rho0*Cp*Te(ind));
caxis([-20 20]);
title('-QI+JI*rho0*Cp*T');
subplot(3,3,7);
pcolPlot(lon,lat,-difts);
caxis([-20 20]);
title('-M');
subplot(3,3,8);
pcolPlot(lon,lat,dHdts - (dVdts - JIs)*rho0*Cp*Te(ind) - difts - QIs);
caxis([-20 20]);
title('I resid mean');
subplot(3,3,9);
pcolPlot(lon,lat,ndifs);
caxis([-20 20]);
title('I mean resid');
colormap(redblue);

