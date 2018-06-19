% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
% MOM-Gyre:
         {'MOM_Gyre',[1]}, ...
         {'MOM_Gyre',[2]}, ...
         {'MOM_Gyre',[3]}, ...
       };

rr = 2;
% $$$ for rr=1:length(RUNS)
outputs = RUNS{rr}{2};
model = RUNS{rr}{1};

clearvars -except base RUNS rr outputs model;
    
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
region = 'Global';

%% Global Calculations:
for i=1:length(outputs)
    
    load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
    
% Fluxes:
M(:,:,i) = GWB.VDF; % Vertical mixing flux (W)
D(:,:,i) = GWB.TEN-GWB.ADV; % Material derivative of T (W)

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
I(:,:,i) = N(:,:,i) - M(:,:,i);

% Non-advective flux into volume:
B(:,:,i) = M(:,:,i)+I(:,:,i);

% WMT from B:
WMTM(:,:,i) = -diff(M(:,:,i),[],1)/dT/rho0/Cp;
WMTI(:,:,i) = -diff(I(:,:,i),[],1)/dT/rho0/Cp;
WMT(:,:,i) = WMTM(:,:,i)+WMTI(:,:,i);

% WMT HB from B:
HWMTM(:,:,i) = rho0*Cp*WMTM(:,:,i).*repmat(T,[1 tL]);
HWMTI(:,:,i) = rho0*Cp*WMTI(:,:,i).*repmat(T,[1 tL]);
HWMT(:,:,i) = HWMTM(:,:,i)+HWMTI(:,:,i);

end

months = 2:length(M(1,:));

%%%%Heat Flux:
% Production fields:
fields = { ...
          {N(:,months,:), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {I(:,months,:), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
          {dHdt(:,months,:), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',2,'--'}, ...
          };

Fscale = 1/1e12;

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
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
ylim([-10 2]);
xlim([0 24]);
box on; 
grid on;
ylabel('Heat flux into fluid warmer than $\Theta$ (TW)');
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%%%WM Transformation / Volume Budget:
fields = { ...
          {dVdt(:,months,:), 'Tendency $\frac{\partial\mathcal{V}}{\partial t}$','m',2,'-'}, ...
          {WMT(:,months,:), 'Total WMT $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
          {WMTM(:,months,:), 'WMT $\mathcal{G}$ from Vertical Mixing','r',2,'-'}, ...
          {WMTI(:,months,:), 'WMT $\mathcal{G}$ from Implicit Mixing','b',2,'-'}, ...
          };

Mscale = 1/1e6;

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
    if (length(fields{i}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Mscale,3),fields{i}{5}, 'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    leg{i} = fields{i}{2};
end
ylim([-2 2]);
xlim([-3 31]);
box on;
grid on;
ylabel('Water Mass Transformation (Sv)');
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%% Spatial Structure:

% $$$ VAR = 'FlMkppish';
VAR = 'FlI';
% $$$ VAR = 'FlSP';
% $$$ VAR = 'WMTP';
% $$$ VAR = 'WMTM';
% $$$ VAR = 'WMTI';
TYPE = 'VertInt';
% $$$ TYPE = 'WMT';
Tl = 18;
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

% $$$ FlMblall = FlM - FlMwave - FlMkppish;
% $$$ FlM = FlMblall;

% CHECK spatial structure sums to total:
% $$$ Tls = [14.75:2.5:27.25]+0.25;
Tls = [0:2:22];
SUM = zeros(size(Tls));
for ii = 1:length(Tls)

    Tl = Tls(ii)
    load([base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']);
    eval([VAR '(isnan(' VAR ')) = 0.0;']);
    eval([VAR 'a = ' VAR ';']);
    for i=2:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']);
        eval([VAR '(isnan(' VAR ')) = 0.0;']);
        eval([VAR 'a = ' VAR 'a + ' VAR ';']);
    end
    eval([VAR ' = ' VAR 'a/length(outputs);']);
    eval([VAR '(' VAR '==0) = NaN;']);
    eval(['FlM = ' VAR ';']);
    tmp = FlM;
    tmp(isnan(tmp)) = 0.0;
    Z = monmean(tmp(:,:,months),3,ndays(months));
    Z(Z == 0) = NaN;
    SUM(ii) = nansum(nansum(area.*Z));
end
% $$$ plot(Tls,SUM/1e6,'Xb','MarkerSize',12,'LineWidth',2); 
plot(Tls,SUM/1e12,'Xb','MarkerSize',12,'LineWidth',2); 

%%% Plot spatial pattern:

LAND = zeros(size(FlM(:,:,1)));

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;

if (rr==1)
    months = [2:12];
else
    months = [2:6];
end

clim = [-10 10];
sp = 0.1;
clim = [-50 50];
sp = 0.5;

cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

cmap = flipud(lbmap(npts-3,'RedBlue'));
cmap = redblue(npts-3);
climn = clim;
    
%Mean of all months:
figure;
set(gcf,'Position',[3          59        1916         914]);
set(gcf,'defaulttextfontsize',20);
set(gcf,'defaultaxesfontsize',20);

subplot(2,2,1);
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
caxis(climn);
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
colormap(cmap);

end
