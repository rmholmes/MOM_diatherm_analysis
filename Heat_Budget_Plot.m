% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

% Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_nipoall/mat_data/';
% $$$ model = 'MOM025_nipoall';
% $$$ outputs = [19];
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [2 3 4 5 6];
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_wombat/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [1978];
% $$$ outputs = [2];
% $$$ base = '/srv/ccrc/data03/z3500785/MOM01_HeatDiag/mat_data/';
% $$$ model = 'MOM01';
% $$$ outputs = [333];

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
region = 'Global';

%% Global Calculations:
for i=1:length(outputs)
    
    load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);

% Fluxes:
P(:,:,i) = PME+RMX; % PME effective heat flux (W)
F(:,:,i) = SWH+VDS+FRZ+ETS; % Surface heat flux (W)
M(:,:,i) = VDF+KNL; % Vertical mixing flux (W)
CNV(:,:,i) = KNL; % KPP Non-local, convection (W)
if (exist('RED'))
    R(:,:,i) = RED+K33; % Redi diffusion (W)
    GM(:,:,i) = NGM; % GM (W)
else
    R(:,:,i) = zeros(size(VDF));
    GM(:,:,i) = zeros(size(VDF));
end
D(:,:,i) = TEN-ADV-SUB-GM(:,:,i); % Material derivative of T (W)
SW(:,:,i) = SWH; % Short-wave heat

% Snapshot fields:
dVdt(:,:,i) = diff(Vsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); % V Change (m3s-1)
dHdt(:,:,i) = diff(Hsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); % H Change (W)

% Water-mass transformation:
G(:,:,i) = dVdt(:,:,i) - SFW; %Water-mass transformation (m3s-1)

% Across-isotherm advective heat flux:
CIA(:,:,i) = G(:,:,i).*repmat(Te,[1 tL])*rho0*Cp;

% Volume Change term:
VCT(:,:,i) = dVdt(:,:,i).*repmat(Te,[1 tL])*rho0*Cp;

% Implicit mixing by residual:
I(:,:,i) = dHdt(:,:,i)-D(:,:,i)-CIA(:,:,i);

% Non-advective flux into volume:
B(:,:,i) = F(:,:,i)+M(:,:,i)+I(:,:,i)+R(:,:,i);

% Interior heat source P:
PI(:,:,i) = P(:,:,i) - SFW.*repmat(Te,[1 tL])*rho0*Cp;

% Total flux:
N(:,:,i) = B(:,:,i) + PI(:,:,i);

% Monthly binned total flux:
Nmon(:,:,i) = TENMON;

% Checks:
% $$$ tmp = -diff(B(:,:,i),[],1)/dT.*repmat(T,[1 tL]);
% $$$ CIAch(:,:,i) = zeros(size(B(:,:,i)));
% $$$ CIAch(2:end-1,:,i) = avg(tmp,1);
% $$$ CIAch(:,:,i) = -diff(B(:,:,i),[],1)/dT.*repmat(T,[1 tL]);
end

%%%%Global flux, Annual Average:
months = [1:12];
label = 'January';

% Production fields:
% $$$ fields = { ...
% $$$           {F(:,months,:)+PI(:,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$           {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {N(:,months,:), 'Total $\mathcal{N}$','m',2,'-'}, ...
% $$$           {M(:,months,:)+I(:,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {Nmon(:,months,:), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$           {SW(:,months,:), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$ % $$$           {P(:,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$ % $$$           {PI(:,months,:), 'P-E+R Interior Heat Source',0.5*[1 1 1],2,':'}, ...
% $$$           };

% Wombat fields:
fields = { ...
          {F(:,months,:)+PI(:,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {R(:,months,:), 'Redi Mixing $\mathcal{R}$',[0.5 0 0.5],2,'-'}, ...
% $$$           {GM(:,months,:), 'GM Transport $\mathcal{GM}$',[0.5 0 0.5],2,'-'}, ...
          {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {N(:,months,:), 'Total $\mathcal{N}$','m',2,'-'}, ...
          {M(:,months,:)+I(:,months,:)+R(:,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
          {Nmon(:,months,:), 'Monthly-Binned Total','m',2,'--'}, ...
          {SW(:,months,:), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$           {P(:,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {PI(:,months,:), 'P-E+R Interior Heat Source',0.5*[1 1 1],2,':'}, ...
          };

Fscale = 1/1e15;
Escale = 1/Cp/rho0/1e6;

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
% $$$     for j=1:length(P(1,1,:));
% $$$         h = plot(Te,monmean(fields{i}{1}(:,:,j),2,ndays(months))*Fscale,fields{i}{5}, 'color',0.7*[1 1 1] ...
% $$$              ,'linewidth',0.5);
% $$$     end
    legh(i) = plot(Te,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    leg{i} = fields{i}{2};
end
ylim([-1.5 1.5]);
xlim([-3 31]);
box on;
grid on;
ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
% $$$ title(['MOM025 Global Heat Budget ' label]);
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%%%Global transformation, Annual Average:
months = [1:12];
label = 'January';

% Wombat fields:
fields = { ...
          {F(:,months,:)+PI(:,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {R(:,months,:), 'Redi Mixing $\mathcal{R}$',[0 0.75 0.75],2,'-'}, ...
% $$$           {GM(:,months,:), 'GM Transport $\mathcal{GM}$',[0 0.25 0.25],2,'-'}, ...
          {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {N(:,months,:), 'Total $\mathcal{N}$','m',2,'-'}, ...
          {M(:,months,:)+I(:,months,:)+R(:,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
          {Nmon(:,months,:), 'Monthly-Binned Total','m',2,'--'}, ...
          {SW(:,months,:), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$           {P(:,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {PI(:,months,:), 'P-E+R Interior Heat Source',0.5*[1 1 1],2,':'}, ...
          };

Fscale = 1/1e15;
Escale = 1/Cp/rho0/1e6;

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
% $$$     for j=1:length(P(1,1,:));
% $$$         h = plot(Te,monmean(fields{i}{1}(:,:,j),2,ndays(months))*Fscale,fields{i}{5}, 'color',0.7*[1 1 1] ...
% $$$              ,'linewidth',0.5);
% $$$     end
    legh(i) = plot(T,-mean(monmean(diff(fields{i}{1},[],1)/dT,2,ndays(months))*Escale,3),fields{i}{5}, 'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    leg{i} = fields{i}{2};
end
ylim([-100 100]);
xlim([-3 31]);
box on;
grid on;
ylabel('Water Mass Transformation (Sv)');
% $$$ title(['MOM025 Global Heat Budget ' label]);
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%% Temperature vs. time:
months = [1:12];
fields = {
% $$$           {N(:,months,:), 'Tendency $\mathcal{N}$','m',2,'-'}, ...
% $$$           {F(:,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$           {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {R(:,months,:), 'Redi Mixing $\mathcal{R}$','b',2,'-'}, ...
          {GM(:,months,:), 'GM $\mathcal{GM}$',[0.5 0 0.5],2,'-'}, ...
% $$$           {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {M(:,months,:)+I(:,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {SW(:,months,:), 'Shortwave Penetration',[0 0.5 0],2,'--'}, ...
% $$$           {SW(:,months,:)+M(:,months,:)+I(:,months,:), 'Shortwave Penetration + $\mathcal{M}$ + $\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {P(:,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {S(:,months,:), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {F(:,months,:)+P(:,months,:), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$           {Nmon(:,months,:), 'Monthly-Binned Tendency','m',2,'--'}, ...
% $$$           {Nsnap(:,months,:), 'Calculated Tendency',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(:,months,:), 'Calculated Tendency 2',0.3*[1 1 1],2,'--'}, ...
          };

% Fluxes:
scale = 1/1e15;label = '(PW)';
% $$$ caxs = [-1 0];x = Te;
% $$$ sp = 0.05;
caxs = [-0.2 0];x = Te;
sp = 0.005;
% $$$ caxs = [-0.8 0];x = Te;
% $$$ sp = 0.05;
cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

figure;
%set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'Position',[3    40   956   963]);
for ii=1:length(fields)
    subplot(1,length(fields),ii);
    V = mean(fields{ii}{1},3)'*scale;
    [X,Y] = ndgrid(1:tL,Te);
    contourf(X,Y,V,cint);%,'linestyle','none');
    cb = colorbar('Location','NorthOutside','FontSize',25);    
% $$$     set(gca,'xtick',0.5:1:11.5);
% $$$     set(gca,'xticklabel',[1:12]);
    set(gca,'ytick',-5:5:35);
    set(gca,'xtick',[1:tL]);
    ylim([-3 31]);
    grid on;
    caxis(caxs);
% $$$     caxis([-0.8 0]);
    xlabel('Month');
    ylabel('Temperature ($^\circ$C)');
    xlabel(cb,['MOM025-WOMBAT' fields{ii}{2} ' ' ...
               label],'FontSize',20);
% $$$     xlabel(cb,'PW','FontSize',25);
    set(gca,'FontSize',25);
end
% $$$ colormap(redblue);
cmap = redblue((length(cint)-3)*2);
cmap = cmap(1:(length(cint)-3),:);
colormap(cmap);
% $$$ caxis([-0.8 0]);
% $$$ hold on;
% $$$ plot([1 tL],[22.5 22.5],'--k','linewidth',2);
% $$$ LabelAxes(gca,3,25,0.008,0.965);

%% Global Seasonal Cycle TS
months = 1:12;
Ts = 22.5;
[tmp ind] = min(abs(Te-Ts));

Fscale = 1/1e15;

figure;
%set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'Position',[34          40        1164         963]);
set(gcf,'defaultlinelinewidth',2);


subplot(2,1,1);
fields = {
          {F(ind,months,:)+PI(ind,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$           {M(ind,:months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(ind,:months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {N(ind,months,:), 'Total $\mathcal{N}$','m',2,'-'}, ...
          {M(ind,months,:)+I(ind,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {P(ind,:months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {S(ind,:months,:), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {F(ind,:months,:)+P(ind,:months,:), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$           {Nmon(ind,months,:), 'Monthly-Binned Tendency','m',2,'--'}, ...
% $$$           {VCT(ind,months,:), '$\rho_0 C_p \Theta \partial \mathcal{V}/\partial t$','y',2,'-'}, ...
% $$$           {dHdt(ind,months,:), '$\partial \mathcal{H}/\partial t$','c',2,'-'}, ...
% $$$           {Nsnap(ind,:months,:), 'Calculated Tendency',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(ind,:months,:), 'Calculated Tendency 2',0.3*[1 1 1],2,'--'}, ...
          };

leg = {};
legh = [];
for i=1:length(fields)
    hold on;
% $$$     for j=1:length(P(1,1,:))
% $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$              'color',0.7*[1 1 1] ...
% $$$              ,'linewidth',0.5);
% $$$     end
    legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
         'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    hold on;
    leg{i} = fields{i}{2};
end
xlabel('Month');
ylabel('PW');
lg = legend(legh,leg);
ylim([-5 5]);
xlim([1 tL]);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'xtick',[1:tL]);
set(gca,'FontSize',25);
grid on;box on;
LabelAxes(gca,1,25,0.003,0.925);

subplot(2,1,2);
fields = {
% $$$           {F(ind,months,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
          {M(ind,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {I(ind,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {R(ind,months,:), 'Redi Mixing $\mathcal{I}$',[0 0.75 0.75],2,'-'}, ...
% $$$           {N(ind,months,:), 'Tendency $\mathcal{N}$','m',2,'-'}, ...
          {M(ind,months,:)+I(ind,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {Mvc(ind,months,:), 'Vertical Mixing VC term $\mathcal{M}$','r',2,'-'}, ...
% $$$           {P(ind,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {S(ind,months,:), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {F(ind,months,:)+P(ind,months,:), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$           {Nmon(ind,months,:), 'Monthly-Binned Tendency','m',2,'--'}, ...
% $$$           {Nsnap(ind,months,:), 'Calculated Tendency',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(ind,months,:), 'Calculated Tendency 2',0.3*[1 1 1],2,'--'}, ...
          };
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
% $$$     for j=1:length(P(1,1,:))
% $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$              'color',0.7*[1 1 1] ...
% $$$              ,'linewidth',0.5);
% $$$     end
    legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
         'color',fields{i}{3} ...
         ,'linewidth',fields{i}{4});
    hold on;
    leg{i} = fields{i}{2};
end
xlabel('Month');
ylabel('PW');
lg = legend(legh,leg);
ylim([-2 0]);
xlim([1 tL]);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'xtick',[1:tL]);
set(gca,'FontSize',25);
grid on;box on;
LabelAxes(gca,2,25,0.003,0.925);

%%% Spatial Structure:

% Do vertical mixing:
% $$$ VAR = 'FlM';
% $$$ VAR = 'FlSP';
% $$$ VAR = 'WMTP';
VAR = 'WMTSP';
% $$$ TYPE = 'VertInt';
TYPE = 'WMT';
Tl = 22.25;
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
FlM = -FlM;

% $$$ %%% Regional time series 
% $$$ 
% $$$ months = [1:12];
% $$$ %Region choice:
% $$$ LATsplit = 10;
% $$$ 
% $$$ reg = [-150 -90 -5 5];
% $$$ EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);
% $$$ 
% $$$ regA = [-25 5 -5 5];
% $$$ EEAinds = lat > regA(3) & lat < regA(4) & lon > regA(1) & lon < regA(2);
% $$$ 
% $$$ regI = [40 60 -2 8];
% $$$ WEIinds = lat > regI(3) & lat < regI(4) & lon > regI(1) & lon < regI(2);
% $$$ 
% $$$ Eqinds = abs(lat)<=LATsplit;
% $$$ Ninds = lat>LATsplit;
% $$$ Sinds = lat<-LATsplit;
% $$$ 
% $$$ for i=1:12
% $$$     Mtmp = FlM(:,:,i);
% $$$     MEq(i) = nansum(nansum(area(Eqinds).*Mtmp(Eqinds),1),2);
% $$$     MN(i) = nansum(nansum(area(Ninds).*Mtmp(Ninds),1),2);
% $$$     MS(i) = nansum(nansum(area(Sinds).*Mtmp(Sinds),1),2);
% $$$     MEEP(i) = nansum(nansum(area(EEPinds).*Mtmp(EEPinds),1),2);
% $$$     MEEA(i) = nansum(nansum(area(EEAinds).*Mtmp(EEAinds),1),2);
% $$$     MWEI(i) = nansum(nansum(area(WEIinds).*Mtmp(WEIinds),1),2);
% $$$     
% $$$     Atmp = area;
% $$$     Atmp(isnan(Mtmp)) = NaN;
% $$$     AREAtotal(i) = nansum(nansum(Atmp));
% $$$     AREAEEP(i) = nansum(nansum(Atmp(EEPinds)));
% $$$     AREAEEA(i) = nansum(nansum(Atmp(EEAinds)));
% $$$     AREAWEI(i) = nansum(nansum(Atmp(WEIinds)));
% $$$     AREAEq(i) = nansum(nansum(Atmp(Eqinds)));
% $$$     AREAN(i) = nansum(nansum(Atmp(Ninds)));
% $$$     AREAS(i) = nansum(nansum(Atmp(Sinds)));
% $$$     
% $$$ % $$$     Ftmp = FlF(:,:,i);
% $$$ % $$$     FEq(i) = nansum(nansum(area(Eqinds).*Ftmp(Eqinds),1),2);
% $$$ % $$$     FN(i) = nansum(nansum(area(Ninds).*Ftmp(Ninds),1),2);
% $$$ % $$$     FS(i) = nansum(nansum(area(Sinds).*Ftmp(Sinds),1),2);
% $$$ % $$$     FEEP(i) = nansum(nansum(area(EEPinds).*Ftmp(EEPinds),1),2);
% $$$ % $$$ 
% $$$ % $$$     Atmp = FlA(:,:,i);
% $$$ % $$$     AEq(i) = nansum(nansum(area(Eqinds).*Atmp(Eqinds),1),2);
% $$$ % $$$     AN(i) = nansum(nansum(area(Ninds).*Atmp(Ninds),1),2);
% $$$ % $$$     AS(i) = nansum(nansum(area(Sinds).*Atmp(Sinds),1),2);
% $$$ % $$$     AEEP(i) = nansum(nansum(area(EEPinds).*Atmp(EEPinds),1),2);
% $$$ end
% $$$ MAll = MEq+MN+MS;
% $$$ % $$$ FAll = FEq+FN+FS;AAll = AEq+AN+AS;
% $$$ 
% $$$ %Display output in terminal for table:
% $$$ str = {['Area fractions model ' model ' Temp ' num2str(Tl)] ; ...
% $$$        sprintf(' NH area = %3.2f',monmean(AREAN,2,ndays(months))/monmean(AREAtotal,2,ndays(months))) ; ...
% $$$        sprintf(' SH area = %3.2f',monmean(AREAS,2,ndays(months))/monmean(AREAtotal,2,ndays(months))) ; ...
% $$$        sprintf(' Eq area = %3.2f',monmean(AREAEq,2,ndays(months))/monmean(AREAtotal,2,ndays(months))) ; ...
% $$$        sprintf(' EEP area = %3.2f',monmean(AREAEEP,2,ndays(months))/monmean(AREAtotal,2,ndays(months))); ...
% $$$        sprintf(' WEI area = %3.2f',monmean(AREAWEI,2,ndays(months))/monmean(AREAtotal,2,ndays(months))); ...
% $$$        sprintf(' EEA area = %3.2f',monmean(AREAEEA,2,ndays(months))/monmean(AREAtotal,2,ndays(months)))}
% $$$ str = {['Annual totals model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$        sprintf(' Total = %3.2f',monmean(MAll,2,ndays(months))/1e15) ; ...
% $$$        sprintf(' NH = %3.2fPW (%3.0f)',monmean(MN,2,ndays(months))/1e15,monmean(MN,2,ndays(months))/monmean(MAll,2,ndays(months))*100) ; ...
% $$$        sprintf(' SH = %3.2fPW (%3.0f)',monmean(MS,2,ndays(months))/1e15,monmean(MS,2,ndays(months))/monmean(MAll,2,ndays(months))*100) ; ...
% $$$        sprintf(' Eq = %3.2fPW (%3.0f)',monmean(MEq,2,ndays(months))/1e15,monmean(MEq,2,ndays(months))/monmean(MAll,2,ndays(months))*100) ; ...
% $$$        sprintf(' EEP = %3.2fPW (%3.0f)',monmean(MEEP,2,ndays(months))/1e15,monmean(MEEP,2,ndays(months))/monmean(MAll,2,ndays(months))*100); ...
% $$$        sprintf(' WEI = %3.2fPW (%3.0f)',monmean(MWEI,2,ndays(months))/1e15,monmean(MWEI,2,ndays(months))/monmean(MAll,2,ndays(months))*100); ...
% $$$        sprintf(' EEA = %3.2fPW (%3.0f)',monmean(MEEA,2,ndays(months))/1e15,monmean(MEEA,2,ndays(months))/monmean(MAll,2,ndays(months))*100)}
% $$$ str = {['SC range model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$        sprintf(' Total = %3.2fPW',(max(MAll)-min(MAll))/1e15) ; ...
% $$$        sprintf(' NH = %3.2fPW (%3.0f)',(max(MN)-min(MN))/1e15,(max(MN)-min(MN))/(max(MAll)-min(MAll))*100) ; ...
% $$$        sprintf(' SH = %3.2fPW (%3.0f)',(max(MS)-min(MS))/1e15,(max(MS)-min(MS))/(max(MAll)-min(MAll))*100) ; ...
% $$$        sprintf(' Eq = %3.2fPW (%3.0f)',(max(MEq)-min(MEq))/1e15,(max(MEq)-min(MEq))/(max(MAll)-min(MAll))*100) ; ...
% $$$        sprintf(' EEP = %3.2fPW (%3.0f)',(max(MEEP)-min(MEEP))/1e15,(max(MEEP)-min(MEEP))/(max(MAll)-min(MAll))*100); ...
% $$$        sprintf(' WEI = %3.2fPW (%3.0f)',(max(MWEI)-min(MWEI))/1e15,(max(MWEI)-min(MWEI))/(max(MAll)-min(MAll))*100); ...
% $$$        sprintf(' EEA = %3.2fPW (%3.0f)',(max(MEEA)-min(MEEA))/1e15,(max(MEEA)-min(MEEA))/(max(MAll)-min(MAll))*100)}
% $$$ 
% $$$ mn1 = 4;mn2 = 7;
% $$$ str = {['SC range Apr-Jul model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$        sprintf(' Total = %3.2fPW',((MAll(mn1))-(MAll(mn2)))/1e15) ; ...
% $$$        sprintf(' NH = %3.2fPW (%3.0f)',((MN(mn1))-(MN(mn2)))/1e15,((MN(mn1))-(MN(mn2)))/((MAll(mn1))-(MAll(mn2)))*100) ; ...
% $$$        sprintf(' SH = %3.2fPW (%3.0f)',((MS(mn1))-(MS(mn2)))/1e15,((MS(mn1))-(MS(mn2)))/((MAll(mn1))-(MAll(mn2)))*100) ; ...
% $$$        sprintf(' Eq = %3.2fPW (%3.0f)',((MEq(mn1))-(MEq(mn2)))/1e15,((MEq(mn1))-(MEq(mn2)))/((MAll(mn1))-(MAll(mn2)))*100) ; ...
% $$$        sprintf(' EEP = %3.2fPW (%3.0f)',((MEEP(mn1))-(MEEP(mn2)))/1e15,((MEEP(mn1))-(MEEP(mn2)))/((MAll(mn1))-(MAll(mn2)))*100); ...
% $$$        sprintf(' WEI = %3.2fPW (%3.0f)',((MWEI(mn1))-(MWEI(mn2)))/1e15,((MWEI(mn1))-(MWEI(mn2)))/((MAll(mn1))-(MAll(mn2)))*100); ...
% $$$        sprintf(' EEA = %3.2fPW (%3.0f)',((MEEA(mn1))-(MEEA(mn2)))/1e15,((MEEA(mn1))-(MEEA(mn2)))/((MAll(mn1))-(MAll(mn2)))*100)}
% $$$ 
% $$$ % The 1% (within Nino 3) area count of the annual mean:
% $$$ M1p = zeros(length(find(EEPinds)),12);
% $$$ A1p = area(EEPinds);
% $$$ for i=1:12
% $$$     Mtmp = FlM(:,:,i);
% $$$     M1p(:,i) = A1p.*Mtmp(EEPinds);
% $$$ end
% $$$ M1pA = monmean(M1p,2,ndays(months));
% $$$ [M1pA,I] = sort(M1pA);
% $$$ M1p = M1p(I,:);
% $$$ Atotal = monmean(AREAtotal,2,ndays(months));
% $$$ [tmp ind] = min(abs(cumsum(A1p(I))/Atotal - 0.005));
% $$$ %plot(cumsum(A1p(I))/Atotal*100,cumsum(M1pA)/monmean(MAll,2,ndays(months))*100)
% $$$ M1pA = sum(M1pA(1:ind));
% $$$ M1pSCR = (sum(M1p(1:ind,mn1))-sum(M1p(1:ind,mn2)));
% $$$ 
% $$$ str = {['1% area within Nino 3 model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$        sprintf(' 1p Annual-Mean = %3.2fPW (%3.0f)',M1pA/1e15,M1pA/monmean(MAll,2,ndays(months))*100); ...
% $$$        sprintf(' 1p SC Apr-Jul = %3.2fPW (%3.0f)',M1pSCR/1e15,M1pSCR/((MAll(mn1))-(MAll(mn2)))*100)}
% $$$ 
% $$$ % $$$ % Average fluxes, isotherm separationl, diffusivity:
% $$$ % $$$ load([base 'MOM025_output002_Ziso_T23C.mat']);
% $$$ % $$$ topiso = ziso;
% $$$ % $$$ load([base 'MOM025_output002_Ziso_T22C.mat']);
% $$$ % $$$ botiso = ziso;
% $$$ % $$$ AvgFlux = MAll./AREAtotal;
% $$$ % $$$ NaNs = isnan(topiso-botiso);
% $$$ % $$$ AvgDZ = zeros(12,1);
% $$$ % $$$ for ii=1:12
% $$$ % $$$     AvgDZ(ii) = nansum(nansum(((topiso(:,:,ii)-botiso(:,:,ii)).*area),1),2)./ ...
% $$$ % $$$                 nansum(area(~isnan((topiso(:,:,ii)-botiso(:,:,ii)))));
% $$$ % $$$ end
% $$$ % $$$ AvgDiff = monmean(AvgFlux/rho0/Cp.*AvgDZ',2,ndays(months))
% $$$ % $$$ AvgFlux = monmean(AvgFlux,2,ndays(months))
% $$$ % $$$ AvgDZ = monmean(AvgDZ',2,ndays(months))
% $$$ % $$$ sprintf('Average flux across the isotherm = %3.1f Wm-2',AvgFlux)
% $$$ % $$$ sprintf('Average diffusivity = %3.5f m2s-1',AvgDiff)
% $$$ 
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ 
% $$$ % $$$ subplot(3,1,1);
% $$$ plot(1:12,MEEP/1e15,'--r','LineWidth',4);
% $$$ hold on;
% $$$ plot(1:12,MEq/1e15,'--k','LineWidth',4);
% $$$ plot(1:12,MN/1e15,':b','LineWidth',4);
% $$$ plot(1:12,MS/1e15,':','color',[0 0.5 0],'LineWidth',4);
% $$$ plot(1:12,(MN+MS)/1e15,':k','LineWidth',4);
% $$$ plot(1:12,MAll/1e15,'-k','LineWidth',4);
% $$$ xlabel('Month');
% $$$ ylabel(['PW']);
% $$$ title(['Vertical mixing heat flux through $' num2str(Tl) '^\circ$C isotherm']);
% $$$ leg = legend('Eastern Equatorial Pacific', ...
% $$$              'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
% $$$              ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
% $$$              ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
% $$$              ['Total']);
% $$$ set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
% $$$ ylim([-0.8 0]);
% $$$ xlim([1 12]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'xtick',[1:12]);
% $$$ set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
% $$$ set(gca,'FontSize',25);
% $$$ grid on;
% $$$ 
% $$$ axes('Position',[0.74 0.2451 0.13 0.6799]);
% $$$ plot([0 1],[min(MAll/1e15) min(MAll/1e15)],'-k','linewidth',2);
% $$$ hold on;
% $$$ plot([0 1],[max(MAll/1e15) max(MAll/1e15)],'-k','linewidth',2);
% $$$ text(0.5,mean([min(MAll/1e15) max(MAll/1e15)]),sprintf('%3.2f ',max(MAll/1e15)-min(MAll/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ plot([1 2],[min(MEq/1e15) min(MEq/1e15)],'--k','linewidth',2);
% $$$ plot([1 2],[max(MEq/1e15) max(MEq/1e15)],'--k','linewidth',2);
% $$$ text(1.5,mean([min(MEq/1e15) max(MEq/1e15)]),sprintf('%3.2f ',max(MEq/1e15)-min(MEq/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ plot([2 3],[min(MEEP/1e15) min(MEEP/1e15)],'--r','linewidth',2);
% $$$ plot([2 3],[max(MEEP/1e15) max(MEEP/1e15)],'--r','linewidth',2);
% $$$ text(2.5,mean([min(MEEP/1e15) max(MEEP/1e15)]),sprintf('%3.2f ',max(MEEP/1e15)-min(MEEP/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','r','FontSize',25);
% $$$ plot([3 4],[min(MN/1e15) min(MN/1e15)],':b','linewidth',2);
% $$$ plot([3 4],[max(MN/1e15) max(MN/1e15)],':b','linewidth',2);
% $$$ text(3.5,mean([min(MN/1e15) max(MN/1e15)]),sprintf('%3.2f ',max(MN/1e15)-min(MN/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','b','FontSize',25);
% $$$ plot([4 5],[min(MS/1e15) min(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$                     0.5 0]);
% $$$ plot([4 5],[max(MS/1e15) max(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$                     0.5 0]);
% $$$ text(4.5,mean([min(MS/1e15) max(MS/1e15)]),sprintf('%3.2f ',max(MS/1e15)- min(MS/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color',[0 0.5 0],'FontSize',25);
% $$$ plot([5 6],[min((MN+MS)/1e15) min((MN+MS)/1e15)],':k','linewidth',2);
% $$$ plot([5 6],[max((MN+MS)/1e15) max((MN+MS)/1e15)],':k','linewidth',2);
% $$$ text(5.5,mean([min((MN+MS)/1e15) max((MN+MS)/1e15)]),sprintf('%3.2f ',max((MN+MS)/1e15)-min((MN+MS)/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ xlim([0 6]);
% $$$ ylim([-0.8 0]);
% $$$ set(gca,'xtick',[]);
% $$$ set(gca,'ytick',[]);
% $$$ box off;
% $$$ grid off;
% $$$ axis off;
% $$$ title('Seasonal Range','FontSize',25);

%%% Plot spatial pattern:

obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
LAND = obj.SST(:,:,1);

%If MOM01, fix NaN's in grid:
if (strfind(model,'01'))
    lon = repmat(lon(:,500),[1 yL]);
    latv = nanmean(lat,1);
    lat = repmat(latv,[xL 1]);
end

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = {[1:12], ...
          [3], ...
          [7], ...
          [11]};
labels = {'(a) Annual', ...
          '(b) March', ...
          '(c) July', ...
          '(d) November'};

%Colormap and continents:
% $$$ clim = [-100 0];
% $$$ sp = 10;
% $$$ clim = [-30 0]; % FOR SWP
% $$$ sp = 0.5; % FOR SWP
clim = [-1 1]*(1e-5)*86400; % FOR WMT
clim = [-0.3 0.3]*(1e-5)*86400; % FOR WMT
sp = (0.3e-6)*86400;

doWMT = 1; % plot WMT instead of flux

cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

if (doWMT)
    cmap = flipud(lbmap(npts-3,'RedBlue'));
else
    cmap = parula(npts-3);
    cmap = parula(npts-3);
    cmap(end,:) = [0.97 0.97 0.8];
    cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
end

tmp = LAND;
tmp(isnan(LAND)) = clim(1)-sp/2;
tmp(~isnan(LAND)) = NaN;
LAND = tmp;
cmap(2:(end+1),:) = cmap;
cmap(1,:) = [0 0 0];

climn = [clim(1)-sp clim(2)];

%Mean of all months:
figure;
set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',20);
poss = [0.1300    0.4553    0.7693    0.4697; ...
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
    if (length(months{i})>1)
        tmp = FlM;
        tmp(isnan(tmp)) = 0.0;
        Z = monmean(tmp(:,:,months{i}),3,ndays(months{i}));
        Z(Z == 0) = NaN;
    else
        Z = FlM(:,:,months{i});
    end
    Z = Z(xvec,yvec);
    if (doWMT)
        Z = Z*86400;
    end
    
    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z.*cos(Y/180*pi),cpts,'linestyle','none');
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    if (i==1)
        cb = colorbar;
        if (~doWMT)
            ylabel(cb,'Wm$^{-2}$');
        else
            ylabel(cb,'m/day');
        end            
        ylim(cb,clim);
    end
    hold on;
% $$$     xlims = get(gca,'xlim');
% $$$     plot(xlims,LATsplit*[1 1],'--k','linewidth',1);
% $$$     plot(xlims,-LATsplit*[1 1],'--k','linewidth',1);
% $$$     plot([reg(1:2) reg(2:-1:1) reg(1)],[reg(3) reg(3:4) reg(4:-1:3)],'--k','linewidth',1);
% $$$     plot([regA(1:2) regA(2:-1:1) regA(1)],[regA(3) regA(3:4) regA(4:-1:3)],'--k');
% $$$     plot([regI(1:2) regI(2:-1:1) regI(1)],[regI(3) regI(3:4) regI(4:-1:3)],'--k');
    if (i>1)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i<=2)
        ylabel('Latitude ($^\circ$N)');
    end
    
    if (i>1)
        text(77,37.5,labels{i},'BackgroundColor','w','Margin',0.5,'HorizontalAlignment','right');
        set(gca,'xtick',[-240:60:60]);
    else
        text(-279,41,labels{i},'BackgroundColor','w','Margin',0.5);
        set(gca,'xtick',[-270:30:60]);
    end        
    set(gca,'Position',[poss(i,:)]);
    ylim([-45 45]);
    set(gca,'ytick',[-45:15:45]);
% $$$     ylim([-60 60]);
% $$$     set(gca,'ytick',[-75:15:75]);
end 
colormap(cmap);
%colormap(parula);%flipud(lbmap(50,'RedBlue')));

%%% Plot spatial pattern of net heat flux and SST:

% Load Variable and calculate mean:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
shfluxa = shflux;
SSTa = SST;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    shfluxa = shfluxa+shflux;
    SSTa = SSTa+SST;
end
shflux = shfluxa/length(outputs);
SST = SSTa/length(outputs);

%If MOM01, fix NaN's in grid:
if (strfind(model,'01'))
    lon = repmat(lon(:,500),[1 yL]);
    latv = nanmean(lat,1);
    lat = repmat(latv,[xL 1]);
end

%Sum of all positives:
shfluxA = mean(shflux,3);
shfluxP = nansum(shfluxA(shfluxA>0).*area(shfluxA>0))/1e15
shfluxM = nansum(shfluxA(shfluxA<0).*area(shfluxA<0))/1e15

%Sum of mean flux above mean isotherm:
SSTA = mean(SST,3);
isot = 23;
shfluxP = nansum(shfluxA(SSTA>=isot).*area(SSTA>=isot))/1e15
shfluxM = nansum(shfluxA(SSTA<isot).*area(SSTA<isot))/1e15

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = {[1:12]};%, ...
                  %          [3], ...
                  %          [7], ...
                  %          [11]};

labels = {'(a) Annual', ...
          '(b) March', ...
          '(c) July', ...
          '(d) November'};

figure;
set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',20);
poss = [0.1300    0.415  0.7693    0.56; ...
        0.1300    0.1    0.2343    0.2680; ...
        0.3951    0.1    0.2343    0.2680; ...
        0.6681    0.1    0.2343    0.2680];
for i=1:length(months)
    if (i == 1)
        subplot(5,3,[1 9]);
    else
        subplot(5,3,[10 13]+(i-2));
    end
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
    Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
    Z = Z(xvec,yvec);
    Z2 = Z2(xvec,yvec);
    
    % Map projection:
% $$$     axesm('MapProjection','robinson','MapLatLimit',[-78 67], ...
% $$$       'Origin',[0 205 0],'Frame','off','Grid','on', ...
% $$$       'MeridianLabel','on','ParallelLabel','on', ...
% $$$       'MLabelLocation',90,'FontName','Latex','FontSize',15)
% $$$     axesm('MapProjection','robinson');
% $$$     [X,Y] = mfwdtran(X,Y);
    contourf(X,Y,Z.*cos(Y/180*pi),[-1e10 -500:20:500 1e10],'linestyle','none');
% $$$     contourf(X,Y,Z,[-1e10 -500:20:500 1e10],'linestyle','none');
    hold on;
% $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$     if (i==1)
    [c,h] = contour(X,Y,Z2,[-3:2:21 25:2:35],'-k');
    clabel(c,h);        
    if (i == 1)
        [c,h] = contour(X,Y,Z2,[23 23],'-k','linewidth',2);
        clabel(c,h);
    else
    end
    caxis([-200 200]);
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
    end
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
    end
    ylim([-75 75]);
    if (i>1)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i<=2)
        ylabel('Latitude ($^\circ$N)');
    end
    if (i>1)
% $$$         text(-277,53,labels{i},'BackgroundColor','w');
        text(-277,64,labels{i},'BackgroundColor','w','margin',0.01);
        set(gca,'xtick',[-240:60:60]);
    else
% $$$         text(-279,55,labels{i},'BackgroundColor','w');
        text(-279,69,labels{i},'BackgroundColor','w','margin',0.01);
        set(gca,'xtick',[-270:30:60]);
    end        
    set(gca,'ytick',[-60:30:60]);
    set(gca,'Position',[poss(i,:)]);
    set(gca,'color','k');
end 
colormap(redblue(20));

%%%%%%%%%%%% Supp FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot SST difference and shflux difference plots:

% Load MOM01 data:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/';
model = 'MOM01';
outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

[xL,yL] = size(lon)
%If MOM01, fix NaN's in grid:
if (strfind(model,'01'))
    lonv = lon(:,500);
    latv = nanmean(lat,1);
    [tmp ind] = min(abs(latv - 60))
    [lon,lat] = ndgrid(lonv,latv(1:ind));
end
MOM01lon = lon;
MOM01lat = lat;

load([base model sprintf('_output%03d_SurfaceVars.mat', ...
                         outputs(1))]);
MOM01shflux = shflux(:,1:ind,:);
MOM01SST = SST(:,1:ind,:);

% Load MOM025:
model = 'MOM025';
outputs = [2 3 4 5 6];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
shfluxa = shflux;
SSTa = SST;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    shfluxa = shfluxa+shflux;
    SSTa = SSTa+SST;
end
shflux = shfluxa/length(outputs);
SST = SSTa/length(outputs);

% Calculate bias:
for i=1:12
    SST(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01SST(:,:,i)',lon,lat,'linear')-SST(:,:,i);
    shflux(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01shflux(:,:,i)',lon,lat,'linear')-shflux(:,:,i);
end

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = {[1:12], ...
          [3], ...
          [7], ...
          [11]};

labels = {'Annual', ...
          'March', ...
          'July', ...
          'November'};

figure;
set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
poss = [0.1300    0.4553    0.7693    0.4697; ...
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
    Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
    Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
    Z = Z(xvec,yvec);
    Z2 = Z2(xvec,yvec);
% $$$     contourf(X,Y,Z2,[-1e10 -5:0.25:5 1e10],'linestyle','none');
    contourf(X,Y,Z,[-1e10 -500:10:500 1e10],'linestyle','none');
    hold on;
% $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$     if (i==1)
% $$$         [c,h] = contour(X,Y,Z,[-100:10:-10],'--k');
% $$$         [c,h] = contour(X,Y,Z,[10:10:100],'-k');
% $$$     end
% $$$     else
% $$$         [c,h] = contour(X,Y,Z2,[-3:4:35],'-k');
% $$$     end
% $$$     clabel(c,h);
    caxis([-100 100]);
    if (i==1)
        cb = colorbar;
% $$$         ylabel(cb,'$^\circ$C');
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

%%% Outcrop area plot:

Aout = zeros(size(Te));
for i=1:length(Te)
    i
    Aout(i) = nansum(area(mean(SST,3)>Te(i)));
end
figure;
plot(Aout/1e6,Te,'-k','linewidth',2);
xlabel('Outcrop area (km$^2$)');
ylabel('Temperature ($^\circ$C)');
plot(-diff(Aout,[],1)/dT/1e6,T,'-k','linewidth',2);
xlabel('Outcrop area (km$^2$/$^\circ$C)');
ylabel('Temperature ($^\circ$C)');

%%% M-L T-t plots: 
months = 1:12;
fields = {
% $$$           {F(:,months,:)+P(:,months,:), 'Surface Forcing $\mathcal{F}+\mathcal{P}$','k',2,'-'}, ...
          {M(:,months,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {I(:,months,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {N(:,months,:), 'Tendency $\mathcal{N}$','m',2,'-'}, ...
% $$$           {M(:,months,:)+I(:,months,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {P(:,months,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {S(:,months,:), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {F(:,months,:)+P(:,months,:), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$           {Nmon(:,months,:), 'Monthly-Binned Tendency','m',2,'--'}, ...
% $$$           {Nsnap(:,months,:), 'Calculated Tendency',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(:,months,:), 'Calculated Tendency 2',0.3*[1 1 1],2,'--'}, ...
          };

% Fluxes:
scale = 1/1e15;label = '(PW)';
caxs = [-15 15];x = Te;

figure;
%set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'Position',[272         107        1302         863]);
for ii=1:length(fields)
    subplot(1,length(fields),ii);
    V = mean(fields{ii}{1},3)'*scale;
    [X,Y] = ndgrid(1:tL,Te);
    cint = [-1e10 -10:0.05:10 1e10];
    contourf(X,Y,V,cint);%,'linestyle','none');
    hold on;
    plot([1 tL],[22.5 22.5],'--k','linewidth',2);
    cb = colorbar('Location','NorthOutside','FontSize',25);    
% $$$     set(gca,'xtick',0.5:1:11.5);
% $$$     set(gca,'xticklabel',[1:12]);
    set(gca,'ytick',-5:5:35);
    set(gca,'xtick',[1:tL]);
    ylim([-3 31]);
    grid on;
    %    caxis(caxs);
    caxis([-1 0]);
    xlabel('Month');
    ylabel('Temperature ($^\circ$C)');
    %    xlabel(cb,['MOM025 Global ' fields{ii}{2} ' '
    %    label],'FontSize',20);
    xlabel(cb,[fields{ii}{2} ' (PW)'],'FontSize',25);
    set(gca,'FontSize',25);
end
% $$$ colormap(redblue);
% $$$ caxis([-2 2]);
cmap = redblue((length(cint)-3)*2);
cmap = cmap(1:(length(cint)-3),:);
colormap(cmap);
% $$$ caxis([-5 5]);
hold on;
plot([1 tL],[22.5 22.5],'--k','linewidth',2);

%%% Spatial Integrated Mixing

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [2 3 4 5 6];
% $$$ outputs = [2];
% $$$ model = 'MOM01';
% $$$ outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

%Region choice:
LATsplit = 10;
reg = [-150 -90 -5 5];
EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);

obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
LAND = obj.SST(:,:,1);

%If MOM01, fix NaN's in grid:
if (strfind(model,'01'))
    lon = repmat(lon(:,500),[1 yL]);
    latv = nanmean(lat,1);
    lat = repmat(latv,[xL 1]);
end

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = [1:12];

%Colormap and continents:
clim = [-50 0];
sp = 2.5;
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

% $$$ cmap = flipud(lbmap(2*(npts-3),'RedBlue'));
% $$$ cmap = cmap(1:(npts-3),:);
% $$$ cmap = redblue(2*(npts-3));
% $$$ cmap = cmap(1:(npts-3),:);
%cmap = summer(npts-3);
cmap = parula(npts-3);

%Custom parula:
cmap = parula(npts-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ mid = round((npts-3)/2);
% $$$ cmap(2:(mid-1),1) = interp1([1 mid],[cmap(1,1) cmap(mid,1)],2:(mid-1),'linear');
% $$$ cmap(2:(mid-1),2) = interp1([1 mid],[cmap(1,2) cmap(mid,2)],2:(mid-1),'linear');
% $$$ cmap(2:(mid-1),3) = interp1([1 mid],[cmap(1,3) cmap(mid,3)],2:(mid-1),'linear');
% $$$ cmap((mid+1):(end-1),1) = interp1([mid npts-3],[cmap(mid,1) cmap(end,1)],(mid+1):(npts-4),'linear');
% $$$ cmap((mid+1):(end-1),2) = interp1([mid npts-3],[cmap(mid,2) cmap(end,2)],(mid+1):(npts-4),'linear');
% $$$ cmap((mid+1):(end-1),3) = interp1([mid npts-3],[cmap(mid,3) cmap(end,3)],(mid+1):(npts-4),'linear');


% $$$ %Custom:
% $$$ cmap = zeros(npts-3,3);
% $$$ cmap(1,:) = [0 0 1];
% $$$ cmap(end,:) = [0.9359 0.9323 0.7947];
% $$$ %cmap(end,:) = [1 1 0.85];
% $$$ %cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(2:(end-1),1) = interp1([1 npts-3],[cmap(1,1) cmap(end,1)],2:(npts-4),'linear');
% $$$ cmap(2:(end-1),2) = interp1([1 npts-3],[cmap(1,2) cmap(end,2)],2:(npts-4),'linear');
% $$$ cmap(2:(end-1),3) = interp1([1 npts-3],[cmap(1,3) cmap(end,3)],2:(npts-4),'linear');

tmp = LAND;
tmp(isnan(LAND)) = clim(1)-sp/2;
tmp(~isnan(LAND)) = NaN;
LAND = tmp;
cmap(2:(end+1),:) = cmap;
cmap(1,:) = [0 0 0];

climn = [clim(1)-sp clim(2)];

%Mean of all months:
figure;
set(gcf,'Position',[3    40   956   963]);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
poss = [0.100    0.6993    0.7750    0.2837; ...
        0.100    0.3996    0.7750    0.2837; ...
        0.100    0.100    0.7750    0.2837];

Tls = [10 17.5 27.5]
labels = {'(a) 10$^\circ$C','(b) 17.5$^\circ$C','(c) 27.5$^\circ$C'};

% $$$ Tls = [-1.5 -1 2.5]
% $$$ labels = {'-1.5$^\circ$C','-1$^\circ$C','2.5$^\circ$C'};

for i=1:length(Tls)
    Tl = Tls(i)
    % Load Variable and calculate mean:
    load([base model sprintf('_output%03d_VertInt_T',outputs(1)) strrep(num2str(Tl),'.','p') 'C.mat']);
    if (strfind(model,'025'))
        FlM(isnan(FlM)) = 0.0;
        FlMa = FlM;
        for ii=2:length(outputs)
            load([base model sprintf('_output%03d_VertInt_T',outputs(ii)) strrep(num2str(Tl),'.','p') 'C.mat']);
            FlM(isnan(FlM)) = 0.0;
            FlMa = FlMa+FlM;
        end
        FlM = FlMa/length(outputs);
        FlM(FlM==0) = NaN;
    end

    subplot(3,1,i);
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    if (length(months)>1)
        tmp = FlM;
        tmp(isnan(tmp)) = 0.0;
        Z = monmean(tmp(:,:,months),3,ndays(months));
        Z(Z == 0) = NaN;
    else
        Z = FlM(:,:,months);
    end
    Z = Z(xvec,yvec);
    
    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z.*cos(Y/180*pi),cpts,'linestyle','none');
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
    ylim(cb,clim);
    hold on;
    xlims = get(gca,'xlim');
% $$$     plot(xlims,LATsplit*[1 1],'--k');
% $$$     plot(xlims,-LATsplit*[1 1],'--k');
% $$$     plot([reg(1:2) reg(2:-1:1) reg(1)],[reg(3) reg(3:4) reg(4:-1:3)],'--k');
    ylim([-75 75]);
    if (i==3)
        xlabel('Longitude ($^\circ$E)');
    end
    ylabel('Latitude ($^\circ$N)');
    text(-276,65,labels{i},'BackgroundColor','w','Margin',0.5);
    set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'ytick',[-45:15:45]);
    set(gca,'xtick',[-360:60:360]);
    if (i~=3)
        set(gca,'xticklabel',[]);
    end
end 
colormap(cmap);
%colormap(parula);%flipud(lbmap(50,'RedBlue')));


%%% Inferred meridional heat flux:

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/';
% $$$ model = 'MOM025';
% $$$ outputs = [2 3 4 5 6];
model = 'MOM01';
outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

% Load Variable and calculate mean:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
shfluxa = shflux;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    shfluxa = shfluxa+shflux;
end
shflux = shfluxa/length(outputs);
shflux = monmean(shflux,3,ndays);

% Calculate meridional heat flux inferred:
latV = linspace(-90,90,181);
V = zeros(size(latV));
for i=1:length(latV)
    inds = lat < latV(i);
    V(i) = nansum(area(inds).*shflux(inds));
end

%Center the flux:
V = V + (V(1)-V(end))/2;

figure;
set(gcf,'Position',[260         339        1055         586]);
set(gcf,'defaulttextfontsize',25);
set(gcf,'defaultaxesfontsize',25);
plot(latV,V/1e15,'-r','linewidth',2);
xlabel('Latitude ($^\circ$N)');
ylabel('Meridional Heat Flux (PW)');
grid on;
box on;
xlim([-90 90]);
ylim([-1 2]);
set(gca,'xtick',[-90:30:90]);

%Calculate fractions of surface heat flux regionally:

%Region choice:
LATsplit = 10;
reg = [-150 -90 -5 5];
EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);

Eqinds = abs(lat)<=LATsplit;
Ninds = lat>LATsplit;
Sinds = lat<-LATsplit;

for i=1:12
    SHFtmp = shflux(:,:,i);
    SHFEq(i) = nansum(nansum(area(Eqinds).*SHFtmp(Eqinds),1),2);
    SHFN(i) = nansum(nansum(area(Ninds).*SHFtmp(Ninds),1),2);
    SHFS(i) = nansum(nansum(area(Sinds).*SHFtmp(Sinds),1),2);
    SHFEEP(i) = nansum(nansum(area(EEPinds).*SHFtmp(EEPinds),1),2);

    Atmp = area;
    Atmp(isnan(SHFtmp)) = NaN;
    AREAtotal(i) = nansum(nansum(Atmp));
    AREAEEP(i) = nansum(nansum(Atmp(EEPinds)));
    AREAEq(i) = nansum(nansum(Atmp(Eqinds)));
    AREAN(i) = nansum(nansum(Atmp(Ninds)));
    AREAS(i) = nansum(nansum(Atmp(Sinds)));
end
SHFAll = SHFEq+SHFN+SHFS;

%Display output in terminal for table:
str = {['Area fractions model ' model] ; ...
sprintf(' NH area = %3.2f',monmean(AREAN,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
sprintf(' SH area = %3.2f',monmean(AREAS,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
sprintf(' Eq area = %3.2f',monmean(AREAEq,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
sprintf(' EEP area = %3.2f',monmean(AREAEEP,2,ndays)/monmean(AREAtotal,2,ndays))}
str = {['Annual totals model ' model]  ; ...
sprintf(' Total = %3.2f',monmean(SHFAll)/1e15) ; ...
sprintf(' NH = %3.2fPW (%3.0f)',monmean(SHFN,2,ndays)/1e15,monmean(SHFN,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
sprintf(' SH = %3.2fPW (%3.0f)',monmean(SHFS,2,ndays)/1e15,monmean(SHFS,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
sprintf(' Eq = %3.2fPW (%3.0f)',monmean(SHFEq,2,ndays)/1e15,monmean(SHFEq,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
sprintf(' EEP = %3.2fPW (%3.0f)',monmean(SHFEEP,2,ndays)/1e15,monmean(SHFEEP,2,ndays)/monmean(SHFAll,2,ndays)*100)}

% Plot regional heat fluxes as a function of season:
figure;
set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'defaultlinelinewidth',2);

% $$$ subplot(3,1,1);
plot(1:12,SHFEEP/1e15,'--r','LineWidth',4);
hold on;
plot(1:12,SHFEq/1e15,'--k','LineWidth',4);
plot(1:12,SHFN/1e15,':b','LineWidth',4);
plot(1:12,SHFS/1e15,':','color',[0 0.5 0],'LineWidth',4);
plot(1:12,(SHFN+SHFS)/1e15,':k','LineWidth',4);
plot(1:12,SHFAll/1e15,'-k','LineWidth',4);
xlabel('Month');
ylabel(['PW']);
title(['Surface Heat Flux']);
leg = legend('Eastern Equatorial Pacific', ...
             'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
             ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
             ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
             ['Total']);
set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
ylim([-20 20]);
xlim([1 12]);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
set(gca,'xtick',[1:12]);
set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
set(gca,'FontSize',25);
grid on;


%%% Latitudinal Slices:

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [2 3 4 5 6];
% $$$ model = 'MOM01';
% $$$ outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

% Load Variable and calculate mean:
lonsl = 110;
load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(1))]);
vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd'};
for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
end
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(i))]);
    for i=1:length(vars)
        eval([vars{i} 'a = ' vars{i} 'a + ' vars{i} ';']);
    end
end
for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a/length(outputs);']);
    eval(['clear ' vars{i} 'a;']);
end

[yL,zL,tL] = size(temp);
TL = length(T);

% Depth of isotherms:
Zi = zeros(yL,TL,tL);
for ti=1:tL
    for yi=1:yL
        tvec = squeeze(temp(yi,:,ti));
        zvec = -Zt(yi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(yi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(yi,:,ti)),1,'last');
        Zi(yi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Yi = repmat(Yt(:,1),[1 TL]);

% $$$ var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux
% $$$ clim = [-250 0];
% $$$ sp = 5;
% $$$ doWMT = 0;

var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
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

% $$$ %Save for schematic:
% $$$ Xl = Yi;
% $$$ Yl = nanmonmean(Zi(:,:,months{i}),3,ndays(months{i}));
% $$$ Zl = nanmonmean(vdif(:,:,months{i}),3,ndays(months{i}));
% $$$ XlC = Yt;
% $$$ YlC = -Zt;
% $$$ ZlC = monmean(temp(:,:,months{i}),3,ndays(months{i}));

% $$$ [tmp Eqind] = min(abs(Yu(:,1)));
% $$$ tauweight = abs(taux(Eqind,:))*200;
% $$$ % Wind stress vectors:
% $$$ sp  =5;
% $$$ yvec = Yu(:,1);
figure;
set(gcf,'Position',[1          36        1920         970]);
for i=1:length(months)
subplot(2,2,i);
contourf(Yi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
[c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
clabel(c,h,[0:2:35]);
[c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3,ndays(months{i})),[23 23],'-k','linewidth',2);
if (strcmp(model,'MOM01'))
    mnu = monthsu01{i};
else
    mnu = months{i};
end
ucol = [0.8706    0.4902         0];
[c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
                'color',ucol);
[c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
                'color',ucol);
plot(Yu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
ylim([-200 0]);
xlim([-10 10]);
cb = colorbar;
if (doWMT)
    ylabel(cb,'m/day');
else
    ylabel(cb,'Wm$^{-2}$');
end
xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
caxis(clim);
text(-9.6,-188,labels{i},'Backgroundcolor','w','FontSize',20);
text(9.6,-188,[num2str(lonsl) '$^\circ$W'],'Backgroundcolor','w','FontSize',20,'HorizontalAlignment','Right');

% $$$ %Add wind-stress vectors:
% $$$ pos = get(gca,'Position')
% $$$ wsh = axes('Position',[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]);
% $$$ % $$$ for ii=1:sp:length(yvec)
% $$$     plot(0,0,'o','MarkerSize',abs(mean(mean(taux(yvec>=-10 & yvec<=10,months{i}),2),1))*200);
% $$$     hold on;
% $$$ % $$$ end
% $$$ %quiver(yvec,zeros(size(yvec)),mean(taux(1:sp:end,months{i}),2),mean(tauy(1:sp:end,months{i}),2));
% $$$ xlim([-10 10]);
% $$$ ylim([-1 1]);
% $$$ box off;axis off;
% $$$ set(wsh,'Position',[[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]]);

LabelAxes(gca,i,20,0.008,0.95);
end
colormap(cmap);

%%% Equatorial Slices:

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [2 3 4 5 6];
% $$$ model = 'MOM01';
% $$$ outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

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
        zvec = -Z(xi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(xi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(xi,:,ti)),1,'last');
        Zi(xi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Xi = repmat(X(:,1),[1 TL]);

% $$$ var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux
% $$$ clim = [-250 0];
% $$$ sp = 5;
% $$$ doWMT = 0;

% $$$ var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
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

% $$$ %Save for schematic:
% $$$ Xe = Xi;
% $$$ Ye = nanmonmean(Zi(:,:,months{i}),3,ndays(months{i}));
% $$$ Ze = nanmonmean(vdif(:,:,months{i}),3,ndays(months{i}));
% $$$ XeC = X;
% $$$ YeC = -Z;
% $$$ ZeC = monmean(temp(:,:,months{i}),3,ndays(months{i}));

% $$$ [tmp Eqind] = min(abs(Yu(:,1)));
% $$$ tauweight = abs(taux(Eqind,:))*200;
% $$$ % Wind stress vectors:
% $$$ sp  =5;
% $$$ yvec = Yu(:,1);
figure;
set(gcf,'Position',[1          36        1920         970]);
for i=1:length(months)
subplot(2,2,i);
contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
[c,h] = contour(X,-Z,monmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
clabel(c,h,[0:2:35]);
[c,h] = contour(X,-Z,monmean(temp(:,:,months{i}),3,ndays(months{i})),[23 23],'-k','linewidth',2);
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


%%% Temperature-latitude heat function:

% Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [7];
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_nipoall/mat_data/';
model = 'MOM025_nipoall';
outputs = [19];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

% Load streamfunction:
load([base model sprintf('_output%03d_Tpsi.mat',outputs(1))]);

% Load SST:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);

maxSST = squeeze(max(max(SST,[],1),[],3));
meanSST = squeeze(nanmean(nanmean(SST,1),3));
minSST = squeeze(min(min(SST,[],1),[],3));

NaNs = PSI == 0;

%cumsum from low temps and correct units:
PSI = cumsum(-PSI/rho0*1e9,2);

%calculate heat function:
H = cumsum(rho0*Cp*PSI*dT,2);

PSI(NaNs) = NaN;
H(NaNs) = NaN;

[YY,TT] = ndgrid(latv,T);

%monthsc = {[1:12],[3],[7],[11]};
monthsc = {[1:12]};%,[3],[7],[11]};

figure;
set(gcf,'Position',[157         359        1742         586]);
%set(gcf,'Position',[3    40   956   963]);
%set(0,'defaulttextfontsize',10);
%set(0,'defaultaxesfontsize',10);
for i=1:length(monthsc)
    months = monthsc{i}
subplot(length(monthsc),2,2*i-1);
contourf(YY,TT,monmean(PSI(:,:,months)/1e6,3,ndays(months)),[-1e10 -100:5:100 1e10],'-k');
hold on;
contour(YY,TT,monmean(PSI(:,:,months)/1e6,3,ndays(months)),[0 0],'-k','linewidth',2);
plot(latv,maxSST,'--r','linewidth',2);
plot(latv,meanSST,'--k','linewidth',2);
plot(latv,minSST,'--b','linewidth',2);
caxis([-40 40]);
colormap(redblue);
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
title('Transport Streamfunction (Sv)');
colorbar;

subplot(length(monthsc),2,2*i);
contourf(YY,TT,monmean(H(:,:,months)/1e15,3,ndays(months)),[-1e10 -5:0.25:5 1e10],'-k');
hold on;
contour(YY,TT,monmean(H(:,:,months)/1e15,3,ndays(months)),[0 0],'-k','linewidth',2);
plot(latv,maxSST,'--r','linewidth',2);
plot(latv,meanSST,'--k','linewidth',2);
plot(latv,minSST,'--b','linewidth',2);
caxis([-3 3]);
colormap(redblue);
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
title('Heat Function (PW)');
colorbar;
end

% $$$ load([base model sprintf('_output%03d_HFunc.mat',outputs(1))]);

% $$$ TENa = HFETS+HFFRZ+HFKNL+HFPME+HFRMX+HFSUB+HFSWH+HFVDF+HFVDS;
% $$$ latv = max(lat,[],1);


%%% Plot Seasonal cycle of wind stress and SST

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [2 3 4 5 6];
% $$$ model = 'MOM01';
% $$$ outputs = [333];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

% Load Variable and calculate mean:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
SSTa = SST;
tauxa = taux;
tauya = tauy;
for i=2:length(outputs)
    i
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    SSTa = SSTa+SST;
    tauxa = tauxa+taux;
    tauya = tauya+tauy;
end
SST = SSTa/length(outputs);
taux = tauxa/length(outputs);
tauy = tauya/length(outputs);

ylims = [-30 30];
xlims = [-260 -60];

[xL,yL] = size(lon);
[tmp ln1] = min(abs(lon(:,1)-xlims(1)));
[tmp ln2] = min(abs(lon(:,1)-xlims(2)));
[tmp lt1] = min(abs(lat(1,:)-ylims(1)));
[tmp lt2] = min(abs(lat(1,:)-ylims(2)));
xvec = (ln1-2):1:(ln2-2);
yvec = (lt1-2):1:(lt2-2);
[xL,yL] = size(lonu);
[tmp ln1] = min(abs(lonu(:,1)-xlims(1)));
[tmp ln2] = min(abs(lonu(:,1)-xlims(2)));
[tmp lt1] = min(abs(latu(1,:)-ylims(1)));
[tmp lt2] = min(abs(latu(1,:)-ylims(2)));
xvec2 = (ln1-2):15:(ln2-2);
yvec2 = (lt1-2):15:(lt2-2);
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = {[1],[3],[5],[7],[9],[11]};
% $$$           [1:12], ...
% $$$           [3], ...
% $$$           [7], ...
% $$$           [11]};

labels = {'Jan', ...
          'Mar', ...
          'May', ...
          'Jul', ...
          'Sep', ...
          'Nov'};
% $$$ 'Annual', ...
% $$$           'March', ...
% $$$           'July', ...
% $$$           'November'};

% $$$ %Colormap:
% $$$ clim = [0 0.1];
% $$$ sp = 0.01;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = parula(npts-3);
% $$$ cmap = parula(npts-3);
% $$$ % $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(end,:) = [1 1 1];
% $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ cmap = flipud(cmap);


clim = [20 30];
sp = 1;
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)
cmap = colormap(redblue(npts-3));
%colormap(flipud(lbmap(51,'RedBlue')));

figure;
%set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
set(gcf,'Position',[322          58        1247         945]);
% $$$ poss = [0.1300    0.4553    0.7693    0.4697; ...
% $$$         0.1300    0.1389    0.2343    0.2680; ...
% $$$         0.3951    0.1389    0.2343    0.2680; ...
% $$$         0.6681    0.1389    0.2343    0.2680];
poss = [0.0867    0.6910    0.4    0.26; ...
        0.5221    0.6910    0.4    0.26; ...
        0.0867    0.3926    0.4    0.26; ...
        0.5221    0.3926    0.4    0.26; ...
        0.0867    0.0899    0.4    0.26; ...
        0.5221    0.0899    0.4    0.26];
for i=1:length(months)
% $$$     if (i == 1)
% $$$         subplot(5,3,[1 9]);
% $$$     else
% $$$         subplot(5,3,[10 13]+(i-2));
% $$$     end
    subplot(3,2,i);
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    Z = monmean(SST(:,:,months{i}),3,ndays(months{i}));
    Z2 = monmean(taux(:,:,months{i}),3,ndays(months{i}));
    Z3 = monmean(tauy(:,:,months{i}),3,ndays(months{i}));
    Z = Z(xvec,yvec);
    contourf(X,Y,Z,cpts,'linestyle','none');
% $$$     contourf(lonu(xvec,yvec),latu(xvec,yvec),sqrt(Z2(xvec,yvec).^2+ ...
% $$$                                                   Z3(xvec,yvec).^2), ...
% $$$              [-1e10 0:0.01:0.5 1e10],'linestyle','none');
    hold on;
% $$$     [c,h] = contour(X,Y,Z,[-3:2:35],'-k');
% $$$     clabel(c,h);
    quiver(lonu(xvec2,yvec2),latu(xvec2,yvec2),Z2(xvec2,yvec2),Z3(xvec2,yvec2),0.75,'-k');
    caxis(clim);
    xlim(xlims);
    ylim(ylims);
    hold on;
    plot([-150 -90 -90 -150 -150],[-5 -5 5 5 -5],'-','color',[0 0.5 ...
                        0],'linewidth',2);
    if (i>=5)
        xlabel('Longitude ($^\circ$E)');
    else
        set(gca,'xticklabel',[]);
    end
    if (i==1 | i == 3 | i ==5 )
        ylabel('Latitude ($^\circ$N)');
    else
        set(gca,'yticklabel',[]);
    end
% $$$     if (i>1)
        text(-257,-25,labels{i},'BackgroundColor','w');
% $$$     else
% $$$         text(-278,35,labels{i},'BackgroundColor','w');
% $$$     end        
    if (i==2 | i == 4 | i ==6)
        cb = colorbar;
        ylabel(cb,'$^\circ$C');
    end
    set(gca,'Position',[poss(i,:)]);
    set(gca,'color','k');
    LabelAxes(gca,i,15,0.008,0.93);
end 
colormap(cmap);
% $$$ colormap(parula);



%%% Heat convergence into layer plots:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
model = 'MOM025';
outputs = [7];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
wname = sprintf('/srv/ccrc/data03/z3500785/MOM_HeatDiag/output%03d/ocean_wmass.nc',outputs(1));

% $$$ Tw = 22.5; %Warm temperature
% $$$ Tc = 15;   %Cool temperature
Tw = 30; %Warm temperature
Tc = 22.5;   %Cool temperature

[tmp iw] = min(abs(Te-Tw));
[tmp ic] = min(abs(Te-Tc));

FlM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlF = NaN*zeros(xL,yL,tL); % surface forcing
FlP = NaN*zeros(xL,yL,tL); % P-E+R
ty  = NaN*zeros(xL,yL,tL); % meridional mass flux
tx  = NaN*zeros(xL,yL,tL); % zonal mass flux
hy  = NaN*zeros(xL,yL,tL); % meridional heat flux
hx  = NaN*zeros(xL,yL,tL); % zonal heat flux
FlA = NaN*zeros(xL,yL,tL); % advection + submeso
FlT = NaN*zeros(xL,yL,tL); % tendency

for ti=1:tL
    ii = iw;
    sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,iw-ii+1,iw-ic+1)
    FlT(:,:,ti) = ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlA(:,:,ti) = ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlP(:,:,ti) = ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlM(:,:,ti) = ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlF(:,:,ti) = ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    ty(:,:,ti)  = ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
    tx(:,:,ti)  = ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
    hy(:,:,ti)  = Te(ii)*Cp*ty(:,:,ti);
    hx(:,:,ti)  = Te(ii)*Cp*tx(:,:,ti);

    for ii=iw-1:-1:ic
        sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,iw-ii+1,iw-ic+1)
    FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
                  ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    ty(:,:,ti)  = ty(:,:,ti) + ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
    tx(:,:,ti)  = tx(:,:,ti) + ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
    hy(:,:,ti)  = hy(:,:,ti) + Te(ii)*Cp*ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
    hx(:,:,ti)  = hx(:,:,ti) + Te(ii)*Cp*ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);

    end
end

save([base model sprintf('_output%03d',outputs(1)) '_VertInt_T' ...
      strrep(num2str(Tc),'.','p') 'C_T' strrep(num2str(Tw),'.','p') ...
      'C.mat'],'FlM','FlT','FlA','FlP','FlF','ty','tx','hy','hx','Tw','Tc');

load([base model sprintf('_output%03d',outputs(1)) '_VertInt_T' ...
      strrep(num2str(Tc),'.','p') 'C_T' strrep(num2str(Tw),'.','p') ...
      'C.mat']);
FlA = nansum(FlA,3)/12;FlA(FlA==0) = NaN;
FlT = nansum(FlT,3)/12;FlT(FlT==0) = NaN;
FlM = nansum(FlM,3)/12;FlM(FlM==0) = NaN;
FlF = nansum(FlF,3)/12;FlF(isnan(FlA)) = NaN;
FlP = nansum(FlP,3)/12;FlP(isnan(FlA)) = NaN;
ty = nansum(ty*1e9,3)/12;ty(ty==0) = NaN;
tx = nansum(tx*1e9,3)/12;tx(tx==0) = NaN;
hy = nansum(hy*1e9,3)/12;hy(hy==0) = NaN;
hx = nansum(hx*1e9,3)/12;hx(hx==0) = NaN;

xvec = 1:1:xL;yvec = 1:1:yL;
xvec2 = 1:5:xL-2;yvec2 = 1:5:yL-2;
xvec3 = 1:8:xL;yvec3 = 1:4:yL;

figure;
set(gcf,'Position',[3    40   956   963]);
subplot(3,1,1);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlF(xvec,yvec)+FlP(xvec,yvec));
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title(['Surface Forcing heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
caxis([-150 150]);
set(gca,'color','k');
ylim([-50 50]);
subplot(3,1,2);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlM(xvec,yvec));
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title(['Vertical Mixing heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
caxis([-150 150]);
set(gca,'color','k');
ylim([-50 50]);
subplot(3,1,3);
pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlA(xvec,yvec)-FlT(xvec,yvec));
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title(['TEN-ADV heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
caxis([-150 150]);
set(gca,'color','k');
ylim([-50 50]);
colormap(redblue);
subplot(4,1,4);
X = avg(avg(lon,1),2);
Y = avg(avg(lat,1),2);
CF = pi*6371000/180; %conversion factor from lat/lon to distance.
% $$$ hx(isnan(hx)) = 0;
% $$$ hy(isnan(hy)) = 0;
% $$$ DIV = avg(diff(hx,[],1)./diff(CF*lon.*cos(pi/180*lat),[],1),2) + ...
% $$$       avg(diff(hy,[],2)./diff(CF*lat,[],2),1);
DIV = avg(diff(hx,[],1)./avg(area,1),2) + ...
      avg(diff(hy,[],2)./avg(area,2),1);
pcolPlot(X(xvec2,yvec2),Y(xvec2,yvec2),DIV(xvec2,yvec2));
hold on;
hand = quiver(lon(xvec3,yvec3),lat(xvec3,yvec3),hx(xvec3,yvec3),hy(xvec3,yvec3),20,'-k');
% $$$ hand = quiver(lon(xvec3,yvec3),lat(xvec3,yvec3),hx(xvec3,yvec3)./sqrt(hx(xvec3,yvec3).^2+hy(xvec3,yvec3).^2),hy(xvec3,yvec3)./sqrt(hx(xvec3,yvec3).^2+hy(xvec3,yvec3).^2),'-k');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title(['Heat Fluxes and divergence between %3.1fC and %3.1fC',Tc,Tw));
caxis([-150 150]);
% $$$ ylim([-50 50]);
% $$$ ylim([-50 50]);
xlim([-260 -70]);
ylim([-40 40]);










