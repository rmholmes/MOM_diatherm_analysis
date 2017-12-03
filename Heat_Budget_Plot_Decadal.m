% This script makes plots of the heat budget in the decadal MOM
% simulations, annually-averaged.

close all;
clear all;

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_nipoall/mat_data/';
model = 'MOM025_nipoall';
outputs = [19];
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [7];

region = 'Pacific';

load([base model sprintf('_output%03d_',outputs(1)) region 'BaseVars.mat']);
ndays = diff(time_snap);

%% Global Calculations:
KEEPvars = {'PME','RMX','SWH','VDS','FRZ','ETS','VDF','KNL', ...
        'JS','P','M','F','A','N','I','Nmon','dVdt','dHdt','G', ...
        'JI','QI', ...
        'JBS','JSP','JITF','QBS','QSP','QITF', ...
        };

BUDGET = struct;

for i=1:length(outputs)
    
    load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
    
    % Calculate combined/extra quantities:
    P = PME+RMX; % PER effective heat flux (W)
    F = SWH+VDS+FRZ+ETS; % Surface heat flux (W)
    M = VDF+KNL; % Vertical mixing flux (W)
    AT= ADV-TEN; % Advection - Tendency to remove noise
    S = SUB; %Submesoscale term
    
    JS = SFW; % Surface volume flux
    
    Nmon = TENMON; %Monthly tendency term

    if (strcmp(region,'Pacific'))
        JI = JBS+JSP+JITF; %Combined volume flux out
        QI = QBS+QSP+QITF; %Combined heat flux out
    else
        QI = zeros(size(JS));
        JI = zeros(size(JS));
    end

    dVdt = diff(Vsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); %Interior volume change
    dHdt = diff(Hsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); %Interior heat content change

    G = dVdt - SFW + JI; % Interior across-isotherm volume flux
    
    Nhsnap = dHdt - dVdt*rho0*Cp.*repmat(Te,[1 tL]); % Hsnap total flux estimate
    N = Nhsnap;
    
    A = AT + N + S; % Total advection term by residual
    
    I = A + QI + (JS-JI).*repmat(Te,[1 tL])*rho0*Cp; % Implicit mixing term by residual

    for ii=1:length(KEEPvars)
        eval(['BUDGET(i).' KEEPvars{ii} ' = ' KEEPvars{ii} ';']);
    end
end

%% Subtract off background simulations:
baseBAK = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
modelBAK = 'MOM025';
%outputsBAK = [2 3 4 5 6];
outputsBAK = [7];

BAKBUDGET = struct;

for i=1:length(outputs)
    
    load([baseBAK modelBAK sprintf('_output%03d_',outputsBAK(i)) region 'HBud.mat']);
    
    % Calculate combined/extra quantities:
    P = PME+RMX; % PER effective heat flux (W)
    F = SWH+VDS+FRZ+ETS; % Surface heat flux (W)
    M = VDF+KNL; % Vertical mixing flux (W)
    AT= ADV-TEN; % Advection - Tendency to remove noise
    S = SUB; %Submesoscale term
    
    JS = SFW; % Surface volume flux
    
    Nmon = TENMON; %Monthly tendency term

    if (strcmp(region,'Pacific'))
        JI = JBS+JSP+JITF; %Combined volume flux out
        QI = QBS+QSP+QITF; %Combined heat flux out
    else
        QI = zeros(size(JS));
        JI = zeros(size(JS));
    end

    dVdt = diff(Vsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); %Interior volume change
    dHdt = diff(Hsnap,[],2)./repmat(diff(time_snap)'*86400,[TL+1 1]); %Interior heat content change

    G = dVdt - SFW + JI; % Interior across-isotherm volume flux
    
    Nhsnap = dHdt - dVdt*rho0*Cp.*repmat(Te,[1 tL]); % Hsnap total flux estimate
    N = Nhsnap;
    
    A = AT + N + S; % Total advection term by residual
    
    I = A + QI + (JS-JI).*repmat(Te,[1 tL])*rho0*Cp; % Implicit mixing term by residual

    for ii=1:length(KEEPvars)
        eval(['BAKBUDGET(i).' KEEPvars{ii} ' = ' KEEPvars{ii} ';']);
    end
end

% Subtract off background budget:
for ii=1:length(KEEPvars)
    eval(['BUDGET(i).' KEEPvars{ii} ' = BUDGET(i).' KEEPvars{ii} ' - ' ...
          'BAKBUDGET(i).' KEEPvars{ii} ';']);
end

% Take annual mean:
for ii=1:length(KEEPvars)
    eval(['BUDGET(i).' KEEPvars{ii} ' = monmean(BUDGET(i).' KEEPvars{ii} ',2,ndays);']);
end

%%%%Fluxes, Annual Average:
%tinds = [1:length(F(1,:))];
tinds = [1];
label = 'Annual Average';

% Interior sources/sinks:
fields = { ...
          {BUDGET.F(:,tinds)+BUDGET.P(:,tinds), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
          {BUDGET.M(:,tinds), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {BUDGET.I(:,tinds), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {BUDGET.N(:,tinds), 'Total Flux $\mathcal{N}$','m',2,'-'}, ...
          {BUDGET.dHdt(:,tinds), '$\partial \mathcal{H}/\partial t$',0.8*[1 1 1],2,'-'}, ...
          {-BUDGET.QI(:,tinds), 'Heat In',[0 0.5 0],2,'-'}, ...
          {-BUDGET.QSP(:,tinds), 'South Pacific',[0 0.5 0],2,'--'}, ...
          {-BUDGET.QITF(:,tinds), 'Indonesian Throughflow',[0 0.5 0],2,':'}, ...
          {-BUDGET.QBS(:,tinds), 'Bering Strait',[0 0.5 0],2,'-.'}, ...
          {BUDGET.G(:,tinds).*Te*rho0*Cp, 'Interior Transformation',0.5*[1 1 1],2,'--'}, ...
% $$$           {SW(:,tinds), 'Short-wave penetration','y',2,'--'}, ...
% $$$           {M(:,tinds)+I(:,tinds), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {P(:,tinds), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {VCT(:,tinds), '$\rho_0 C_p\Theta \partial \mathcal{V}/\partial t$',0.5*[1 1 1],2,'-'}, ...
% $$$           {S(:,tinds), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {Nmon(:,tinds), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$           {Nsnap(:,tinds), 'Calculated Total',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(:,tinds), 'Calculated Total 2',0.3*[1 1 1],2,'--'}, ...
          };

tint = 0; %Time-integrate and take last value (assumed yearly time slots)
dtr = 0; %Take derivative and look at transformations.
if (tint)
    if (dtr)
        Fscale = 1/1e12/rho0/Cp; % 1e12 m3
    else
        Fscale = 1/1e24; %yotta jules
    end
else
    if (dtr)
        Fscale = 1/1e6/rho0/Cp; %Sv
    else
        Fscale = 1/1e15; %PW
    end
end


figure;
set(gcf,'Position',[207          97        1609         815]);
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
    if (dtr)
        TT = T;
        ff = -diff(fields{i}{1},[],1)/dT;
    else
        TT = Te;
        ff = fields{i}{1};
    end        
    if (tint)
        % Time-integral:
        int = cumsum(ff*365*86400,2);
        legh(i) = plot(TT,int(:,end)*Fscale,fields{i}{5}, 'color',fields{i}{3} ...
                       , 'linewidth',fields{i}{4});
    else
        for j=1:length(tinds);
            legh(i) = plot(TT,ff(:,j)*Fscale,fields{i}{5}, 'color',fields{i}{3} ...
                           ,'linewidth',fields{i}{4});
        end
    end
    if (tint)
        leg{i} = ['$\int$ ' fields{i}{2} ' dt'];
    else
        leg{i} = fields{i}{2};
    end
end
if (tint)
    if (dtr)
    
    else
        ylim([-0.3 0.3]);
    end
else
    if (dtr)
    else
        ylim([-1.5 1.5]);
    end
end
xlim([-3 31]);
box on;
grid on;
if (tint)
    if (dtr)
        ylabel('Volume accumulated in fluid warmer than $\Theta$ (1e12 m$^3$)');
    else
        ylabel('Energy accumulated in fluid warmer than $\Theta$ (YJ)');
    end
else
    if (dtr)        
        ylabel('Transformation (Sv)');
    else
        ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
    end
end
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg,'FontSize',15);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%% Pacific Transports in and out:
figure;
set(gcf,'Position',[3    40   956   963]);
subplot(2,1,1);
plot(Te,TN/1e6);
hold on;
plot(Te,TS/1e6,'-b');
plot(Te,TW/1e6,'-r');
plot(Te,PER/1e6,'-','color',[0 0.5 0]);
plot(Te,(TW+TS+TN)/1e6,'--k');
xlabel('Temperature ($^\circ$C)');
ylabel('Volume Transport above $\Theta$ (Sv)');
legend('North Out','South Out','ITF','P-E+R','Total Out');

subplot(2,1,2);
plot(Te,HN/1e15);
hold on;
plot(Te,HS/1e15,'-b');
plot(Te,HW/1e15,'-r');
plot(Te,(HW+HS+HN)/1e15,'--k');
xlabel('Temperature ($^\circ$C)');
ylabel('Heat Transport above $\Theta$ (PW)');


%%% T-t plots: 
tinds = [1:length(F(1,:))];
fields = {
          {N(:,tinds,:), 'Total $\mathcal{N}$','m',2,'-'}, ...
          {F(:,tinds,:), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$           {M(:,tinds,:), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(:,tinds,:), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
          {M(:,tinds,:)+I(:,tinds,:), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {SW(:,tinds,:), 'Shortwave Penetration',[0 0.5 0],2,'--'}, ...
% $$$           {SW(:,tinds,:)+M(:,tinds,:)+I(:,tinds,:), 'Shortwave Penetration + $\mathcal{M}$ + $\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {P(:,tinds,:), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$           {S(:,tinds,:), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$           {F(:,tinds,:)+P(:,tinds,:), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$           {Nmon(:,tinds,:), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$           {Nsnap(:,tinds,:), 'Calculated Total',0.3*[1 1 1],2,':'}, ...
% $$$          {N2(:,tinds,:), 'Calculated Total 2',0.3*[1 1 1],2,'--'}, ...
          };

% Fluxes:
scale = 1/1e15;label = '(PW)';
% $$$ caxs = [-3 0];x = Te;
% $$$ sp = 0.1;
caxs = [-0.32 0.32];x = Te;
sp = 0.01;
cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

figure;
set(gcf,'Position',get(0,'ScreenSize'));
%set(gcf,'Position',[3    40   956   963]);
for ii=1:length(fields)
    subplot(1,length(fields),ii);
    V = fields{ii}{1}'*scale;
    [X,Y] = ndgrid(tinds,Te);
    contourf(X,Y,V,cint);%,'linestyle','none');
    cb = colorbar('Location','NorthOutside','FontSize',25);    
    set(gca,'ytick',-5:5:35);
    set(gca,'xtick',[0:5:tinds(end)]);
    set(gca,'xticklabel',[0:5:tinds(end)]);
    xlim([0 tinds(end)]);
    ylim([-3 31]);
    grid on;
    caxis(caxs);
    xlabel('Year');
    ylabel('Temperature ($^\circ$C)');
    xlabel(cb,['MOM025 ' fields{ii}{2} ' ' ...
               label],'FontSize',20);
    set(gca,'FontSize',25);
end
colormap(redblue);
% $$$ cmap = redblue((length(cint)-3)*2);
% $$$ cmap = cmap(1:(length(cint)-3),:);
% $$$ colormap(cmap);
% $$$ caxis([-0.8 0]);
% $$$ hold on;
% $$$ plot([1 tL],[22.5 22.5],'--k','linewidth',2);
% $$$ LabelAxes(gca,3,25,0.008,0.965);


% H and V distributions:
figure;
subplot(1,2,1);
for i=tinds
    hold on;
% $$$     plot(VsnapAA(:,i),Te,'-k');
    plot(-diff(VsnapAA(:,i),[],1)/dT,T,'-k');
end
box on;
grid on;
% $$$ xlabel(' m$^3$');
xlabel(' m$^3$ / $^\circ$C');
ylabel('Temperature ($^\circ$C)');
%title('$V(\Theta,t)-V_{steady}(\Theta)$');
title('$V(\Theta,t)$');
subplot(1,2,2);
for i=tinds
    hold on;
% $$$     plot(HsnapAA(:,i),Te,'-k');
    plot(-diff(HsnapAA(:,i),[],1)/dT,T,'-k');
end
% $$$ xlabel('J');
xlabel('J / $^\circ$C');
ylabel('Temperature ($^\circ$C)');
box on;
grid on;
title('$H(\Theta,t)$');

%% Plot difference first to last year of shflux/SST/taux/PSI(lat,T):

% This script makes plots of the heat budget in the decadal MOM
% simulations, annually-averaged.

close all;
clear all;

% Load BAK:
model = 'MOM025';
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
outputs = [2 3 4 5 6];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
shfluxa = shflux;
SSTa = SST;
tauxa = taux;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    shfluxa = shfluxa+shflux;
    SSTa = SSTa+SST;
    tauxa = tauxa+taux;
end
shflux = shfluxa/length(outputs);
SST = SSTa/length(outputs);
taux = tauxa/length(outputs);

shfluxBAK = mean(shflux,3);
SSTBAK = mean(SST,3);
tauxBAK = mean(taux,3);
clear shfluxa SSTa tauxa shflux SST taux;

% Load BAK PSI:
model = 'MOM025';
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
outputs = [7];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
load([base model sprintf('_output%03d_Tpsi.mat',outputs(1))]);
NaNs = PSI == 0;
PSI = cumsum(-PSI/rho0*1e9,2);
PSI(NaNs) = NaN;
PSIBAK = mean(PSI,3);


% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_nipoall/mat_data/';
model = 'MOM025_nipoall';
outputs = [19];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
shflux = mean(shflux,3)-shfluxBAK;
SST = mean(SST,3)-SSTBAK;
taux = mean(taux,3)-tauxBAK;

% Load PSI:
model = 'MOM025_nipoall';
base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_nipoall/mat_data/';
outputs = [19];
load([base model sprintf('_output%03d_Tpsi.mat',outputs(1))]);
NaNs = PSI == 0;
PSI = cumsum(-PSI/rho0*1e9,2);
PSI(NaNs) = NaN;
PSI = mean(PSI,3)-PSIBAK;

xvec = 1:4:xL;
yvec = 1:4:yL;

figure;
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);
contourf(lon(xvec,yvec),lat(xvec,yvec),shflux(xvec,yvec),[-1e10 -100:5:100 1e10],'linestyle','none');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title('Year 20 - Year 0 net surface heat flux (Wm$^{-2}$)');
set(gca,'color','k');
caxis([-80 80]);
cb = colorbar;

subplot(2,2,2);
contourf(lon(xvec,yvec),lat(xvec,yvec),taux(xvec,yvec),[-1e10 -0.05:0.001:0.05 1e10],'linestyle','none');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title('Year 20 - Year 0 zonal wind stress (Nm$^{-2}$)');
set(gca,'color','k');
caxis([-0.05 0.05]);
cb = colorbar;

subplot(2,2,3);
contourf(lon(xvec,yvec),lat(xvec,yvec),SST(xvec,yvec),[-1e10 -3:0.05:3 1e10],'linestyle','none');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title('Year 20 - Year 0 SST ($^\circ$C)');
set(gca,'color','k');
caxis([-3 3]);
cb = colorbar;

subplot(2,2,4);
[YY,TT] = ndgrid(latv,T);
contourf(YY,TT,PSI/1e6,[-1e10 -10:1:10 1e10]);
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
title('Year 20 - Year 0 $\Psi(y,\Theta)$ (Sv)');
caxis([-6 6]);
cb = colorbar;
colormap(redblue);



