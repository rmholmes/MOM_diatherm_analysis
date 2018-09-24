% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

% $$$ % Load Base Variables:
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];

% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = 94;
% $$$ 
model = 'MOM025';
outputs = [15:19];

% $$$ model = 'MOM025_btide';
% $$$ outputs = [21];
% $$$ outputs = [12]
% $$$ 
% $$$ 
% $$$ model = 'MOM025_kb1em6';
% $$$ outputs = 30;

% $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ outputs = 36;

% $$$ model = 'MOM01';
% $$$ outputs = [222];

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
region = 'Global';
% $$$ region = 'Pacific';

%% Make Vars

direc = 1; % 0 = all warmer than Theta (native)
           % 1 = all colder than Theta
load([base model sprintf('_output%03d_',outputs(1)) 'ZAHBud.mat']);
names = fieldnames(ZA);

for i=1:length(outputs)
    load([base model sprintf('_output%03d_',outputs(i)) 'ZAHBud.mat']);
    
    for ii=1:length(names)
        eval(['z' names{ii} '(:,:,:,i) = ZA.' names{ii} ';']);
        if (direc)
            eval(['z' names{ii} '(:,:,:,i) = repmat(z' names{ii} '(:,1,:,i),[1 length(z' ...
                  names{ii} '(1,:,1)) 1 1]) - z' names{ii} '(:,:,:,i);']);
        end
    end
    if (direc)
        % Flip signs of fluxes:
        zF(:,:,:,i) = -zF(:,:,:,i);zM(:,:,:,i) = -zM(:,:,:,i);zP(:,:,:,i) = -zP(:,:,:,i);
        zSWH(:,:,:,i) = -zSWH(:,:,:,i);zJS(:,:,:,i) = -zJS(:,:,:,i);
        if (isfield(ZA,'Mkppiw')) % Vertical mixing components
            zMkppiw(:,:,:,i) = -zMkppiw(:,:,:,i);zMkppish(:,:,:,i) = -zMkppish(:,:,:,i);
            zMwave(:,:,:,i) = -zMwave(:,:,:,i);zMkppbl(:,:,:,i) = -zMkppbl(:,:,:,i);
            zMoth(:,:,:,i) = -zMoth(:,:,:,i);
        end
        if (isfield(ZA,'RED')) % Redi mixing
            zRED(:,:,:,i) = -zRED(:,:,:,i);zK33(:,:,:,i) = -zK33(:,:,:,i);
        end
    end
    % extras:
    if (direc)
        zAI(:,:,:,i) = cumsum(-rho0*Cp*zPSI(:,:,:,i)*dT,2); % Heat Function
    else
        zAI(:,:,:,i) = cumsum(rho0*Cp*zPSI(:,:,:,i)*dT,2,'reverse'); % Heat Function
    end
    
    tmp = -diff(zAI(:,:,:,i),[],1)./repmat(diff(latv',[],1),[1 TL+1 12]); %Diathermal component of heat function
    zJdia(:,:,:,i) = cat(1,cat(1,zeros(1,TL+1,tL),avg(tmp,1)),zeros(1,TL+1,tL));
    zJSH(:,:,:,i) = zJS(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;% Surface Volume flux base flux:
    zPI(:,:,:,i) = zP(:,:,:,i) - zJSH(:,:,:,i); % Interior heat source PI
    zN(:,:,:,i) = zdHdt(:,:,:,i) - zdVdt(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
    zI(:,:,:,i) = zJdia(:,:,:,i)-zN(:,:,:,i)-zF(:,:,:,i)-zPI(:,:,:,i)-zM(:,:,:,i);
    if (isfield(ZA,'RED'))
        zI(:,:,:,i) = zI(:,:,:,i) - zRED(:,:,:,i) - zK33(:,:,:,i);
    end
end
months = [1:length(zP(1,1,:,1))];

names = {names{:},'AI','Jdia','JSH','PI','N','I'};
for i=1:length(names)
    % Take mean across years:
    eval(['z' names{i} ' = mean(z' names{i} ',4);']);
% $$$ 
% $$$     % Apply latitude smoothing:
% $$$     eval(['z' names{i} ' = permute(filter_field(permute(mean(z' names{i} ...
% $$$           ',4),[3 2 1]),11,''-t''),[3 2 1]);']);    
end

if (direc)
    NaNs = monmean(zPSI,3,ndays) == repmat(monmean(zPSI(:,end,:),3,ndays),[1 TL+1 1]);
else
    NaNs = monmean(zPSI,3,ndays) == 0;
end

%% Plot lat-temp

% Load SST for plotting:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'SST');
SSTa = SST;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))],'SST');
    SSTa = SSTa+SST;
end
SST = SSTa/length(outputs);
SST(SST==0) = NaN;
meanSST = squeeze(nanmean(monmean(SST,3,ndays),1));
minSST = squeeze(min(monmean(SST,3,ndays),[],1));

% $$$ % Plot Streamfunction and Heat Function:
% $$$ fields = { ...
% $$$           {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-30 30],2,'Sv'}, ...
% $$$           {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ };

% Plot overall fields:
fields = { ...
          {-rho0*Cp*zPSI(:,:,months)/1e12, 'Isothermal Heat Flux $\mathcal{J}_{iso} = -\rho_0C_p\Psi$',[-150 150],15,'TW / $^\circ$C'}, ...
          {zJdia(:,:,months)/1e12, 'Diathermal Heat Flux $\mathcal{J}_{dia} = -\frac{\partial\mathcal{A}_I}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {(zF(:,:,months)+zPI(:,:,months))/1e12, 'Diathermal Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {zM(:,:,months)/1e12, 'Diathermal Vertical Mixing $\frac{\partial\mathcal{M}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zNUM(:,:,months)/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {zI(:,:,months)/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {(zRED(:,:,months)+zK33(:,:,months))/1e12, 'Diathermal Redi Mixing $\frac{\partial\mathcal{R}}{\partial\phi}$',[-10 10],1,'TW / $^\circ$'}, ...
          {zN(:,:,months)/1e12, 'Diathermal Internal Tendency $\frac{\partial\mathcal{N}}{\partial\phi}$',[-10 10],1,'TW / $^\circ$'}, ...
          };

% $$$ % Plot latitudinal integral:
% $$$ fields = { ...
% $$$           {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-30 30],2,'Sv'}, ...
% $$$           {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl((zF(:,:,months)+zPI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Surface Forcing $-(\mathcal{F}+\mathcal{P}_I)$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl(zM(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Vertical Mixing $-\mathcal{M}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl(zI(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Numerical Mixing $-\mathcal{I}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl(zN(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Tendency $-\mathcal{N}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           };

% $$$ % Plot mixing components:
% $$$ fields = { ...
% $$$           {zM(:,:,months)/1e12, 'Diathermal Vertical Mixing $\frac{\partial\mathcal{M}}{\partial\phi}$',[-40 0],2,'TW / $^\circ$'}, ...
% $$$           {zMkppiw(:,:,months)/1e12, 'Background Mixing',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zMkppish(:,:,months)/1e12, 'Interior Shear Instability',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zMkppbl(:,:,months)/1e12, 'KPP Boundary Layer',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zMwave(:,:,months)/1e12, 'Topographic Internal Wave',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           };

% $$$ % Plot perturbations from MOM025 control:
% $$$ load(['/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/MOM025_kb3seg_latT_save.mat']);
% $$$ fields = { ...
% $$$           {(zF(:,:,months)+zPI(:,:,months)-repmat(zFc,[1 1 length(months)]))/1e12, 'Diathermal Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-20 20],2,'TW / $^\circ$'}, ...
% $$$           {(zM(:,:,months)-repmat(zMc,[1 1 length(months)]))/1e12, 'Diathermal Vertical Mixing $\frac{\partial\mathcal{M}}{\partial\phi}$',[-20 20],2,'TW / $^\circ$'}, ...
% $$$           {(zI(:,:,months)-repmat(zIc,[1 1 length(months)]))/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-20 20],2,'TW / $^\circ$'}, ...
% $$$           };
% $$$ zAI = zAI - repmat(zAIc,[1 1 12]);

cpts = cell(1,length(fields));
for i=1:length(fields)
    cpts{i} = [-1e10 fields{i}{3}(1):fields{i}{4}:fields{i}{3}(2) 1e10];
end
npts = length(cpts{1});
clab = [1 0 0 0 0 0];

cmap = redblue(npts-3);
% $$$ cmap = parula(npts-3);
% $$$ cmap = parula(npts-3);
% $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

% $$$ AIsp = 0.25;
AIsp = 0.05;

% $$$ latfilt = 5;
latfilt = 11;

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
for i=1:length(fields)
    subplot(2,3,i)
    if (length(fields{i}{1}(1,:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    [Yg,Tg] = ndgrid(latv,x);

    VAR = monmean(fields{i}{1},3,ndays(months));
    VAR(VAR==0) = NaN;
    VAR(NaNs) = NaN;
    VAR = filter_field(VAR',latfilt,'-t')';
    contourf(Yg,Tg,VAR,cpts{i},'linestyle','none');
    hold on;
    col = [0 0 0];
% $$$     plot(latv,filter_field(meanSST,latfilt,'-t'),'--','color',col,'linewidth',2);
% $$$     plot(latv,filter_field(minSST,latfilt,'-t'),'--','color',col,'linewidth',2);
% $$$     plot(latv,21.5*ones(size(latv)),'--k');
    
    tmp = monmean(zAI,3,ndays)/1e15;
    tmp(NaNs) = NaN;
    tmp = filter_field(tmp',latfilt,'-t')';
% $$$     [c,h] = contour(Yg,Tg,tmp,[-3:AIsp:-AIsp],'--k');
% $$$     if (clab(i))
% $$$         clabel(c,h);
% $$$     end
% $$$     [c,h] =  contour(Yg,Tg,tmp,[AIsp:AIsp:3],'-k');
% $$$     if (clab(i))
% $$$         clabel(c,h);
% $$$     end
    tmp = filter_field(monmean(zPSI,3,ndays)'/1e6,latfilt,'-t')';
    contour(Yg,Tg,tmp,[-14.5 -14.5],'-','color',[0 0.5 0],'linewidth',2);
    contour(Yg,Tg,tmp,[1.5 1.5],'--','color',[0 0.5 0],'linewidth',2);
    
    ylim([-3 32]);
    xlim([-80 70]);
    caxis(fields{i}{3});
    box on; 
    grid on;
    ylabel('Temperature $\Theta$ ($^\circ$C)');
    xlabel('Latitude ($^\circ$N)');
    cb = colorbar;
    ylabel(cb,fields{i}{5});
    title(fields{i}{2});
    LabelAxes(gca,i,15,0.006,0.95);
end
colormap(cmap);

%% Calculate heat transports in different cells:

[Yg,Tg] = ndgrid(latv,Te);
psi = monmean(zPSI,3,ndays)/1e6;
psi(isnan(psi)) = 0;
% MOM025 Control:
regions = { ...
    {'Southern Subtropical Cell',3,psi>=3 & Tg>=9,'-',[-50 0],'r'}, ...
    {'Northern Subtropical Cell',-14.5,psi<=-14.5 & Tg>=12.75,'-',[0 40],'b'}, ...
    {'AABW Cell',3,psi>=3 & Tg<=9 & Yg <=-10,'-',[-72 -45],[0 0.5 0]}, ...
    {'NADW Cell',-14.5,psi<=-14.5 & Tg<=12.75 & Yg >=20,'-',[26 65],'m'}, ...
          };
% kb1em5:
regions = { ...
    {'Southern Subtropical Cell',3,psi>=3 & Tg>=8.9,'-',[-50 0],'r'}, ...
    {'Northern Subtropical Cell',-15.5,psi<=-15.5 & Tg>=12.75,'-',[0 40],'b'}, ...
    {'AABW Cell',3,psi>=3 & Tg<=8.9 & Yg <=-5,'-',[-72 -45],[0 0.5 0]}, ...
    {'NADW Cell',-15.5,psi<=-15.5 & Tg<=12.75 & Yg >=20,'-',[26 65],'m'}, ...
          };
% kb0:
regions = { ...
    {'Southern Subtropical Cell',1.5,psi>=1.5 & Tg>=9,'-',[-50 0],'r'}, ...
    {'Northern Subtropical Cell',-14.5,psi<=-14.5 & Tg>=12.75,'-',[0 40],'b'}, ...
    {'AABW Cell',1.5,psi>=1.5 & Tg<=9 & Yg <=-8,'-',[-72 -45],[0 0.5 0]}, ...
    {'NADW Cell',-14.5,psi<=-14.5 & Tg<=12.75 & Yg >=-25,'-',[26 65],'m'}, ...
          };

load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'mhflux');
mhfluxa = mhflux;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))],'mhflux');
    mhfluxa = mhfluxa+mhflux;
end
mhflux = monmean(mhfluxa/length(outputs)/1e15,2,ndays);

% Plot latitudinal integral:
fields = { ...
          {zAI(:,:,months)/1e15, 'Internal Heat Transport $\Delta\mathcal{A}_I$',[0 0.5 0],4}, ...
% $$$           {-cumsum(Repl((zF(:,:,months)+zPI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Surface Forcing $-\Delta(\mathcal{F}+\mathcal{P}_I)$','k',2}, ...
% $$$           {-cumsum(Repl(zM(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Vertical Mixing $-\Delta\mathcal{M}$','r',2}, ...
% $$$           {-cumsum(Repl(zI(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Numerical Mixing $-\Delta\mathcal{I}$','b',2}, ...
% $$$           {-cumsum(Repl(zN(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Tendency $-\Delta\mathcal{N}$','m',2}, ...
          };

% $$$ fields = { ...
% $$$           {-cumsum(Repl(zM(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Vertical Mixing $-\Delta\mathcal{M}$','r',2}, ...
% $$$           {-cumsum(Repl(zMkppiw(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Background Mixing','m',2}, ...
% $$$           {-cumsum(Repl(zMkppish(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Shear Instability','c',2}, ...
% $$$           {-cumsum(Repl(zMwave(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Topographic Internal Wave',[0 0.5 0],2}, ...
% $$$           {-cumsum(Repl(zMkppbl(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Boundary Layer Mixing','b',2}, ...
% $$$           };

% Get OHT due to each term:
OHT = cell(length(regions),length(fields));
for i=1:length(regions)
% $$$     subplot(2,2,i);
    maxT = NaN*zeros(size(latv));
    minT = NaN*zeros(size(latv));
    for k=1:length(latv)
        if (length(find(regions{i}{3}(k,:),1,'last')))
            maxTi = find(regions{i}{3}(k,:),1,'last');
            minTi = find(regions{i}{3}(k,:),1,'first');
            maxT(k) = interp1(psi(k,maxTi:maxTi+1),Te(maxTi:maxTi+1),regions{i}{2},'linear');
            minT(k) = interp1(psi(k,minTi-1:minTi),Te(minTi-1:minTi),regions{i}{2},'linear');
        end
    end

    for ii=1:length(fields)
        OHT{i,ii} = NaN*zeros(size(latv));
        extT = NaN*zeros(size(latv));
        VAR = monmean(fields{ii}{1},3,ndays);
        
        for k=1:length(latv)
            if (~isnan(maxT(k)))
                OHT{i,ii}(k) = interp1(Te,VAR(k,:),maxT(k),'linear')-interp1(Te,VAR(k,:),minT(k),'linear');
                extT(k) = (maxT(k)-minT(k))*(rho0*Cp*regions{i}{2}*1e6)/1e15;
                OHT{i,ii}(k) = OHT{i,ii}(k)+extT(k);
            end
        end
    end
end

% Total transport:
total = zeros(size(latv));
for ii=1:length(regions)
    total = total + Repl(OHT{ii,1},NaN,0);
end
mixed = mhflux - total';

% $$$ save('MOM025Control_Mhflux.mat','OHT','mhflux','mixed');

% Subtract off MOM025Control:
OHTa = OHT;
mhfluxa = mhflux;
mixeda = mixed;
load('MOM025Control_Mhflux.mat');
for ii=1:length(regions)
    for i=1:length(fields)
        OHT{ii,i} = OHTa{ii,i} - Repl(OHT{ii,i},NaN,0);
    end
end
mhflux = mhfluxa-Repl(mhflux,NaN,0);
mixed = mixeda-Repl(mixed,NaN,0);


latfilt = 25;


lg = {};
lgh = [];
figure;
set(gcf,'Position',[207          97        1609         815]);
set(gcf,'defaulttextfontsize',25);
set(gcf,'defaultaxesfontsize',25);

lgh(1) = plot(latv,filter_field(mhflux,latfilt,'-t'),'-k','linewidth',3);
hold on;
lg{1} = 'Total';
for i=1:length(regions)
    lgh(i+1) = plot(latv,Repl(filter_field(OHT{i,1},latfilt,'-t'),NaN,0),regions{i}{4},'color',regions{i}{6},'linewidth',3);
    lg{i+1} = regions{i}{1};
end
hold on;
lgh(i+2) = plot(latv,filter_field(mixed,latfilt,'-t'),'-','color',[0.6 ...
                    0.2 0],'linewidth',3);
lg{i+2} = 'Mixed';

xlabel('Latitude ($^\circ$N)');
ylabel('Meridional Heat Transport (PW)');
xlim([-80 80]);
grid on;
legend(lgh,lg);


    

