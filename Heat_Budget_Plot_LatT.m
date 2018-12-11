% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

% $$$ % $$$ % Load Base Variables:
model = 'MOM025_kb3seg';
outputs = [86:90];
output = 86;

% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = [95:99];
% $$$ 
% $$$ model = 'MOM025';
% $$$ outputs = [15:19];

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

%% Make Vars
load([base model sprintf('_output%03d_',outputs(1)) 'ZAHBud.mat']);
names = fieldnames(ZA);

for i=1:length(outputs)
    load([base model sprintf('_output%03d_',outputs(i)) 'ZAHBud.mat']);
    
    % Flip sign of fluxes and save:
    for ii=1:length(names)
        eval(['z' names{ii} '(:,:,:,i) = -ZA.' names{ii} ';']);
    end
    
    % Non-fluxes (don't flip sign):
    flbck = {'dVdt','dHdt','PSI','AHD'};
    for iii = 1:length(flbck)
        eval(['z' flbck{iii} '(:,:,:,i) = ZA.' flbck{iii} ';']);
    end

    % Heat Function:
    zAIpsi(:,:,:,i) = -rho0*Cp*cumsum(zPSI(:,:,:,i)*dT,2); % Heat Function (defined on Tc, since Psi is on Te, and defined on v-points)
    zAI(:,:,:,i) = zAHD(:,:,:,i) - rho0*Cp*repmat(Te',[yL 1 tL]).*zPSI(:,:,:,i); % (defined on Te and v-points)  USE THIS!
    
    zJdia(:,:,:,i) = -diff(cat(1,zeros(1,TL+1,tL),zAI(:,:,:,i)),[],1)./repmat(yuo,[1 TL+1 tL]);
    zJSH(:,:,:,i) = zJS(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
    zPI(:,:,:,i) = zP(:,:,:,i) - zJSH(:,:,:,i);
    zN(:,:,:,i) = zdHdt(:,:,:,i) - zdVdt(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
    zI(:,:,:,i) = zJdia(:,:,:,i)-zN(:,:,:,i)-zF(:,:,:,i)-zPI(:,:,:,i)-zM(:,:,:,i);
    if (isfield(ZA,'RED'))
        zI(:,:,:,i) = zI(:,:,:,i) - zRED(:,:,:,i) - zK33(:,:,:,i);
    end
    
    zZT(:,:,i) = ZAtemp;
    zT(:,:,i) = tempZA;
    zZTx(:,:,i) = ZMtemp;
    zTx(:,:,i) = tempZM;
    zZTxa(:,:,i) = ZMAtemp;
    zTxa(:,:,i) = tempZMA;
    
    MHT(:,:,i) = zAHD(:,end,:,i);
end
months = [1:length(zP(1,1,:,1))];

% Take mean across years:
zZT = mean(zZT,3);
zT = mean(zT,3);
zZTx = mean(zZTx,3);
zTx = mean(zTx,3);
zZTxa = mean(zZTxa,3);
zTxa = mean(zTxa,3);
MHT = mean(MHT,3);
names = {names{:},'AI','AIpsi','Jdia','JSH','PI','N','I'};
for i=1:length(names)
    eval(['z' names{i} ' = mean(z' names{i} ',4);']);
end

% Generate a NaNs array:
tmp = monmean(zPSI,3,ndays);
NaNs = zeros(size(tmp));
for i = 1:(TL+1)
    NaNs(:,i) = tmp(:,i) == tmp(:,end);
end
% Max SST (not monthly, ever) line and add to NaNs (note: on Te):
maxT = zeros(yL,1);
maxTi = zeros(yL,1);
psiT = zeros(yL,1);
for i = 1:yL 
    ind = find(NaNs(i,:),1,'first');
    maxTi(i) = ind;
    maxT(i) = Te(ind);
    psiT(i) = tmp(i,ind);
    NaNs(i,ind) = 0;
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
% $$$ maxSST = squeeze(max(monmean(SST,3,ndays),[],1));

% $$$ % Plot Streamfunction and Heat Function:
% $$$ fields = { ...
% $$$           {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-30 30],2,'Sv'}, ...
% $$$           {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {zAIpsi(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$ from $\Psi$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {zAHD(:,:,months)/1e15, 'Heat Function $\mathcal{A}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {(zAI(:,:,months)-zAIpsi(:,:,months))/1e15, 'Heat Function diff $\mathcal{A}$',[-0.01 0.01],0.005,'PW'}, ...
% $$$ };

% Plot overall fields:
% $$$ fields = { ...
% $$$           {-rho0*Cp*zPSI(:,:,months)/1e12, 'Isothermal Heat Flux $\mathcal{J}_{iso} = -\rho_0C_p\Psi$',[-150 150],15,'TW / $^\circ$C'}, ...
% $$$           {zJdia(:,:,months)/1e12, 'Diathermal Heat Flux $\mathcal{J}_{dia} = -\frac{\partial\mathcal{A}_I}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {(zF(:,:,months)+zPI(:,:,months))/1e12, 'Diathermal Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zM(:,:,months)/1e12, 'Diathermal Vertical Mixing $\frac{\partial\mathcal{M}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zI(:,:,months)/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zN(:,:,months)/1e12, 'Diathermal Internal Tendency $\frac{\partial}{\partial\phi}\frac{\partial\mathcal{H}_I}{\partial t}$',[-10 10],1,'TW / $^\circ$'}, ...
% $$$           };
fields = { ...
% $$$           {-rho0*Cp*zPSI(:,:,months)/1e12, 'Isothermal Heat Flux $\mathcal{J}_{iso} = -\rho_0C_p\Psi$',[-150 150],15,'TW / $^\circ$C'}, ...
          {zJdia(:,:,months)/1e12, 'Diathermal Heat Flux $\mathcal{J}_{dia} = -\frac{\partial\mathcal{A}_I}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {zN(:,:,months)/1e12, 'Tendency $\frac{\partial}{\partial\phi}\frac{\partial\mathcal{H}_I}{\partial t}$',[-25 25],1.25,'TW / $^\circ$'}, ...
          {(zF(:,:,months)+zPI(:,:,months))/1e12, 'Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {zM(:,:,months)/1e12+zI(:,:,months)/1e12, 'Mixing $\frac{\partial(\mathcal{M}+\mathcal{I})}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zNUM(:,:,months)/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zI(:,:,months)/1e12, 'Diathermal Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {(zRED(:,:,months)+zK33(:,:,months))/1e12, 'Diathermal Redi Mixing $\frac{\partial\mathcal{R}}{\partial\phi}$',[-10 10],1,'TW / $^\circ$'}, ...
          };

% $$$ % Plot latitudinal integral:
% $$$ fields = { ...
% $$$ % $$$           {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-30 30],2,'Sv'}, ...
% $$$           {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl((zF(:,:,months)+zPI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Surface Forcing $-(\mathcal{F}+\mathcal{P}_I)$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {-cumsum(Repl((zM(:,:,months)+zI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Total Mixing $-(\mathcal{M}+\mathcal{I})$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {-cumsum(Repl(zI(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Numerical Mixing $-\mathcal{I}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {-cumsum(Repl(zM(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Vertical Mixing $-\mathcal{M}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {-cumsum(Repl(zI(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Numerical Mixing $-\mathcal{I}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$ % $$$           {-cumsum(Repl(zN(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Diathermal Tendency $-\frac{\partial\mathcal{H}_I}{\partial t}$',[-0.25 0.25],0.025,'PW'}, ...
% $$$           };
% $$$ 
% $$$ % Plot latitidinal integral from zero AI line:
% $$$ tmp = monmean(zAI(:,:,months),3,ndays);
% $$$ tmp(latv>60,:,:) = 1e15;
% $$$ zinds = zeros(TL+1,1);
% $$$ lats = zeros(TL+1,1);
% $$$ for i=1:TL+1
% $$$     ind = find(tmp(:,i)<0,1,'last');
% $$$     if (length(ind)>0)
% $$$         zinds(i) = ind;
% $$$     else
% $$$         zinds(i) = 1;
% $$$     end
% $$$     lats(i) = latv(zinds(i));
% $$$ end
% $$$ for i=2:length(fields)
% $$$     for ii=1:TL+1
% $$$     fields{i}{1}(:,ii,:) = fields{i}{1}(:,ii,:) - repmat(fields{i}{1}(zinds(ii),ii,:),[yL ...
% $$$                         1 1]);
% $$$     end
% $$$ end
% $$$ 
% $$$ % Pull out totals:
% $$$ MHTtot = cell(length(fields)-1,2);
% $$$ for i=1:(length(fields))
% $$$     tmp = monmean(fields{i}{1},3,ndays);
% $$$     MHTtot{i,1} = zeros(yL,1);
% $$$     for ii=1:yL
% $$$         MHTtot{i,1}(ii) = tmp(ii,maxTi(ii));
% $$$     end
% $$$     MHTtot{i,2} = fields{i}{2};
% $$$ end     

% $$$ % Plot mixing components:
% $$$ fields = { ...
% $$$           {zM(:,:,months)/1e12, 'Total Vertical Mixing $\frac{\partial\mathcal{M}}{\partial\phi}$',[-40 0],2,'TW / $^\circ$'}, ...
% $$$           {zMkppiw(:,:,months)/1e12, 'Background Mixing',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zMkppish(:,:,months)/1e12, 'Interior Shear Instability',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zI(:,:,months)/1e12, 'Numerical Mixing $\frac{\partial\mathcal{I}}{\partial\phi}$',[-40 0],2,'TW / $^\circ$'}, ...
% $$$           {zMkppbl(:,:,months)/1e12, 'KPP Boundary Layer',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           {zMwave(:,:,months)/1e12, 'Internal Tide',[-20 0],0.5,'TW / $^\circ$'}, ...
% $$$           };

% $$$ % Plot perturbations from MOM025 control:
% $$$ zFc = monmean(zF+zPI,3,ndays(months));
% $$$ zMc = monmean(zM,3,ndays(months));
% $$$ zIc = monmean(zI,3,ndays(months));
% $$$ zAIc = monmean(zAI,3,ndays(months));
% $$$ save('/srv/ccrc/data03/z3500785/mom/mat_data/MOM025_kb3seg_86to90_latT_save.mat','zFc','zMc','zIc','zAIc');
% $$$ load('/srv/ccrc/data03/z3500785/mom/mat_data/MOM025_kb3seg_86to90_latT_save.mat');
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
% $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

AIsp = 0.25;
% $$$ AIsp = 0.025;

% $$$ latfilt = 5;
latfilt = 1;

doZAremap = 0; % remap to depth space

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
set(gcf,'defaulttextfontsize',15);
 set(gcf,'defaultaxesfontsize',15);

% $$$ set(gcf,'Position',[57         115        1603         654]);
% $$$ set(gcf,'defaulttextfontsize',20);
% $$$ set(gcf,'defaultaxesfontsize',20);

% $$$ set(gcf,'Position',[1728          93        1101         815]);

for i=1:length(fields)
    if (i >=3)
        subplot(2,3,i+1);
    else
        subplot(2,3,i);
% $$$         subplot(3,1,i);
    end
    if (length(fields{i}{1}(1,:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    if (doZAremap)
        Yg = repmat(latv',[1 TL+1]);
        Tg = zZTxa;
    else
        [Yg,Tg] = ndgrid(latv,x);
    end

    VAR = monmean(fields{i}{1},3,ndays(months));
    VAR(VAR==0) = NaN;
    VAR(NaNs==1) = NaN;
    VAR = filter_field(VAR',latfilt,'-t')';
    contourf(Yg,Tg,VAR,cpts{i},'linestyle','none');
    hold on;
    col = [0 0 0];

    if (~doZAremap)
        plot(latv,filter_field(meanSST,latfilt,'-t'),'--','color',col,'linewidth',2);
        plot(latv,filter_field(minSST,latfilt,'-t'),'--','color',col,'linewidth',2);
    plot(latv,filter_field(maxT,latfilt,'-t'),':k');
% $$$     plot(latv,21.5*ones(size(latv)),'--k');
% $$$         plot(lats,Te,'-k','linewidth',2);%'color',[0 0.5 0]);

% $$$     tmp = filter_field(monmean(zPSI,3,ndays)'/1e6,latfilt,'-t')';
% $$$     tmp(Yg<20 & (Tg > 4.9 & Tg < 8.75)) = NaN;
% $$$     contour(Yg,Tg,tmp,[-14.5 -14.5],'-','color',[0 0.5 0],'linewidth',2);
% $$$     contour(Yg,Tg,tmp,[5 5],'--','color',[0 0.5 0],'linewidth',2);
    else
        [tt,zz] = ndgrid(latv,-z);
        [c,h] = contour(tt,zz,zTxa,[0:4:34],'-k');
        clabel(c,h);
    end

% $$$     if (i==1)
    tmp = monmean(zAI,3,ndays)/1e15;
    tmp(NaNs==1) = NaN;
    tmp = filter_field(tmp',latfilt,'-t')';
    [c,h] = contour(Yg,Tg,tmp,[-3:AIsp:-AIsp],'--k');
    if (clab(i))
        clabel(c,h);
    end
    [c,h] =  contour(Yg,Tg,tmp,[AIsp:AIsp:3],'-k');
    if (clab(i))
        clabel(c,h);
    end
% $$$     end

    if (doZAremap)
        ylabel('Depth (m)');
        ylim([-1050 0]);
    else
        ylabel('Temperature $\Theta$ ($^\circ$C)');
        ylim([-3 34]);
    end
    xlim([-70 70]);
    caxis(fields{i}{3});
    box on; 
    grid on;
    xlabel('Latitude ($^\circ$N)');
    cb = colorbar;
    ylabel(cb,fields{i}{5});
    title(fields{i}{2});
    LabelAxes(gca,i,15,0.006,0.95);
end
colormap(cmap);

%% Plot total MHT contributions:
% Plot latitudinal integral:
fields = { ...
          {zAI(:,:,months)/1e15, 'Internal Heat Transport $\mathcal{A}_I^{max}$',[-1.5 1.5],0.1,'PW'}, ...
          {-cumsum(Repl((zF(:,:,months)+zPI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Surface Forcing $-(\mathcal{F}+\mathcal{P}_I)^{max}$',[-1.5 1.5],0.1,'PW'}, ...
          {-cumsum(Repl((zM(:,:,months)+zI(:,:,months)).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Mixing $-(\mathcal{M}+\mathcal{I})^{max}$',[-1.5 1.5],0.1,'PW'}, ...
          {-cumsum(Repl(zN(:,:,months).*repmat(yto,[1 TL+1 tL]),NaN,0),1)/1e15,'Tendency $-\frac{\partial\mathcal{H}_I}{\partial t}^{max}$',[-0.25 0.25],0.025,'PW'}, ...
          };

% Plot latitidinal integral from zero AI line:
tmp = monmean(zAI(:,:,months),3,ndays);
tmp(latv>60,:,:) = 1e15;
zinds = zeros(TL+1,1);
lats = zeros(TL+1,1);
for i=1:TL+1
    ind = find(tmp(:,i)<0,1,'last');
    if (length(ind)>0)
        zinds(i) = ind;
    else
        zinds(i) = 1;
    end
    lats(i) = latv(zinds(i));
end
for i=2:length(fields)
    for ii=1:TL+1
    fields{i}{1}(:,ii,:) = fields{i}{1}(:,ii,:) - repmat(fields{i}{1}(zinds(ii),ii,:),[yL ...
                        1 1]);
    end
end

% Pull out totals:
MHTtot = cell(length(fields)-1,2);
for i=1:(length(fields))
    tmp = monmean(fields{i}{1},3,ndays);
    MHTtot{i,1} = zeros(yL,1);
    for ii=1:yL
        MHTtot{i,1}(ii) = tmp(ii,maxTi(ii));
    end
    MHTtot{i,2} = fields{i}{2};
end     

figure;
colors = {'-k','-b','-r','-m','-c'};
lfilt = 11;
plot(latv,filter_field(MHT,lfilt,'-t'),'--k','linewidth',3);
hold on; 
plot(latv,filter_field(rho0*Cp*psiT.*maxT,lfilt,'-t')/1e15,':k','linewidth',3);
for i=1:4
    plot(latv,filter_field(MHTtot{i,1},lfilt,'-t'),colors{i},'linewidth',3);
end
% $$$ plot(latv,filter_field(MHTtot{2,1}+MHTtot{3,1}+MHTtot{4,1},1,'-t'),'-c');
% $$$ plot(latv,filter_field(MHTtot{1,1},1,'-t')+rho0*Cp*psiT.*maxT/1e15,'-y');
legend({'Total Heat Transport $\mathcal{A}^{max}$','External Heat Transport $\rho_0C_p\Psi^{max}\Theta^{max}$',MHTtot{:,2},});
xlabel('Latitude ($^\circ$N)');
ylabel('PW');
xlim([-90 90]);
grid on;




%% Calculate heat transports in different cells:

[Yg,Tg] = ndgrid(latv,Te);
psi = filter_field(monmean(zPSI,3,ndays)'/1e6,latfilt,'-t')';
psi(isnan(psi)) = 0;
% MOM025 Control:
% $$$ regions = { ... %OLD years:
% $$$     {'Southern Subtropical Cell',3,psi>=3 & Tg>=9,'-',[-50 0],'r'}, ...
% $$$     {'Northern Subtropical Cell',-14.5,psi<=-14.5 & Tg>=12.75,'-',[0 40],'b'}, ...
% $$$     {'AABW Cell',3,psi>=3 & Tg<=9 & Yg <=-10,'-',[-72 -45],[0 0.5 0]}, ...
% $$$     {'NADW Cell',-14.5,psi<=-14.5 & Tg<=12.75 & Yg >=20,'-',[26 65],'m'}, ...
% $$$           };
% $$$ regions = { ... % NEW years:
% $$$     {'SSTC',5,psi>=5 & Tg>=8.75,'-',[-50 0],'r'}, ...
% $$$     {'NSTC',-14.5,psi<=-14.5 & Tg>=12,'-',[0 40],'b'}, ...
% $$$     {'AABW',5,psi>=5 & Tg<8.75 & Yg <=-10,'-',[-72 -45],[0 0.5 0]}, ...
% $$$     {'NADW',-14.5,psi<=-14.5 & Tg<12 & Yg >=20,'-',[26 65],'m'}, ...
% $$$           };
% $$$ % kb1em5:
% $$$ regions = { ... %OLD years:
% $$$     {'Southern Subtropical Cell',3,psi>=3 & Tg>=8.9,'-',[-50 0],'r'}, ...
% $$$     {'Northern Subtropical Cell',-15.5,psi<=-15.5 & Tg>=12.75,'-',[0 40],'b'}, ...
% $$$     {'AABW Cell',3,psi>=3 & Tg<=8.9 & Yg <=-5,'-',[-72 -45],[0 0.5 0]}, ...
% $$$     {'NADW Cell',-15.5,psi<=-15.5 & Tg<=12.75 & Yg >=20,'-',[26 65],'m'}, ...
% $$$           };
% $$$ regions = { ... % NEW years:
% $$$     {'SSTC',5,psi>=5 & Tg>=8.75,'-',[-50 0],'r'}, ...
% $$$     {'NSTC',-14.5,psi<=-14.5 & Tg>=12,'-',[0 40],'b'}, ...
% $$$     {'AABW',5,psi>=5 & Tg<8.75 & Yg <=-10,'-',[-72 -45],[0 0.5 0]}, ...
% $$$     {'NADW',-14.5,psi<=-14.5 & Tg<12 & Yg >=20,'-',[26 65],'m'}, ...
% $$$           };
% $$$ % kb0:
% $$$ regions = { ... %OLD years:
% $$$     {'Southern Subtropical Cell',1.5,psi>=1.5 & Tg>=9,'-',[-50 0],'r'}, ...
% $$$     {'Northern Subtropical Cell',-14.5,psi<=-14.5 & Tg>=12.75,'-',[0 40],'b'}, ...
% $$$     {'AABW Cell',1.5,psi>=1.5 & Tg<=9 & Yg <=-8,'-',[-72 -45],[0 0.5 0]}, ...
% $$$     {'NADW Cell',-14.5,psi<=-14.5 & Tg<=12.75 & Yg >=-25,'-',[26 65],'m'}, ...
% $$$           };
regions = { ... % NEW years:
    {'SSTC',5,psi>=5 & Tg>=8.75,'-',[-50 0],'r'}, ...
    {'NSTC',-14.5,psi<=-14.5 & Tg>=12,'-',[0 40],'b'}, ...
    {'AABW',5,psi>=5 & Tg<8.75 & Yg <=-10,'-',[-72 -45],[0 0.5 0]}, ...
    {'NADW',-14.5,psi<=-14.5 & Tg<12 & Yg >=3.8,'-',[26 65],'m'}, ...
          };

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

mhflux = monmean(MHT,2,ndays);

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


latfilt = 5;


lg = {};
lgh = [];
figure;
set(gcf,'Position',[207          97        1609         815]);
set(gcf,'defaulttextfontsize',30);
set(gcf,'defaultaxesfontsize',30);

cols = {};
lgh(1) = plot(latv,filter_field(mhflux,latfilt,'-t'),'-k','linewidth',3);
hold on;
lg{1} = 'Total';
cols{1} = 'k';
for i=1:length(regions)
    lgh(i+1) = plot(latv,Repl(filter_field(OHT{i,1},latfilt,'-t'),NaN,0),regions{i}{4},'color',regions{i}{6},'linewidth',3);
    lg{i+1} = regions{i}{1};
    cols{i+1} = regions{i}{6};
end
hold on;
lgh(i+2) = plot(latv,filter_field(mixed,latfilt,'-t'),'-','color',[0.6 ...
                    0.2 0],'linewidth',3);
lg{i+2} = 'Mixed';
cols{i+2} = [0.6 0.2 0];

xlabel('Latitude ($^\circ$N)');
ylabel('Meridional Heat Transport (PW)');
xlim([-80 80]);
grid on;
for i=1:length(lg)
    text(0,0,lg{i},'color',cols{i});
end

legend(lgh,lg);


%%%% Heat Uptake Experiment:

% Plot SO zonal averages:
[Y,Z] = ndgrid(latv,-z);
figure;
[c,h] = contourf(Y,Z,tempZA,[-2:1:34],'linestyle','none');
clabel(c,h);
hold on;
[c,h] = contour(Y,Z,rho0,[1020:0.2:1028],'-k','linewidth',2);
clabel(c,h);
xlim([-60 -10]);
ylim([-1200 0]);
xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
caxis([0 24]);
colorbar;
colormap(jet);

% Calculate volume in T-layers:
load([base model sprintf('_output%03d_',outputs(1)) 'VHZAf.mat']);
V = cat(2,zeros(yL,1,12),cumsum(V,2));
HD = cat(2,zeros(yL,1,12),cumsum(rho0*Cp*repmat(T',[yL 1 tL]).*diff(V,[],2)/dT*dT,2));
H = cat(2,zeros(yL,1,12),cumsum(H,2));
HE = rho0*Cp*repmat(Te',[yL 1 tL]).*V;
HI = H-He;

reg = latv>-40 & latv<-15;
Vreg = nansum(nanmean(V(reg,:,:),3),1)';

Hreg = nansum(nanmean(H(reg,:,:),3),1)';
HDreg = nansum(nanmean(HD(reg,:,:),3),1)';
HEreg = nansum(nanmean(HE(reg,:,:),3),1)';
HIreg = nansum(nanmean(HI(reg,:,:),3),1)';

wc = 14;
wid = 2;
Vper = 0.1e17*(Te-(wc-wid))/(wid.^2).*exp(-(Te-(wc-wid)).^2/(wid.^2));

VPreg = Vreg+cumsum(Vper,1);
HPreg = cat(1,[0],cumsum(rho0*Cp*T.*diff(VPreg,[],1)/dT*dT,1));
HEPreg = rho0*Cp*Te.*VPreg;
HIPreg = HPreg-HEPreg;
% $$$ HIPreg = -rho0*Cp*cumsum(VPreg*dT,1);

% $$$ figure;
clf;
subplot(1,2,1);
plot(Vreg/1e15,Te,'-k','linewidth',2);
hold on;
plot(VPreg/1e15,Te,'--k','linewidth',2);
ylim([0 20]);
xlim([0 300]);
ylabel('Temperature $\Theta$ ($^\circ$C)');
xlabel('Volume Below $\Theta$ ($10^{15}$m$^3$)');
title('Volume $40^\circ$S-$15^\circ$S');
legend('MOM025','Perturbed');
% $$$ subplot(1,3,2);
% $$$ plot(Hreg,Te,'-k');
% $$$ hold on;
% $$$ plot(HEreg,Te,'-m');
% $$$ plot(HIreg,Te,'-c');
% $$$ ylim([0 25]);
% $$$ % $$$ plot(HDreg,Te,'-b');
% $$$ plot(HPreg,Te,'--k');
% $$$ hold on;
% $$$ plot(HEPreg,Te,'--m');
% $$$ plot(HIPreg,Te,'--c');
% $$$ ylim([0 25]);
subplot(1,2,2);
plot((HPreg-Hreg)/1e21,Te,'-k','linewidth',2);
hold on;
plot((HEPreg-HEreg)/1e21,Te,'-b','linewidth',2);
plot((HIPreg-HIreg)/1e21,Te,'-r','linewidth',2);
ylim([0 20]);
legend('Total Heat $\mathcal{H}$','External Heat $\mathcal{H}_E=\rho_0C_p\Theta\mathcal{V}$','Internal Heat $\mathcal{H}-\mathcal{H}_E$');
title('Heat $40^\circ$S-$15^\circ$S');
ylabel('Temperature ($^\circ$C)');
xlabel('Change in Heat Below $\Theta$ ($10^{21}$J)');




