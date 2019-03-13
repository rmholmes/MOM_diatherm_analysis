
% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

% Load Base Variables:
model = 'MOM025_kb3seg';
outputs = [86:90];
% $$$ outputs = 86;

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
% $$$ outputs = [444];

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);
regions = {'Atlantic2BAS','IndoPacific2BAS','Global'};
regLets = {'A','P',''};
for reg = 1:length(regions)
    region = regions{reg}
    regLet = regLets{reg};

%% Make Vars
type = 'ZA';
load([base model sprintf('_output%03d_',outputs(1)) region '_' type 'HBud_nomet.mat']);
eval(['names = fieldnames(' type ');']);

for i=1:length(outputs)
    load([base model sprintf('_output%03d_',outputs(i)) region '_' type 'HBud_nomet.mat']);
    if (strcmp(type,'MA'))
        ZA = MA;
        clear MA;
        yL = xL;
    end
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
    
    % NaN lats on zAI for convergence calculation:
    NANlats = zAI(:,end,1) == 0;
    zAI(NANlats,:,:,i) = NaN;

    zJSH(:,:,:,i) = zJS(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
    zPI(:,:,:,i) = zP(:,:,:,i) - zJSH(:,:,:,i);
    zN(:,:,:,i) = zdHdt(:,:,:,i) - zdVdt(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;

% $$$ zJ(:,:,:,i) = -diff(cat(1,zeros(1,TL+1,tL),zAI(:,:,:,i)),[],1);
% Note that this will not work for one grid cell of the Atlantic
% calculation for the Bering strait...  Hence to fix, moving the
% calculation of numerical mixing to after the time-averaging and
% using both Atlantic and Pacific values:
    
% $$$     zI(:,:,:,i) = zJ(:,:,:,i)-zN(:,:,:,i)-zF(:,:,:,i)-zPI(:,:,:,i)-zM(:,:,:,i);
% $$$     if (isfield(ZA,'RED'))
% $$$         zI(:,:,:,i) = zI(:,:,:,i) - zRED(:,:,:,i) - zK33(:,:,:,i);
% $$$     end
    
% $$$     zZT(:,:,i) = ZAtemp;
% $$$     zT(:,:,i) = tempZA;
% $$$     zZTx(:,:,i) = ZMtemp;
% $$$     zTx(:,:,i) = tempZM;
% $$$     zZTxa(:,:,i) = ZMAtemp;
% $$$     zTxa(:,:,i) = tempZMA;
    
    MHT(:,:,i) = zAHD(:,end,:,i);
end
months = [1:length(zP(1,1,:,1))];

% Take mean across years and months:
% $$$ zZT = mean(zZT,3);
% $$$ zT = mean(zT,3);
% $$$ zZTx = mean(zZTx,3);
% $$$ zTx = mean(zTx,3);
% $$$ zZTxa = mean(zZTxa,3);
% $$$ zTxa = mean(zTxa,3);
MHT = mean(MHT,3);
names = {names{:},'AI','AIpsi','JSH','PI','N'};
for i=1:length(names)
    eval([names{i} regLet ' = monmean(mean(z' names{i} ',4),3,ndays);']);
end

% Generate NaNst and NANsu arrays:
tmpt = monmean(mean(zF,4),3,ndays);
tmpu = monmean(mean(zPSI,4),3,ndays);
NaNst = zeros(size(tmpt));
NaNsu = zeros(size(tmpu));
for i = 1:(TL+1)
    NaNst(:,i) = tmpt(:,i) == tmpt(:,end);
    NaNsu(:,i) = tmpu(:,i) == tmpu(:,end);
end
% Max SST (not monthly, ever) line and add to NaNs (note: on Te):
maxTt = zeros(yL,1);maxTit = zeros(yL,1);
maxTu = zeros(yL,1);maxTiu = zeros(yL,1);psiTu = zeros(yL,1);
for i = 1:yL 
    indt = find(NaNst(i,:),1,'first');
    indu = find(NaNsu(i,:),1,'first');
    maxTit(i) = indt;
    maxTt(i) = Te(indt);
    NaNst(i,indt) = 0;
    maxTiu(i) = indu;
    maxTu(i) = Te(indu);
    NaNsu(i,indu) = 0;
    psiTu(i) = tmpu(i,indu);
end

% Latitude difference vector for plotting per-degree:
dy = [yu(2)-yu(1); diff(yu)]; % (First-element is done by hand -
                              % but dy is equal to second.

names = {'NaNst','NaNsu','MHT','maxTt','maxTit','psiTu','maxTu','maxTiu'};
for i=1:length(names)
    eval([names{i} regLet ' = ' names{i} ';']);
end

end %end region loop
    
%% Calculate total diathermal transport and numerical mixing (special for Atlantic):

% Global:
J = -diff(cat(1,zeros(1,TL+1),AI),[],1);
I = J-M-N-F-PI;

% Indo-Pacific:
JP = -diff(cat(1,zeros(1,TL+1),AIP),[],1);
IP = JP-MP-NP-FP-PIP;

% Atlantic:
ind = find(~isnan(AIP(:,end)),1,'last'); % 66N indicy (u-points)
JA = -diff(cat(1,zeros(1,TL+1),AIA),[],1);
JA(ind,:) = JA(ind,:) + AIP(ind-1,:); % Correct for missing bit because main Atlantic overlaps with Bering Strait cutoff
IA = JA-MA-NA-FA-PIA;

% $$$ % Check (1e-5 PW yes! without correction it's 0.1PW):
% $$$ Jd = J - JA - JP;
% $$$ Id = I - IA - IP;
% $$$ max(max(abs(Jd/1e15)))
% $$$ max(max(abs(Id/1e15)))

% $$$ %% Set above max-SST to be equal to max-SST
% $$$ % (An easy way to restrict the top temp to be bounded above by
% $$$ % max-SST):
% $$$ for reg = 1:length(regions)
% $$$     regLet = regLets{reg};
% $$$     fields = {'F','PI','M','I','N','J'};
% $$$     for i=1:length(fields)
% $$$         for yi = 1:yL
% $$$             eval([fields{i} regLet '(yi,(maxTit' regLet '(yi)+1):end) = ' fields{i} regLet ...
% $$$                   '(yi,maxTit' regLet '(yi));']);
% $$$         end
% $$$     end
% $$$     for yi=1:yL
% $$$             eval(['AI' regLet '(yi,(maxTiu' regLet '(yi)+1):end) = ' 'AI' regLet ...
% $$$                   '(yi,maxTiu' regLet '(yi));']);
% $$$     end        
% $$$ end

%% Indo-Pacific/Atlantic Temperature budgets:
ltminW = -45;
ltmin = -34;
ltmax = 90;
[tmp ltminI] = min(abs(yu-ltmin)); % u-point index
[tmp ltminWI] = min(abs(yu-ltminW)); % 45S
[tmp ltmaxI] = min(abs(yu-ltmax)); % u-point index
indBS = find(~isnan(AIP(:,end)),1,'last');
[X,Y] = ndgrid(yt,Te);

% Below a given temperature:
% Atl:
FAb = nansum(FA(ltminI+1:ltmaxI,:)+PIA(ltminI+1:ltmaxI,:),1)/1e15;
MAb = nansum(MA(ltminI+1:ltmaxI,:)+IA(ltminI+1:ltmaxI,:),1)/1e15;
NAb = nansum(NA(ltminI+1:ltmaxI,:),1)/1e15;
AA34Sb = AIA(ltminI,:)/1e15;
% Pac:
FPb = nansum(FP(ltminI+1:ltmaxI,:)+PIP(ltminI+1:ltmaxI,:),1)/1e15;
MPb = nansum(MP(ltminI+1:ltmaxI,:)+IP(ltminI+1:ltmaxI,:),1)/1e15;
NPb = nansum(NP(ltminI+1:ltmaxI,:),1)/1e15;
AP34Sb = AIP(ltminI,:)/1e15;
% Glo and BS:
A34Sb = AI(ltminI,:)/1e15;
ABSb = AIP(indBS,:)/1e15;

% SO:
FSOb = nansum(F(1:ltminI,:)+PI(1:ltminI,:),1)/1e15;
MSOb = nansum(M(1:ltminI,:)+I(1:ltminI,:),1)/1e15;
NSOb = nansum(N(1:ltminI,:),1)/1e15;

% $$$ % SO and warm route:
FWRb = nansum(F((ltminWI+1):ltminI,:)+PI((ltminWI+1):ltminI,:),1)/1e15;
MWRb = nansum(M((ltminWI+1):ltminI,:)+I((ltminWI+1):ltminI,:),1)/1e15;
NWRb = nansum(N((ltminWI+1):ltminI,:),1)/1e15;

A45Sb = AI(ltminWI,:)/1e15;
% $$$ FSOb = nansum(F(1:ltminWI,:)+PI(1:ltminWI,:),1)/1e15;
% $$$ MSOb = nansum(M(1:ltminWI,:)+I(1:ltminWI,:),1)/1e15;
% $$$ NSOb = nansum(N(1:ltminWI,:),1)/1e15;

% $$$ % Do integral the other way:
% $$$ vars = {'FAb','MAb','NAb','AA34Sb','FPb','MPb','NPb','AP34Sb', ...
% $$$         'A34Sb','ABSb'};
% $$$ for i = 1:length(vars)
% $$$     eval([vars{i} ' = -' vars{i} ' + ' vars{i} '(end);']);
% $$$ end
% Note: Diathermal fluxes are positive upwards, AI's are positive northwards

% Summary numbers at temps:
%sT = [-3 5 10 15 20 34];
sT = [-3 20 34];
for Ti=(length(sT)-1):-1:1
    [tmp indL] = min(abs(Te-sT(Ti)));
    [tmp indU] = min(abs(Te-sT(Ti+1)));
    
    disp(' ')
    disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
    disp(['-------------------------------------------------------------'])
    disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL),FWRb(indU)-FWRb(indL)))
    disp(sprintf('Mixing      (flux at bottom) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',MPb(indL),MAb(indL),MSOb(indL),MWRb(indL)))
    disp(sprintf('Transport 34S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
    disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
    disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
    disp(sprintf('Transport BS    (into layer)          = %5.2f',ABSb(indU)-ABSb(indL)))
end

% Global:
indL = 1;indU=TL+1;
for i=1:1
    disp(' ')
    disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
    disp(['-------------------------------------------------------------'])
    disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL),FWRb(indU)-FWRb(indL)))
    disp(sprintf('Mixing      (flux at bottom) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',MPb(indL),MAb(indL),MSOb(indL),MWRb(indL)))
    disp(sprintf('Transport 32S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
    disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
    disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
    disp(sprintf('Transport BS    (into layer)          = %5.2f',ABSb(indU)-ABSb(indL)))
end

% All temperatures figure:
figure;
% Atlantic Budget:
subplot(3,1,1);
hold on;
grid on;
ylabel('Heat Flux (PW)');
set(gca,'xticklabel',[]);
box on;
xlim([-2 34]);
ylim([-0.2 0.9]);
plot(Te,FAb,'-k','linewidth',2);
plot(Te,-MAb,'-r','linewidth',2);
plot(Te,-NAb,'-m','linewidth',2);
plot(Te,AA34Sb,'-b','linewidth',2);
plot(Te,ABSb,'-.b','linewidth',2);
plot(Te,AA34Sb+ABSb-MAb-NAb,'--k');
legend('$\mathcal{F}_A$','-$\mathcal{M}_A$','-$\frac{\partial\mathcal{H}_I}{\partial t}_A$','$\mathcal{A}_A^{34^\circ S}$','$\mathcal{A}^{BS}$');%,'$\mathcal{F}_A^{NH}$');
% $$$ title('Atlantic Budget (PW)');
set(gca,'FontSize',15);
set(gca,'Position',[0.1300    0.6894    0.7750    0.2665]);

% 34S Transports:
subplot(3,1,2);
hold on;
grid on;
ylabel('Heat Flux (PW)');
set(gca,'xticklabel',[]);
box on;
xlim([-2 34]);
ylim([-1 1]);
plot(Te,A34Sb,'--b','linewidth',2);
plot(Te,AA34Sb,'-b','linewidth',2);
plot(Te,-AP34Sb,':b','linewidth',2);
legend('Total $\mathcal{A}_G^{34^\circ S}$','Atlantic $\mathcal{A}_A^{34^\circ S}$','Pacific $-\mathcal{A}_P^{34^\circ S}$');
set(gca,'FontSize',15);
set(gca,'Position',[0.1300    0.41    0.7750    0.2665]);

% Pacific budget:
subplot(3,1,3);
hold on;
grid on;
ylabel('Heat Flux (PW)');
xlabel('Temperature ($^\circ$C)');
box on;
xlim([-2 34]);
ylim([-0.4 1.2]);
plot(Te,-AP34Sb,':b','linewidth',2);
plot(Te,-FPb,'-k','linewidth',2);
plot(Te,-MPb,'-r','linewidth',2);
plot(Te,-NPb,'-m','linewidth',2);
plot(Te,-ABSb,'-.b','linewidth',2);
plot(Te,-(FPb+MPb+NPb+ABSb),'--k','linewidth',2);
legend('Pacific $-\mathcal{A}_P^{34^\circ S}$','$-\mathcal{F}_P$','$-\mathcal{M}_P$','$-\frac{\partial\mathcal{H}_I}{\partial t}_P$','$-\mathcal{A}^{BS}$');
set(gca,'FontSize',15);
set(gca,'Position',[0.1300    0.11    0.7750    0.2665]);

figure;
% Southern Ocean Budget:
subplot(3,1,1);
hold on;
grid on;
xlabel('Temperature ($^\circ$C)');
ylabel('Heat Flux (PW)');
box on;
xlim([-2 34]);
ylim([-0.5 0.5]);
plot(Te,A34Sb,'--b','linewidth',2);
plot(Te,-FSOb,'-k','linewidth',2);
plot(Te,-MSOb,'-r','linewidth',2);
plot(Te,-NSOb,'-m','linewidth',2);
% $$$ plot(Te,-FSOb-MSOb-NSOb,'--k');
legend('$\mathcal{A}_G^{34^\circ S}$','$-\mathcal{F}_{SO}$','$-\mathcal{M}_{SO}$','$-\frac{\partial\mathcal{H}_I}{\partial t}_{SO}$');%,'$\mathcal{F}_A^{NH}$');
% $$$ title('Atlantic Budget (PW)');
set(gca,'FontSize',15);
set(gca,'Position',[0.1300    0.6894    0.7750    0.2665]);

%% Schematic plot:
ThetaB = 20;
figure;
set(gcf,'Position',[2883          98         560         905]);
subplot(2,1,1);
tmp = AIA/1e15;
tmp(NaNsuA==1) = NaN;
[Xu,Yu] = ndgrid(yu,Te);
contourf(Xu,Yu,tmp,[-3:0.1:3],'linestyle','none');
hold on;
[c,h] = contour(Xu,Yu,tmp,[0.2:0.2:3],'-k');
clabel(c,h);
[c,h] = contour(Xu,Yu,tmp,[-3:0.2:-0.2],'--k');
clabel(c,h);
plot(yu,maxTuA,':k');
plot(ltmin*[1 1],[Te(1) Te(end)],'--k');
ylabel('Temperature ($^\circ$C)');
xlabel('Latitude ($^\circ$N)');
cb = colorbar;
ylabel(cb,'PW');
colormap(redblue);
plot([-34 -34 max(yu) max(yu) -34],[Te(1) ThetaB ThetaB Te(1) Te(1)],'-m','linewidth',2);
xlim([-45 90]);
caxis([-1.5 1.5]);
ylim([-2 34]);
title('Atlantic');
set(gca,'FontSize',15);
subplot(2,1,2);
tmp = AIP/1e15;
tmp(NaNsuP==1) = NaN;
[Xu,Yu] = ndgrid(yu,Te);
contourf(Xu,Yu,tmp,[-3:0.1:3],'linestyle','none');
hold on;
[c,h] = contour(Xu,Yu,tmp,[0.2:0.2:3],'-k');
clabel(c,h);
[c,h] = contour(Xu,Yu,tmp,[-3:0.2:-0.2],'--k');
clabel(c,h);
plot(yu,maxTuP,':k');
plot(ltmin*[1 1],[Te(1) Te(end)],'--k');
ylabel('Temperature ($^\circ$C)');
xlabel('Latitude ($^\circ$N)');
cb = colorbar;
ylabel(cb,'PW');
colormap(redblue);
plot([-34 -34 66 66 -34],[Te(1) ThetaB ThetaB Te(1) Te(1)],'-m','linewidth',2);
xlim([-45 90]);
ylim([-2 34]);
caxis([-1.5 1.5]);
title('Indo-Pacific');
set(gca,'FontSize',15);



% Region:
[tmp Bind] = min(abs(Te-ThetaB));
Xreg = [max(max(X)) ltmin ltmin];
Yreg = [Te(1) Te(1) ThetaB];
cnt = 4;
for yi=1:yL
    if (yi>=ltminI & maxTi(yi) <= Bind)
        Xreg = [Xreg yt(yi)];
        Yreg = [Yreg maxT(yi)];
    end
end

Xf = [
mask = NaN*zeros(size(AI));
for yi=1:yL
    if (yi>=ltminI)
        mask(yi,1:maxTi(yi)) = 1;
    end
end
mask(:,Bind+1:end) = NaN;



%% Plot lat-temp

% Load SST for plotting:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'SST');
SSTa = SST;
ofor i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))],'SST');
    SSTa = SSTa+SST;
end
SST = SSTa/length(outputs);
if (~(strcmp(region,'') | strcmp(region,'Global')))
    [maskREG,~,~,~,~,~,~] = Heat_Budget_Mask(region(1:end-1),'','','',base, ...
                                             model);
else
    maskREG = ones(size(SST(:,:,1)));
end
SST = SST.*repmat(maskREG,[1 1 tL]);
SST(SST==0) = NaN;
if (strcmp(type,'ZA'))
    meanSST = squeeze(nanmean(monmean(SST,3,ndays),1));
    minSST = squeeze(min(monmean(SST,3,ndays),[],1));
% $$$ maxSST = squeeze(max(monmean(SST,3,ndays),[],1));
else
    meanSST = squeeze(nanmean(monmean(SST,3,ndays),2));
    minSST = squeeze(min(monmean(SST,3,ndays),[],2));
    maxSST = squeeze(max(monmean(SST,3,ndays),[],2));
end    

% Plot Streamfunction and Heat Function:
fields = { ...
% $$$           {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-30 30],1,'Sv'}, ...
% $$$           {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-1.5 1.5],0.05,'PW'}, ...
          {zPSI(:,:,months)/1e6, 'Streamfunction $\Psi$',[-180 180],2,'Sv'}, ...
          {zAI(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$',[-20 20],0.05,'PW'}, ...
% $$$           {zAIpsi(:,:,months)/1e15, 'Heat Function $\mathcal{A}_I$ from $\Psi$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {zAHD(:,:,months)/1e15, 'Heat Function $\mathcal{A}$',[-1.5 1.5],0.1,'PW'}, ...
% $$$           {(zAI(:,:,months)-zAIpsi(:,:,months))/1e15, 'Heat Function diff $\mathcal{A}$',[-0.01 0.01],0.005,'PW'}, ...
};

% Plot overall fields:
fields = { ...
% $$$           {-rho0*Cp*zPSI(:,:,months)/1e12, 'Isothermal Heat Flux $\mathcal{J}_{iso} = -\rho_0C_p\Psi$',[-150 150],15,'TW / $^\circ$C'}, ...
% $$$           {zJ(:,:,months)/1e12, 'Diathermal Heat Flux $\mathcal{J}_{dia} = -\frac{\partial\mathcal{A}_I}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zN(:,:,months)/1e12, 'Tendency $\frac{\partial}{\partial\phi}\frac{\partial\mathcal{H}_I}{\partial t}$',[-10 10],0.5,'TW / $^\circ$'}, ...
% $$$           {(zF(:,:,months)+zPI(:,:,months))/1e12, 'Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
% $$$           {zM(:,:,months)/1e12+zI(:,:,months)/1e12, 'Mixing $\frac{\partial(\mathcal{M}+\mathcal{I})}{\partial\phi}$',[-50 50],5,'TW / $^\circ$'}, ...
          {(zF(:,:,months)+zPI(:,:,months))/1e12, 'Surface Forcing $\frac{\partial\left(\mathcal{F}+\mathcal{P}_I\right)}{\partial\phi}$',[-20 20],0.5,'TW / $^\circ$'}, ...
% $$$           {zM(:,:,months)/1e12+zI(:,:,months)/1e12, 'Mixing $\frac{\partial(\mathcal{M}+\mathcal{I})}{\partial\phi}$',[-20 20],0.5,'TW / $^\circ$'}, ...
          {zM(:,:,months)/1e12, 'Vertical Mixing $\frac{\partial(\mathcal{M}+\mathcal{I})}{\partial\phi}$',[-10 10],0.25,'TW / $^\circ$'}, ...
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
% $$$ tmp(yt>60,:,:) = 1e15;
% $$$ zinds = zeros(TL+1,1);
% $$$ lats = zeros(TL+1,1);
% $$$ for i=1:TL+1
% $$$     ind = find(tmp(:,i)<0,1,'last');
% $$$     if (strcmp(region,'IndoPacificNZ_'))
% $$$         ind = find(tmp(:,i)>0,1,'first');
% $$$     elseif (strcmp(region,'AtlanticNZ_'))
% $$$         ind = find(tmp(:,i)>0,1,'first')+1;
% $$$     end
% $$$     if (length(ind)>0)
% $$$         zinds(i) = ind;
% $$$     else
% $$$         zinds(i) = 1;
% $$$     end
% $$$     lats(i) = yt(zinds(i));
% $$$ end
% $$$ for i=1:length(fields)
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
clab = [1 1 1 1 1 1];
% $$$ clab = [1 0 0 0 0 0];

cmap = redblue(npts-3);
% $$$ cmap = parula(npts-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

% $$$ AIsp = 0.25;
AIsp = 1;

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

if (strcmp(type,'ZA'))
    yvec = yt';
else
    yvec = lon(:,1);
end
 
for i=1:length(fields)
% $$$     if (i >=3)
% $$$         subplot(2,3,i+1);
% $$$     else
        subplot(2,3,i);
% $$$         subplot(3,1,i);
% $$$     end
    if (length(fields{i}{1}(1,:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    if (doZAremap)
        Yg = repmat(yvec,[1 TL+1]);
        Tg = zZTxa;
    else
        [Yg,Tg] = ndgrid(yvec,x);
    end

    VAR = monmean(fields{i}{1},3,ndays(months));
    VAR(VAR==0) = NaN;
    VAR(NaNs==1) = NaN;
    VAR = filter_field(VAR',latfilt,'-t')';
    contourf(Yg,Tg,VAR,cpts{i},'linestyle','none');
    hold on;
    col = [0 0 0];

    if (~doZAremap)
% $$$         plot(yvec,filter_field(meanSST,latfilt,'-t'),':','color',col);
        plot(yvec,filter_field(minSST,latfilt,'-t'),':','color',col);
    plot(yvec,filter_field(maxT,latfilt,'-t'),':k');
% $$$     plot(yvec,21.5*ones(size(yvec)),'--k');
% $$$         plot(lats,Te,'-k','linewidth',2);%'color',[0 0.5 0]);

% $$$     tmp = filter_field(monmean(zPSI,3,ndays)'/1e6,latfilt,'-t')';
% $$$     tmp(Yg<20 & (Tg > 4.9 & Tg < 8.75)) = NaN;
% $$$     contour(Yg,Tg,tmp,[-14.5 -14.5],'-','color',[0 0.5 0],'linewidth',2);
% $$$     contour(Yg,Tg,tmp,[5 5],'--','color',[0 0.5 0],'linewidth',2);
    else
        [tt,zz] = ndgrid(yvec,-z);
        [c,h] = contour(tt,zz,zTxa,[0:4:34],'-k');
        clabel(c,h);
    end

% $$$     if (i==1)
    tmp = monmean(zAI,3,ndays)/1e15;
    tmp(NaNs==1) = NaN;
    tmp = filter_field(tmp',latfilt,'-t')';
    [c,h] = contour(Yg,Tg,tmp,[-50:AIsp:-AIsp],'--k');
    if (clab(i))
        clabel(c,h);
    end
    [c,h] =  contour(Yg,Tg,tmp,[AIsp:AIsp:50],'-k');
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
    caxis(fields{i}{3});
    box on; 
    grid on;
    if (strcmp(type,'ZA'))
        xlim([-70 70]);
        xlabel('Latitude ($^\circ$N)');
    else
        xlim([-280 80]);
        xlabel('Longitude ($^\circ$E)');
    end
    cb = colorbar;
    ylabel(cb,fields{i}{5});
    title([fields{i}{2}]);
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
tmp(yt>60,:,:) = 1e15;
zinds = zeros(TL+1,1);
lats = zeros(TL+1,1);
for i=1:TL+1
    ind = find(tmp(:,i)<0,1,'last');
    if (strcmp(region,'IndoPacificNZ_'))
        ind = find(tmp(:,i)>0,1,'first');
    elseif (strcmp(region,'AtlanticNZ_'))
        ind = find(tmp(:,i)>0,1,'first')+1;
    end
    if (length(ind)>0)
        zinds(i) = ind;
    else
        zinds(i) = 1;
    end
    lats(i) = yt(zinds(i));
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
plot(yt,filter_field(monmean(MHT,2,ndays)/1e15,lfilt,'-t'),'--k','linewidth',3);
hold on; 
plot(yt,filter_field(rho0*Cp*psiT.*maxT,lfilt,'-t')/1e15,':k','linewidth',3);
for i=1:4
    plot(yt,filter_field(MHTtot{i,1},lfilt,'-t'),colors{i},'linewidth',3);
end
% $$$ plot(yt,filter_field(MHTtot{2,1}+MHTtot{3,1}+MHTtot{4,1},1,'-t'),'-c');
% $$$ plot(yt,filter_field(MHTtot{1,1},1,'-t')+rho0*Cp*psiT.*maxT/1e15,'-y');
legend({'Total Heat Transport $\mathcal{A}^{max}$','External Heat Transport $\rho_0C_p\Psi^{max}\Theta^{max}$',MHTtot{:,2},});
xlabel('Latitude ($^\circ$N)');
ylabel('PW');
xlim([-90 90]);
grid on;

% Add observations:
obsfile = '../Observations/ANNUAL_TRANSPORTS_1985_1989.ascii.txt';
data = importdata(obsfile);
obs_lat = data.data(:,1)/100;
obs_total_ncep = data.data(:,7)/100;
obs_atl_ncep = data.data(:,8)/100;
obs_pac_ncep = data.data(:,9)/100;
obs_ind_ncep = data.data(:,10)/100;
obs_total_ecmwf = data.data(:,15)/100;
obs_atl_ecmwf = data.data(:,16)/100;
obs_pac_ecmwf = data.data(:,17)/100;
obs_ind_ecmwf = data.data(:,18)/100;

obs_total_ncep(abs(obs_total_ncep)>5) = 0;
obs_atl_ncep(abs(obs_atl_ncep)>5) = 0;
obs_pac_ncep(abs(obs_pac_ncep)>5) = 0;
obs_ind_ncep(abs(obs_ind_ncep)>5) = 0;
obs_total_ecmwf(abs(obs_total_ecmwf)>5) = 0;
obs_atl_ecmwf(abs(obs_atl_ecmwf)>5) = 0;
obs_pac_ecmwf(abs(obs_pac_ecmwf)>5) = 0;
obs_ind_ecmwf(abs(obs_ind_ecmwf)>5) = 0;

hold on;
plot(obs_lat,obs_total_ncep,'xk');
plot(obs_lat,obs_total_ecmwf,'ok');
plot(obs_lat,obs_atl_ncep,'xr');
plot(obs_lat,obs_pac_ncep+obs_ind_ncep,'xb');
plot(obs_lat,obs_atl_ecmwf,'or');
plot(obs_lat,obs_pac_ecmwf+obs_ind_ecmwf,'ob');
% $$$ plot(obs_lat,data.data(:,5)/100,'om');

obsfile = '../Observations/GW2000_MeridionalHeatTransportData.mat';
load(obsfile);
hold on;
errorbar(GLOlat,GLOHT,GLOHTe,'dk','linewidth',2,'MarkerSize',10);
errorbar(IPAClat,IPACHT,IPACHTe,'db','linewidth',2,'MarkerSize',10);
errorbar(ATLlat,ATLHT,ATLHTe,'dr','linewidth',2,'MarkerSize',10);





%% Calculate heat transports in different cells:

[Yg,Tg] = ndgrid(yt,Te);
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
    maxT = NaN*zeros(size(yt));
    minT = NaN*zeros(size(yt));
    for k=1:length(yt)
        if (length(find(regions{i}{3}(k,:),1,'last')))
            maxTi = find(regions{i}{3}(k,:),1,'last');
            minTi = find(regions{i}{3}(k,:),1,'first');
            maxT(k) = interp1(psi(k,maxTi:maxTi+1),Te(maxTi:maxTi+1),regions{i}{2},'linear');
            minT(k) = interp1(psi(k,minTi-1:minTi),Te(minTi-1:minTi),regions{i}{2},'linear');
        end
    end

    for ii=1:length(fields)
        OHT{i,ii} = NaN*zeros(size(yt));
        extT = NaN*zeros(size(yt));
        VAR = monmean(fields{ii}{1},3,ndays);
        
        for k=1:length(yt)
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
total = zeros(size(yt));
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
lgh(1) = plot(yt,filter_field(mhflux,latfilt,'-t'),'-k','linewidth',3);
hold on;
lg{1} = 'Total';
cols{1} = 'k';
for i=1:length(regions)
    lgh(i+1) = plot(yt,Repl(filter_field(OHT{i,1},latfilt,'-t'),NaN,0),regions{i}{4},'color',regions{i}{6},'linewidth',3);
    lg{i+1} = regions{i}{1};
    cols{i+1} = regions{i}{6};
end
hold on;
lgh(i+2) = plot(yt,filter_field(mixed,latfilt,'-t'),'-','color',[0.6 ...
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
[Y,Z] = ndgrid(yt,-z);
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

reg = yt>-40 & yt <-15;
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




