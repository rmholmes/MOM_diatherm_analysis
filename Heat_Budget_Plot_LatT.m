

% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

% Load Base Variables:
model = 'MOM025_kb3seg';
outputs = [101:110];

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
end
% Latitude difference vector for plotting per-degree:
dy = [yu(2)-yu(1); diff(yu)]; % (First-element is done by hand - but dy is equal to second).

regions = {'Atlantic2BAS','IndoPacific2BAS','Global'};%,'SO_Atlantic','SO_IndoPacific'};
regLets = {'A','P','G'};%,'SA','SP'};

for reg = 1:length(regions)
    region = regions{reg}
    regLet = regLets{reg};

    %% Make Vars
    type = 'ZA';
    load([base model sprintf('_output%03d_',outputs(1)) region '_' type 'HBud.mat']);
    eval(['names = fieldnames(' type ');']); % names in this file
    namesEX = {'K33','RED','AHDR','GM','PSIGM','AHDGM','MDS','SIG'}; % extra names that might not exist
    namesALL = {names{:},namesEX{:}};
            
    for i=1:length(outputs)
        i
        load([base model sprintf('_output%03d_',outputs(i)) region '_' type 'HBud.mat']);
        for ii=1:length(names)
            eval(['ZAR.' names{ii} '(:,:,:,i) = ZA.' names{ii} ';']);
        end
        for ii=1:length(namesEX) % add extras as zeros if they don't exist
            if (~isfield(ZA,namesEX{ii}))
                eval(['ZAR.' namesEX{ii} '(:,:,:,i) = zeros(size(ZAR.F(:,:,:,i)));']);
            end
        end
        clear ZA;
    
        % Heat Function:
        ZAR.AIadv(:,:,:,i) = ZAR.AHD(:,:,:,i) - rho0*Cp*repmat(Te',[yL 1 tL]).*ZAR.PSI(:,:,:,i); % advective-AI (defined on Te and v-points)
% $$$     zAIpsi(:,:,:,i) = -rho0*Cp*cumsum(zPSI(:,:,:,i)*dT,2); % Heat Function (defined on Tc, since Psi is on Te, and defined on v-points)
        ZAR.AI(:,:,:,i) = ZAR.AIadv(:,:,:,i)+ZAR.AHDGM(:,:,:,i)+ZAR.AHDSUB(:,:,:,i)+ZAR.AHDR(:,:,:,i); % total internal heat content transport

        % NaN lats on zAI for convergence calculation (for sub-regions):
        NANlats = ZAR.AI(:,end,1) == 0;
        ZAR.AI(NANlats,:,:,i) = NaN;

        ZAR.JSH(:,:,:,i) = ZAR.JS(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
        ZAR.PI(:,:,:,i)  = ZAR.P(:,:,:,i) - ZAR.JSH(:,:,:,i);
        ZAR.N(:,:,:,i) = ZAR.dHdt(:,:,:,i) - ZAR.dVdt(:,:,:,i).*repmat(Te',[yL 1 tL])*rho0*Cp;
    end

    % Take annual and mean across outputs:
    names = fieldnames(ZAR);
    for ii=1:length(names)
        eval(['ZAR.' names{ii} ' = monmean(mean(ZAR.' names{ii} ',4),3,ndays);']);
    end

    % Generate NaNst and NANsu arrays:
    ZAR.NaNst = zeros(size(ZAR.F));
    ZAR.NaNsu = zeros(size(ZAR.PSI));
    for i = 1:(TL+1)
        ZAR.NaNst(:,i) = ZAR.F(:,i) == ZAR.F(:,end);
        ZAR.NaNsu(:,i) = ZAR.PSI(:,i) == ZAR.PSI(:,end);
    end

    % Max SST (not monthly, ever) line and add to NaNs (note: on Te):
    ZAR.maxTt = zeros(yL,1);ZAR.maxTit = zeros(yL,1);
    ZAR.maxTu = zeros(yL,1);ZAR.maxTiu = zeros(yL,1);ZAR.psiTu = zeros(yL,1);
    for i = 1:yL 
        indt = find(ZAR.NaNst(i,:),1,'first');
        indu = find(ZAR.NaNsu(i,:),1,'first');
        ZAR.maxTit(i) = indt;
        ZAR.maxTt(i) = Te(indt);
        ZAR.NaNst(i,indt) = 0;
        ZAR.maxTiu(i) = indu;
        ZAR.maxTu(i) = Te(indu);
        ZAR.NaNsu(i,indu) = 0;
        ZAR.psiTu(i) = ZAR.PSI(i,indu);
    end

    % Min SST (monthly - from surface forcing):
    % Load SST for plotting:
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'SST');
    SSTa = SST;
    for i=2:length(outputs)
        load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))],'SST');
        SSTa = SSTa+SST;
    end
    SST = SSTa/length(outputs);
    if (~(strcmp(region,'') | strcmp(region,'Global')))
        [maskREG,~,~,~,~,~,~] = Heat_Budget_Mask(region,'','',base,model);
    else
        maskREG = ones(size(SST(:,:,1)));
    end
    SST = SST.*repmat(maskREG,[1 1 tL]);
    SST(SST==0) = NaN;
    ZAR.minSST = squeeze(min(monmean(SST,3,ndays),[],1)');
    if (max(ZAR.minSST)>100); ZAR.minSST = ZAR.minSST-273.15; end

    % Total MHTs:
    ZAR.MHTSUB = ZAR.AHDSUB(:,end);
    ZAR.MHTADV = ZAR.AI(:,end);
    ZAR.MHTGM  = ZAR.AHDGM(:,end);
    ZAR.MHTR   = ZAR.AHDR(:,end);
    ZAR.MHT    = ZAR.MHTADV + ZAR.MHTSUB + ZAR.MHTGM + ZAR.MHTR;

    eval(['ZA_' regLet ' = ZAR;']);
    clear ZAR;
end %end region loop
    
%% Calculate total diathermal transport and numerical mixing (special for Atlantic):
[X,Y] = ndgrid(yt,Te);

% Add KPPNL for old processing files:
if (~isfield(ZA_G,'KPPNL'))
    ZA_G.KPPNL = 0*ZA_G.M;
end
if (~isfield(ZA_P,'KPPNL'))
    ZA_P.KPPNL = 0*ZA_P.M;
end
if (~isfield(ZA_A,'KPPNL'))
    ZA_A.KPPNL = 0*ZA_A.M;
end

% Global:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_G.AI-ZA_G.AHDR),[],1); % convergence of that transport
ZA_G.Jdia = -ZA_G.N-dAI_mR_dphi; % total diathermal transport
ZA_G.I = -(ZA_G.Jdia+ZA_G.M+ZA_G.KPPNL+ZA_G.F+ZA_G.PI+ZA_G.RED+ZA_G.K33+ZA_G.MDS+ZA_G.SIG); % numerical mixing (both advective and submesoscale)
if (isfield(ZA_G,'NUM_SUBlf'))
    ZA_G.NUM = ZA_G.NUM_SUBlf;
end

% Indo-Pacific:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_P.AI-ZA_P.AHDR),[],1); % convergence of that transport
ZA_P.Jdia = -ZA_P.N-dAI_mR_dphi; % total diathermal transport
ZA_P.I = -(ZA_P.Jdia+ZA_P.M+ZA_P.KPPNL+ZA_P.F+ZA_P.PI+ZA_P.RED+ZA_P.K33+ZA_P.MDS+ZA_P.SIG); % numerical mixing (both advective and submesoscale)
if (isfield(ZA_P,'NUM_SUBlf'))
    ZA_P.NUM = ZA_P.NUM_SUBlf;
end

% Atlantic:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_A.AI-ZA_A.AHDR),[],1); % convergence of that transport
ZA_A.Jdia = -ZA_A.N-dAI_mR_dphi; % total diathermal transport
ind = find(~isnan(ZA_P.AI(:,end)),1,'last'); % 66N indicy (u-points)
ZA_A.Jdia(ind,:) = ZA_A.Jdia(ind,:) + ZA_P.AI(ind-1,:); % Correct for missing bit because main Atlantic overlaps with Bering Strait cutoff
ZA_A.I = -(ZA_A.Jdia+ZA_A.M+ZA_A.KPPNL+ZA_A.F+ZA_A.PI+ZA_A.RED+ZA_A.K33+ZA_A.MDS+ZA_A.SIG); % numerical mixing (both advective and submesoscale)
if (isfield(ZA_A,'NUM_SUBlf'))
    ZA_A.NUM = ZA_A.NUM_SUBlf;
end

if (exist('ZA_SA'))
% SO Atlantic:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_SA.AI-ZA_SA.AHDR),[],1); % convergence of that transport
ZA_SA.Jdia = -ZA_SA.N-dAI_mR_dphi; % total diathermal transport
ZA_SA.I = -(ZA_SA.Jdia+ZA_SA.M+ZA_SA.KPPNL+ZA_SA.F+ZA_SA.PI+ZA_SA.RED+ZA_SA.K33+ZA_SA.MDS+ZA_SA.SIG); % numerical mixing (both advective and submesoscale)
end

if (exist('ZA_SP'))
% SO IndoPacific:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_SP.AI-ZA_SP.AHDR),[],1); % convergence of that transport
ZA_SP.Jdia = -ZA_SP.N-dAI_mR_dphi; % total diathermal transport
ZA_SP.I = -(ZA_SP.Jdia+ZA_SP.M+ZA_SP.KPPNL+ZA_SP.F+ZA_SP.PI+ZA_SP.RED+ZA_SP.K33+ZA_SP.MDS+ZA_SP.SIG); % numerical mixing (both advective and submesoscale)
end

groups = {'G','P','A'}%,'SA','SP'};
for gi=1:length(groups)
    eval(['ZA_' groups{gi} '.AIF = 0*ZA_G.F;']); % Heat lost max SST line from AI

    % Fix above max SST to calculate budgets between theta or max SST
    % and -2C:
    for yi=1:yL
        %PI:
        eval(['ZA_' groups{gi} '.PI(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
              '= ZA_' groups{gi} '.PI(yi,ZA_' groups{gi} '.maxTit(yi));']);
        %tendency:
        eval(['ZA_' groups{gi} '.N(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
              '= ZA_' groups{gi} '.N(yi,ZA_' groups{gi} '.maxTit(yi));']);
        %Num-mix:
        eval(['ZA_' groups{gi} '.I(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
              '= ZA_' groups{gi} '.I(yi,ZA_' groups{gi} '.maxTit(yi));']);
        %AI:
        eval(['ZA_' groups{gi} '.AI(yi,ZA_' groups{gi} '.maxTiu(yi):end) ' ...
              '= ZA_' groups{gi} '.AI(yi,ZA_' groups{gi} '.maxTiu(yi));']);
    end
    
    % Correct surface forcing by residual:
    eval(['ZA_' groups{gi} '.AIF = diff(cat(1,zeros(1,TL+1),ZA_' groups{gi} '.AI-ZA_' ...
          groups{gi} '.AHDR),[],1)-ZA_' groups{gi} '.PI-ZA_' groups{gi} '.F-ZA_' ...
          groups{gi} '.M-ZA_' groups{gi} '.I-ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.N;']);

    % Fix Atlantic Bering Strait part:
    if (strcmp(groups{gi},'A'))
        ind = find(~isnan(ZA_P.AI(:,end)),1,'last'); % 66N indicy (u-points)
        eval(['ZA_' groups{gi} '.AIF(ind,:) = ZA_' groups{gi} '.AIF(ind,:)-ZA_P.AI(ind-1,:);']);
    end
    
    eval(['ZA_' groups{gi} '.Fall = ZA_' groups{gi} '.F + ZA_' groups{gi} '.PI+ZA_' groups{gi} '.AIF;']);
    eval(['ZA_' groups{gi} '.Mall = ZA_' groups{gi} '.M + ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.I;']);
end

%%%% Summary schematics:
ltminW = -45;
ltmin = -34;
ltmax = 90;
[tmp ltminI] = min(abs(yu-ltmin)); % u-point index
[tmp ltminWI] = min(abs(yu-ltminW)); % 45S
[tmp ltmaxI] = min(abs(yu-ltmax)); % u-point index
[tmp lt50I] = min(abs(yu-50)); % u-point index
indBS = find(~isnan(ZA_P.AI(:,end)),1,'last');
[X,Y] = ndgrid(yt,Te);

% Below a given temperature:
% Atl:
FAb = nansum(ZA_A.Fall(ltminI+1:ltmaxI,:),1)/1e15;
MAb = nansum(ZA_A.Mall(ltminI+1:ltmaxI,:),1)/1e15;
NAb = nansum(ZA_A.N(ltminI+1:ltmaxI,:),1)/1e15;
AA34Sb = ZA_A.AI(ltminI,:)/1e15;
AA50Nb = ZA_A.AI(lt50I,:)/1e15;
% Pac:
FPb = nansum(ZA_P.Fall(ltminI+1:ltmaxI,:),1)/1e15;
MPb = nansum(ZA_P.Mall(ltminI+1:ltmaxI,:),1)/1e15;
NPb = nansum(ZA_P.N(ltminI+1:ltmaxI,:),1)/1e15;
AP34Sb = ZA_P.AI(ltminI,:)/1e15;
% Glo and BS:
A34Sb = ZA_G.AI(ltminI,:)/1e15;
ABSPb = ZA_P.AI(indBS,:)/1e15;

% SO:
FSOb = nansum(ZA_G.Fall(1:ltminI,:),1)/1e15;
MSOb = nansum(ZA_G.Mall(1:ltminI,:),1)/1e15;
NSOb = nansum(ZA_G.N(1:ltminI,:),1)/1e15;

% SO and warm route:
FWRb = nansum(ZA_G.Fall((ltminWI+1):ltminI,:),1)/1e15;
MWRb = nansum(ZA_G.Mall((ltminWI+1):ltminI,:),1)/1e15;
NWRb = nansum(ZA_G.N((ltminWI+1):ltminI,:),1)/1e15;

if (exist('ZA_SA'))
% SO Atlantic and SO Pacific surface forcing:
FSAb = nansum(ZA_SA.Fall(1:ltminI,:),1)/1e15;
FSPb = nansum(ZA_SP.Fall(1:ltminI,:),1)/1e15;
FWRAb = nansum(ZA_SA.Fall((ltminWI+1):ltminI,:),1)/1e15;
FWRPb = nansum(ZA_SP.Fall((ltminWI+1):ltminI,:),1)/1e15;
else
FSAb = zeros(size(Te));
FSPb = FSAb;FWRAb = FSAb;FWRPb=FSAb;    
end

A45Sb = ZA_G.AI(ltminWI,:)/1e15;

% Residuals by basin:
RESP = FPb+MPb+AP34Sb-NPb-ABSPb;
RESA = FAb+MAb+AA34Sb-NAb+ABSPb;
RESS = FSOb+MSOb-A34Sb-NSOb; % All checks out.

% Summary numbers at temps:
%sT = [-3 5 10 15 20 34];
sT = [-3 15 20 34];
for Ti=(length(sT)-1):-1:1
    [tmp indL] = min(abs(Te-sT(Ti)));
    [tmp indU] = min(abs(Te-sT(Ti+1)));
    
    disp(' ')
    disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
    disp(['-------------------------------------------------------------'])
    disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f (Atl = %5.2f, Pac = %5.2f), WR = %5.2f (Atl = %5.2f, Pac = %5.2f)',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL),FSAb(indU)-FSAb(indL),FSPb(indU)-FSPb(indL),FWRb(indU)-FWRb(indL),FWRAb(indU)-FWRAb(indL),FWRPb(indU)-FWRPb(indL)))
    disp(sprintf('Mixing      (flux at bottom) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',MPb(indL),MAb(indL),MSOb(indL),MWRb(indL)))
    disp(sprintf('Transport 34S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
    disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
    disp(sprintf('Transport 50N   (into layer) Atlantic = %5.2f',AA50Nb(indU)-AA50Nb(indL)))
    disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
    disp(sprintf('Transport BS    (into layer)            %5.2f.',ABSPb(indU)-ABSPb(indL)))
    disp(sprintf('Residuals                    Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',RESP(indU)-RESP(indL),RESA(indU)-RESA(indL),RESS(indU)-RESS(indL)))
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
    disp(sprintf('Transport 50N   (into layer) Atlantic = %5.2f',AA50Nb(indU)-AA50Nb(indL)))
    disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
    disp(sprintf('Transport BS    (into layer)          = %5.2f',ABSPb(indU)-ABSPb(indL)))
    disp(sprintf('Residuals                    Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',RESP(indU)-RESP(indL),RESA(indU)-RESA(indL),RESS(indU)-RESS(indL)))
end

%% Plot latitude - temperature plane for different basins:
dy = diff(yu);
dy = [dy(1); dy];

% PsiAI nice plot:
fields = { ...
          {'PSI',1/1e6,'$\Psi$',[-30 30],2,'Sv'}, ...
          {'AI',1/1e15,'$\mathcal{A}_I$',[-1.25 1.25],0.05,'PW'}, ...
};

% Diathermal components nice plot:
fields = { ...
          {'Fall',-1./repmat(dy,[1 TL+1])/1e12,'Surface Forcing',[-50 50],5,'TW/$^\circ$latitude'}, ...
          {'Mall',-1./repmat(dy,[1 TL+1])/1e12,'Mixing',[-50 50],5,'TW/$^\circ$latitude'}, ...
% $$$           {'N',1./repmat(dy,[1 TL+1])/1e12,'',[-50 50],5,'TW/$^\circ$latitude'}, ...
% $$$           {'Jdia',1./repmat(dy,[1 TL+1])/1e12,'Total',[-50 50],5,'TW/$^\circ$latitude'}, ...
};

rego = [3 1 2];
names = {'Global','Atlantic','Indo-Pacific'};
clab = [1 1 1 1 1 1];

cpts = cell(1,length(fields));
for i=1:length(fields)
    cpts{i} = [-1e10 fields{i}{4}(1):fields{i}{5}:fields{i}{4}(2) 1e10];
end
npts = length(cpts{1});
clab = [1 1 1 1 1 1];

cmap = redblue(npts-3);
% $$$ cmap = parula(npts-3);
% $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

AIsp = 0.25;
% $$$ AIsp = 1;
latfilt = 1;

% $$$ doZAremap = 0; % remap to depth space

%Fluxes only:
figure;
set(gcf,'Position',[2125          11        1680         960]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

% 2x3:
poss = [0.11     0.5949    0.25      0.3301; ...
        0.375    0.5949    0.1875    0.3301; ...
        0.58     0.5949    0.1719    0.3301; ...
        0.11     0.2300    0.25      0.3301; ...
        0.375    0.2300    0.1875    0.3301; ...
        0.58     0.2300    0.1719    0.3301];
        
letlabs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for i=1:length(fields)
    for r = 1:length(regLets)
        subplot(length(fields),3,3*(i-1)+r);
        reg = 1;%rego(r);
% $$$         if (doZAremap)
% $$$             Yg = repmat(yvec,[1 TL+1]);
% $$$             Tg = zZTxa;
% $$$         else
        [X,Y] = ndgrid(yt,Te);
% $$$         end

        eval(['VAR = ZA_' regLets{reg} '.' fields{i}{1} '.*fields{i}{2};']);
        VAR(VAR==0) = NaN;
        eval(['VAR(ZA_' regLets{reg} '.NaNst==1) = NaN;']);
        VAR = filter_field(VAR',latfilt,'-t')';
        contourf(X,Y,VAR,cpts{i},'linestyle','none');
        hold on;
        col = [0 0 0];

% $$$         if (~doZAremap)
% $$$         plot(yvec,filter_field(meanSST,latfilt,'-t'),':','color',col);
        eval(['plot(yt,filter_field(ZA_' regLets{reg} '.minSST,latfilt,''-t''),'':'',''color'',col);']);
        eval(['plot(yt,filter_field(ZA_' regLets{reg} '.maxTt,latfilt,''-t''),'':k'');']);

% $$$     else
% $$$         [tt,zz] = ndgrid(yvec,-z);
% $$$         [c,h] = contour(tt,zz,zTxa,[0:4:34],'-k');
% $$$         clabel(c,h);
% $$$     end

        if (i>=1)
            eval(['VAR = ZA_' regLets{reg} '.AI/1e15;']);        
            eval(['VAR(ZA_' regLets{reg} '.NaNst==1) = NaN;']);
            VAR = filter_field(VAR',latfilt,'-t')';
            [c,h] = contour(X,Y,VAR,[-50:AIsp:-AIsp],'--k');
            if (clab(i))
                clabel(c,h);
            end
            [c,h] =  contour(X,Y,VAR,[AIsp:AIsp:50],'-k');
            if (clab(i))
                clabel(c,h);
            end
            
% $$$             if (i==2 & r == 3)
% $$$             % 0-contour on Indo-Pacific mixing:
% $$$             VAR(X>39) = 0.25;
% $$$             [c,h] = contour(X,Y,VAR,[0 0],'--','color',[0 0.5 0]);
% $$$             clabel(c,h);
% $$$             end
        end
        
        ylim([-3 34]);
        plot([-34 -34],[-3 34],'--k');
        if (r > 1)
            plot([66 66],[-3 20],'--k');
        end
% $$$         plot([-45 -45],[-3 34],'--k');
        caxis(fields{i}{4});
        box on; 
        grid on;
        letno = 3*(i-1)+r;
% $$$         if (strcmp(fields{i}{1},'Mall'))
% $$$             letno = letno+3;
% $$$         end
        if (r == 1)
            xlim([-80 80]);
            text(-79,32.15,[letlabs{letno} ' ' fields{i}{3}]);
        elseif (r == 2)
            xlim([-40 80]);
            text(-33,32.15,[letlabs{letno} ' ' fields{i}{3}]);
        else
            xlim([-40 70]);
            text(-33,32.15,[letlabs{letno} ' ' fields{i}{3}]);
        end
        set(gca,'xtick',[-90:30:90]);
% $$$         if (i > 1)
            xlabel('Latitude ($^\circ$N)');
% $$$         else
% $$$             set(gca,'xticklabel',[]);
% $$$         end
        if (r == 1)
            ylabel('Temperature $\Theta$ ($^\circ$C)');
        else
            set(gca,'yticklabel',[]);
        end
        if (r == 3)
            cb = colorbar;
            ylabel(cb,fields{i}{6});
        end
        if (i == 1)
            title(names{r});
        end
        set(gca,'Position',poss(3*(i-1)+r,:));
    end
end
colormap(cmap);

%% Plot global mixing components:
dy = diff(yu);
dy = [dy(1); dy];

% Mixing components plot:
fields = { ...
          {'M',-1./repmat(dy,[1 TL+1])/1e12,'Total Vertical Mixing',[-40 0],2,'TW/$^\circ$latitude'}, ...
          {'Mkppiw',-1./repmat(dy,[1 TL+1])/1e12,'Background Mixing',[-20 0],1,'TW/$^\circ$latitude'}, ...
          {'Mkppish',-1./repmat(dy,[1 TL+1])/1e12,'Shear Instability',[-20 0],1,'TW/$^\circ$latitude'}, ...
          {'I',-1./repmat(dy,[1 TL+1])/1e12,'Numerical Mixing',[-40 0],2,'TW/$^\circ$latitude'}, ...
          {'Mkppbl',-1./repmat(dy,[1 TL+1])/1e12,'KPP Boundary Layer',[-20 0],1,'TW/$^\circ$latitude'}, ...
          {'Mwave',-1./repmat(dy,[1 TL+1])/1e12,'Internal Tide',[-20 0],1,'TW/$^\circ$latitude'}, ...
};
clab = [1 0 0 0 0 0];

cpts = cell(1,length(fields));
for i=1:length(fields)
    cpts{i} = [-1e10 fields{i}{4}(1):fields{i}{5}:fields{i}{4}(2) 1e10];
end
npts = length(cpts{1});
clab = [1 0 0 0 0 0];

cmap = parula(npts-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

AIsp = 0.25;
latfilt = 1;

%Fluxes only:
figure;
set(gcf,'Position',[2125          11        1680         960]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

% 2x3:
poss = [0.05     0.5949    0.24    0.3301; ...
        0.36    0.5949    0.24    0.3301; ...
        0.615     0.5949    0.24    0.3301; ...
        0.05     0.2300    0.24    0.3301; ...
        0.36    0.2300    0.24    0.3301; ...
        0.615     0.2300    0.24    0.3301];
        
letlabs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
rego = 'G';
for i=1:length(fields)
    subplot(2,3,i);
    [X,Y] = ndgrid(yt,Te);

    eval(['VAR = ZA_' rego '.' fields{i}{1} '.*fields{i}{2};']);
    VAR(VAR==0) = NaN;
    eval(['VAR(ZA_' rego '.NaNst==1) = NaN;']);
    VAR = filter_field(VAR',latfilt,'-t')';
    contourf(X,Y,VAR,cpts{i},'linestyle','none');
    hold on;
    col = [0 0 0];

    eval(['plot(yt,filter_field(ZA_' rego '.minSST,latfilt,''-t''),'':'',''color'',col);']);
    eval(['plot(yt,filter_field(ZA_' rego '.maxTt,latfilt,''-t''),'':k'');']);

    eval(['VAR = ZA_' rego '.AI/1e15;']);        
    eval(['VAR(ZA_' rego '.NaNst==1) = NaN;']);
    VAR = filter_field(VAR',latfilt,'-t')';
    [c,h] = contour(X,Y,VAR,[-50:AIsp:-AIsp],'--k');
    if (clab(i))
        clabel(c,h);
    end
    [c,h] =  contour(X,Y,VAR,[AIsp:AIsp:50],'-k');
    if (clab(i))
        clabel(c,h);
    end
        
    ylim([-3 34]);
    caxis(fields{i}{4});
    box on; 
    grid on;
    xlim([-80 80]);
    text(-79,32.15,[letlabs{i} ' ' fields{i}{3}]);
    set(gca,'xtick',[-90:30:90]);
    if (i >= 4)
        xlabel('Latitude ($^\circ$N)');
    else
        set(gca,'xticklabel',[]);
    end
    if (i == 1 | i == 4)
        ylabel('Temperature $\Theta$ ($^\circ$C)');
    else
        set(gca,'yticklabel',[]);
    end
    if (i == 1 | i == 3 | i == 4 | i == 6)
        cb = colorbar;
        ylabel(cb,fields{i}{6});
        set(cb,'ytick',[-40:5:0]);
    end
    set(gca,'Position',poss(i,:));
end
colormap(cmap);

%% Plot MHT:
figure;
set(gcf,'Position',[1           1        1920         962]);
colors = {'-k','-b','-r','-m','-c'};
lfilt = 1;
% $$$ MHTE = rho0*Cp*ZA_G.psiTu.*ZA_G.maxTu;
% $$$ MHTI = ZA_G.MHT - avg([0; MHTE]);
ZA_A.MHT(ZA_A.MHT==0) = NaN;
ZA_P.MHT(ZA_P.MHT==0) = NaN;
% $$$ [tmp ind] = min(abs(Te-15));
% $$$ MHTA15 = ZA_A.AI(:,ind);
% $$$ MHTA15(MHTA15==0) = NaN;
plot(yt,filter_field(ZA_A.MHT/1e15,lfilt,'-t'),':r','linewidth',3);
hold on; 
plot(yt,filter_field(ZA_P.MHT/1e15,lfilt,'-t'),':b','linewidth',3);
plot(yt,filter_field(ZA_G.MHT/1e15,lfilt,'-t'),':k','linewidth',3);
% $$$ plot(yt,filter_field(MHTI/1e15,lfilt,'-t'),'--k','linewidth',2);
% $$$ plot(yu,filter_field(MHTE/1e15,lfilt,'-t'),':k','linewidth',2);
% $$$ plot(yu,filter_field(MHTA15/1e15,lfilt,'-t'),'--r','linewidth',2);
legend({'Atlantic','Indo-Pacific','Global'});%'Global Internal','Global External'});%,'Atlantic Internal $<15^\circ$C'});
xlabel('Latitude ($^\circ$N)');
ylabel('PW');
xlim([-80 80]);
ylim([-1.3 1.8]);
set(gca,'xtick',[-90:30:90]);
grid on;
set(gca,'Position',[0.2196    0.3809    0.4298    0.3662]);

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

