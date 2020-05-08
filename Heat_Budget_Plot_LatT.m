% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

% Load Base Variables:
model = 'MOM025_kb3seg';
outputs = [101110];
%outputs = [121];
% $$$ outputs = [111120];

% $$$ model = 'MOM025_nipoall';
% $$$ outputs = [10:19];

% $$$ model = 'MOM025_RCP45';
% $$$ outputs = [0:39];

% $$$ model = 'MOM025_SOUP15';
% $$$ outputs = [10:19];

% model = 'MOM01';
% outputs = [4567];

% $$$ model = 'ACCESS-OM2_1deg_jra55_iaf';
% $$$ outputs = [39];

model = 'ACCESS-OM2_1deg_ryf_4dt';
outputs = [51];

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
end
% Latitude difference vector for plotting per-degree:
dy = [yu(2)-yu(1); diff(yu)]; % (First-element is done by hand - but dy is equal to second).

regions = {'Global'};
regLets = {'G'};
% $$$ regions = {'Atlantic2BAS','IndoPacific2BAS','Global'};%,'SO_Atlantic','SO_IndoPacific'};
% $$$ regLets = {'A','P','G'};%,'SA','SP'};

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
% $$$ if (~isfield(ZA_P,'KPPNL'))
% $$$     ZA_P.KPPNL = 0*ZA_P.M;
% $$$ end
% $$$ if (~isfield(ZA_A,'KPPNL'))
% $$$     ZA_A.KPPNL = 0*ZA_A.M;
% $$$ end

% Global:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_G.AI-ZA_G.AHDR),[],1); % convergence of that transport
ZA_G.Jdia = -ZA_G.N-dAI_mR_dphi; % total diathermal transport
ZA_G.I = -(ZA_G.Jdia+ZA_G.M+ZA_G.KPPNL+ZA_G.F+ZA_G.PI+ZA_G.RED+ZA_G.K33+ZA_G.MDS+ZA_G.SIG); % numerical mixing (both advective and submesoscale)
if (isfield(ZA_G,'NUM_SUBlf'))
    ZA_G.NUM = ZA_G.NUM_SUBlf;
end
% XXX: Inclusion of Redi here is confusing but correct I think - it
% should be included in AI (as it is), but to calculate numerical
% mixing you subtract it off and use the 3D convergence instead (which
% includes both diathermal and meridional convergences). I've never
% plotted the diathermal heat flux due to redi diffusion (I guess this
% is the difference between the meridional flux convergence and the 3D
% convergence - or is it???)

% $$$ % Indo-Pacific:
% $$$ dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_P.AI-ZA_P.AHDR),[],1); % convergence of that transport
% $$$ ZA_P.Jdia = -ZA_P.N-dAI_mR_dphi; % total diathermal transport
% $$$ ZA_P.I = -(ZA_P.Jdia+ZA_P.M+ZA_P.KPPNL+ZA_P.F+ZA_P.PI+ZA_P.RED+ZA_P.K33+ZA_P.MDS+ZA_P.SIG); % numerical mixing (both advective and submesoscale)
% $$$ if (isfield(ZA_P,'NUM_SUBlf'))
% $$$     ZA_P.NUM = ZA_P.NUM_SUBlf;
% $$$ end
% $$$ 
% $$$ % Atlantic:
% $$$ dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_A.AI-ZA_A.AHDR),[],1); % convergence of that transport
% $$$ ZA_A.Jdia = -ZA_A.N-dAI_mR_dphi; % total diathermal transport
% $$$ ind = find(~isnan(ZA_P.AI(:,end)),1,'last'); % 66N indicy (u-points)
% $$$ ZA_A.Jdia(ind,:) = ZA_A.Jdia(ind,:) + ZA_P.AI(ind-1,:); % Correct for missing bit because main Atlantic overlaps with Bering Strait cutoff
% $$$ ZA_A.I = -(ZA_A.Jdia+ZA_A.M+ZA_A.KPPNL+ZA_A.F+ZA_A.PI+ZA_A.RED+ZA_A.K33+ZA_A.MDS+ZA_A.SIG); % numerical mixing (both advective and submesoscale)
% $$$ if (isfield(ZA_A,'NUM_SUBlf'))
% $$$     ZA_A.NUM = ZA_A.NUM_SUBlf;
% $$$ end
% $$$ 
% $$$ if (exist('ZA_SA'))
% $$$ % SO Atlantic:
% $$$ dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_SA.AI-ZA_SA.AHDR),[],1); % convergence of that transport
% $$$ ZA_SA.Jdia = -ZA_SA.N-dAI_mR_dphi; % total diathermal transport
% $$$ ZA_SA.I = -(ZA_SA.Jdia+ZA_SA.M+ZA_SA.KPPNL+ZA_SA.F+ZA_SA.PI+ZA_SA.RED+ZA_SA.K33+ZA_SA.MDS+ZA_SA.SIG); % numerical mixing (both advective and submesoscale)
% $$$ end
% $$$ 
% $$$ if (exist('ZA_SP'))
% $$$ % SO IndoPacific:
% $$$ dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_SP.AI-ZA_SP.AHDR),[],1); % convergence of that transport
% $$$ ZA_SP.Jdia = -ZA_SP.N-dAI_mR_dphi; % total diathermal transport
% $$$ ZA_SP.I = -(ZA_SP.Jdia+ZA_SP.M+ZA_SP.KPPNL+ZA_SP.F+ZA_SP.PI+ZA_SP.RED+ZA_SP.K33+ZA_SP.MDS+ZA_SP.SIG); % numerical mixing (both advective and submesoscale)
% $$$ end

groups = {'G'};%,'P','A'}%,'SA','SP'};
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
          groups{gi} '.M-ZA_' groups{gi} '.I-ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.N' ...
          '-ZA_' groups{gi} '.RED-ZA_' groups{gi} '.K33-ZA_' groups{gi} '.MDS-ZA_' groups{gi} '.SIG;']);

    % Fix Atlantic Bering Strait part:
    if (strcmp(groups{gi},'A'))
        ind = find(~isnan(ZA_P.AI(:,end)),1,'last'); % 66N indicy (u-points)
        eval(['ZA_' groups{gi} '.AIF(ind,:) = ZA_' groups{gi} '.AIF(ind,:)-ZA_P.AI(ind-1,:);']);
    end
    
    eval(['ZA_' groups{gi} '.Fall = ZA_' groups{gi} '.F + ZA_' groups{gi} '.PI+ZA_' groups{gi} '.AIF;']);
    eval(['ZA_' groups{gi} '.Mall = ZA_' groups{gi} '.M + ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.I;']);
end

% $$$ % Check Bering Strait/net psi:
% $$$ 
% $$$ figure;
% $$$ subplot(2,2,1);
% $$$ pcolPlot(avg(X,1),avg(Y,1),-rho0*Cp*diff(ZA_P.PSI,[],1));
% $$$ hold on;
% $$$ plot(yt,ZA_P.maxTt,'--k');
% $$$ caxis([-1 1]*1e11);
% $$$ subplot(2,2,2);
% $$$ pcolPlot(avg(X,2),avg(Y,2),diff(ZA_P.F,[],2)/dT);
% $$$ hold on;
% $$$ plot(yt,ZA_P.maxTt,'--k');
% $$$ caxis([-1 1]*1e11);
% $$$ subplot(2,2,3);
% $$$ pcolPlot(avg(X,2),avg(Y,2),diff(ZA_P.PI,[],2)/dT);
% $$$ hold on;
% $$$ plot(yt,ZA_P.maxTt,'--k');
% $$$ caxis([-1 1]*1e11);
% $$$ subplot(2,2,4);
% $$$ pcolPlot(avg(X,2),avg(Y,2),diff(ZA_P.M,[],2)/dT);
% $$$ hold on;
% $$$ plot(yt,ZA_P.maxTt,'--k');
% $$$ caxis([-1 1]*1e11);
% $$$ colormap(redblue);
% $$$ 
% $$$ figure;
% $$$ ind = find(~isnan(ZA_P.AI(:,end)),1,'last');
% $$$ plot(ZA_P.AI(ind,:)/1e15,Te);
% $$$ %caxis([
% $$$ 
% $$$ %%% Difference two runs:
% $$$ % $$$ ZA_Gc = ZA_G;ZA_Pc = ZA_P;ZA_Ac = ZA_A;
% $$$ 
% $$$ fields = fieldnames(ZA_G);
% $$$ for f = 1:length(fields)
% $$$     if (~(strcmp(fields{f},'NaNst') | strcmp(fields{f},'NaNsu')))
% $$$         eval(['sz = size(ZA_G.' fields{f} ');']);
% $$$         if (sz(2)>1 | strcmp(fields{f}(1:3),'MHT'))
% $$$             eval(['ZA_G.' fields{f} ' = ZA_G.' fields{f} ' - ZA_Gc.' fields{f} ';']);
% $$$             eval(['ZA_A.' fields{f} ' = ZA_A.' fields{f} ' - ZA_Ac.' fields{f} ';']);
% $$$             eval(['ZA_P.' fields{f} ' = ZA_P.' fields{f} ' - ZA_Pc.' fields{f} ';']);
% $$$         end
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ %%%% Summary schematics:
% $$$ ltminW = -45;
% $$$ ltmin = -34;
% $$$ ltmax = 90;
% $$$ [tmp ltminI] = min(abs(yu-ltmin)); % u-point index
% $$$ [tmp ltminWI] = min(abs(yu-ltminW)); % 45S
% $$$ [tmp ltmaxI] = min(abs(yu-ltmax)); % u-point index
% $$$ [tmp lt50I] = min(abs(yu-50)); % u-point index
% $$$ indBS = find(~isnan(ZA_P.AI(:,end)),1,'last');
% $$$ [X,Y] = ndgrid(yt,Te);
% $$$ 
% $$$ % Below a given temperature:
% $$$ % Atl:
% $$$ FAb = nansum(ZA_A.Fall(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ MAb = nansum(ZA_A.Mall(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ NAb = nansum(ZA_A.N(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ AA34Sb = ZA_A.AI(ltminI,:)/1e15;
% $$$ AA50Nb = ZA_A.AI(lt50I,:)/1e15;
% $$$ % Pac:
% $$$ FPb = nansum(ZA_P.Fall(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ MPb = nansum(ZA_P.Mall(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ NPb = nansum(ZA_P.N(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ AP34Sb = ZA_P.AI(ltminI,:)/1e15;
% $$$ % Glo and BS:
% $$$ A34Sb = ZA_G.AI(ltminI,:)/1e15;
% $$$ ABSPb = ZA_P.AI(indBS,:)/1e15;
% $$$ 
% $$$ % SO:
% $$$ FSOb = nansum(ZA_G.Fall(1:ltminI,:),1)/1e15;
% $$$ MSOb = nansum(ZA_G.Mall(1:ltminI,:),1)/1e15;
% $$$ NSOb = nansum(ZA_G.N(1:ltminI,:),1)/1e15;
% $$$ 
% $$$ % SO and warm route:
% $$$ FWRb = nansum(ZA_G.Fall((ltminWI+1):ltminI,:),1)/1e15;
% $$$ MWRb = nansum(ZA_G.Mall((ltminWI+1):ltminI,:),1)/1e15;
% $$$ NWRb = nansum(ZA_G.N((ltminWI+1):ltminI,:),1)/1e15;
% $$$ 
% $$$ if (exist('ZA_SA'))
% $$$ % SO Atlantic and SO Pacific surface forcing:
% $$$ FSAb = nansum(ZA_SA.Fall(1:ltminI,:),1)/1e15;
% $$$ FSPb = nansum(ZA_SP.Fall(1:ltminI,:),1)/1e15;
% $$$ FWRAb = nansum(ZA_SA.Fall((ltminWI+1):ltminI,:),1)/1e15;
% $$$ FWRPb = nansum(ZA_SP.Fall((ltminWI+1):ltminI,:),1)/1e15;
% $$$ else
% $$$ FSAb = zeros(size(Te));
% $$$ FSPb = FSAb;FWRAb = FSAb;FWRPb=FSAb;    
% $$$ end
% $$$ 
% $$$ A45Sb = ZA_G.AI(ltminWI,:)/1e15;
% $$$ 
% $$$ % Residuals by basin:
% $$$ RESP = FPb+MPb+AP34Sb-NPb-ABSPb;
% $$$ RESA = FAb+MAb+AA34Sb-NAb+ABSPb;
% $$$ RESS = FSOb+MSOb-A34Sb-NSOb; % All checks out.
% $$$ 
% $$$ % Summary numbers at temps:
% $$$ %sT = [-3 5 10 15 20 34];
% $$$ sT = [-3 15 20 34];
% $$$ for Ti=(length(sT)-1):-1:1
% $$$     [tmp indL] = min(abs(Te-sT(Ti)));
% $$$     [tmp indU] = min(abs(Te-sT(Ti+1)));
% $$$     
% $$$     disp(' ')
% $$$     disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
% $$$     disp(['-------------------------------------------------------------'])
% $$$     disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f (Atl = %5.2f, Pac = %5.2f), WR = %5.2f (Atl = %5.2f, Pac = %5.2f)',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL),FSAb(indU)-FSAb(indL),FSPb(indU)-FSPb(indL),FWRb(indU)-FWRb(indL),FWRAb(indU)-FWRAb(indL),FWRPb(indU)-FWRPb(indL)))
% $$$     disp(sprintf('Mixing      (flux at bottom) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',MPb(indL),MAb(indL),MSOb(indL),MWRb(indL)))
% $$$     disp(sprintf('Transport 34S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
% $$$     disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
% $$$     disp(sprintf('Transport 50N   (into layer) Atlantic = %5.2f',AA50Nb(indU)-AA50Nb(indL)))
% $$$     disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
% $$$     disp(sprintf('Transport BS    (into layer)            %5.2f.',ABSPb(indU)-ABSPb(indL)))
% $$$     disp(sprintf('Residuals                    Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',RESP(indU)-RESP(indL),RESA(indU)-RESA(indL),RESS(indU)-RESS(indL)))
% $$$ end
% $$$ 
% $$$ % Global:
% $$$ indL = 1;indU=TL+1;
% $$$ for i=1:1
% $$$     disp(' ')
% $$$     disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
% $$$     disp(['-------------------------------------------------------------'])
% $$$     disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL),FWRb(indU)-FWRb(indL)))
% $$$     disp(sprintf('Mixing      (flux at bottom) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',MPb(indL),MAb(indL),MSOb(indL),MWRb(indL)))
% $$$     disp(sprintf('Transport 32S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
% $$$     disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
% $$$     disp(sprintf('Transport 50N   (into layer) Atlantic = %5.2f',AA50Nb(indU)-AA50Nb(indL)))
% $$$     disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f, WR = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL),NWRb(indU)-NWRb(indL)))
% $$$     disp(sprintf('Transport BS    (into layer)          = %5.2f',ABSPb(indU)-ABSPb(indL)))
% $$$     disp(sprintf('Residuals                    Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',RESP(indU)-RESP(indL),RESA(indU)-RESA(indL),RESS(indU)-RESS(indL)))
% $$$ end

% $$$ % Global 0C-ref:
% $$$ FAb = nansum(ZA_A.F(ltminI+1:ltmaxI,:)+ZA_A.P(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ NAb = nansum(ZA_A.N(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ AA34Sb = ZA_A.AHD(ltminI,:)/1e15;
% $$$ AA50Nb = ZA_A.AHD(lt50I,:)/1e15;
% $$$ 
% $$$ % Pac:
% $$$ FPb = nansum(ZA_P.F(ltminI+1:ltmaxI,:)+ZA_P.P(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ NPb = nansum(ZA_P.N(ltminI+1:ltmaxI,:),1)/1e15;
% $$$ AP34Sb = ZA_P.AHD(ltminI,:)/1e15;
% $$$ % Glo and BS:
% $$$ A34Sb = ZA_G.AHD(ltminI,:)/1e15;
% $$$ ABSPb = ZA_P.AHD(indBS,:)/1e15;
% $$$ A45Sb = ZA_G.AHD(ltminWI,:)/1e15;
% $$$ 
% $$$ % SO:
% $$$ FSOb = nansum(ZA_G.F(1:ltminI,:)+ZA_G.P(1:ltminI,:),1)/1e15;
% $$$ NSOb = nansum(ZA_G.N(1:ltminI,:),1)/1e15;
% $$$ 
% $$$ indL = 1;indU=TL+1;
% $$$ for i=1:1
% $$$     disp(' ')
% $$$     disp(['Temperatures ' num2str(Te(indL)) 'C - ' num2str(Te(indU)) 'C:'])
% $$$     disp(['-------------------------------------------------------------'])
% $$$     disp(sprintf('Surface Forcing (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',FPb(indU)-FPb(indL),FAb(indU)-FAb(indL),FSOb(indU)-FSOb(indL)))
% $$$     disp(sprintf('Transport 32S   (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',AP34Sb(indU)-AP34Sb(indL),AA34Sb(indU)-AA34Sb(indL),A34Sb(indU)-A34Sb(indL)))
% $$$     disp(sprintf('Transport 45S   (into layer) SO = %5.2f',A45Sb(indU)-A45Sb(indL)))
% $$$     disp(sprintf('Transport 50N   (into layer) Atlantic = %5.2f',AA50Nb(indU)-AA50Nb(indL)))
% $$$     disp(sprintf('Tendency        (into layer) Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',NPb(indU)-NPb(indL),NAb(indU)-NAb(indL),NSOb(indU)-NSOb(indL)))
% $$$     disp(sprintf('Transport BS    (into layer)          = %5.2f',ABSPb(indU)-ABSPb(indL)))
% $$$ % $$$     disp(sprintf('Residuals                    Indo-Pac = %5.2f, Atlantic = %5.2f, SO = %5.2f',RESP(indU)-RESP(indL),RESA(indU)-RESA(indL),RESS(indU)-RESS(indL)))
% $$$ end

%% Plot latitude - temperature plane for different basins:
dy = diff(yu);
dy = [dy(1); dy];

% PsiAI nice plot:
fields = { ...
          {'PSI',1/1e6,'$\Psi$',[-30 30],2,'Sv'}, ...
          {'AI',1/1e15,'$\mathcal{A}_I$',[-1.25 1.25],0.05,'PW'}, ...
% Perturbations:
% $$$           {'PSI',1/1e6,'$\Psi$',[-5 5],0.25,'Sv'}, ...
% $$$ % $$$           {'AI',1/1e15,'$\mathcal{A}_I$',[-0.25 0.25],0.002,'PW'}, ...
% $$$           {'AI',1/1e15,'$\mathcal{A}_I$',[-0.15 0.15],0.001,'PW'}, ...
};

% $$$ % $$$ % Diathermal components nice plot:
fields = { ...
% $$$           {'Fall',-1./repmat(dy,[1 TL+1])/1e12,'Surface Forcing',[-50 50],5,'TW/$^\circ$latitude'}, ...
% $$$           {'Mall',-1./repmat(dy,[1 TL+1])/1e12,'Mixing',[-50 50],5,'TW/$^\circ$latitude'}, ...
          {'N',1./repmat(dy,[1 TL+1])/1e12,'',[-50 50],5,'TW/$^\circ$latitude'}, ...
          {'Jdia',1./repmat(dy,[1 TL+1])/1e12,'Total',[-50 50],5,'TW/$^\circ$latitude'}, ...
% Perturbations:
% $$$           {'Fall',-1./repmat(dy,[1 TL+1])/1e12,'Surface Forcing',[-20 20],0.5,'TW/$^\circ$latitude'}, ...
% $$$           {'Mall',-1./repmat(dy,[1 TL+1])/1e12,'Mixing',[-20 20],0.5,'TW/$^\circ$latitude'}, ...
% $$$           {'Jdia',1./repmat(dy,[1 TL+1])/1e12,'Total',[-20 20],0.5,'TW/$^\circ$latitude'}, ...
% $$$           {'N',1./repmat(dy,[1 TL+1])/1e12,'Tendency',[-20 20],0.5,'TW/$^\circ$latitude'}, ...
};

% NumMix:
fields = { ...
          {'I',-1./repmat(dy,[1 TL+1])/1e12,'Numerical Mixing',[-20 20],1,'TW/$^\circ$latitude'}, ...
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

% MHT contributions:
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

%% Plot lat-T diffusivity:
load([base model sprintf('_output%03d_diff_cbt_T.mat',121)]);%outputs(1))]);

diff_cbt_G = nanmonmean(diff_cbt_T_G,3,ndays);
diff_cbt_P = nanmonmean(diff_cbt_T_P,3,ndays);
diff_cbt_A = nanmonmean(diff_cbt_T_A,3,ndays);
diff_cbt_G(diff_cbt_G==0) = NaN;
diff_cbt_P(diff_cbt_P==0) = NaN;
diff_cbt_A(diff_cbt_A==0) = NaN;

[X,Y] = ndgrid(yt,T);

load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
SST_G(SST_G==0) = NaN;
minSST_G = squeeze(min(monmean(SST_G,3,ndays),[],1)');
if (max(minSST_G)>100); minSST_G = minSST_G-273.15; end
maxSST_G = NaN*zeros(size(minSST_G));
maxSSTi_G = NaN*maxSST_G;
for i=1:yL
    fnds = find(~isnan(diff_cbt_G(i,:)),1,'last');
    if length(fnds)>0
        maxSSTi_G(i) = fnds;
        maxSST_G(i) = T(fnds);
    end
end
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
SST_P(SST_P==0) = NaN;
minSST_P = squeeze(min(monmean(SST_P,3,ndays),[],1)');
if (max(minSST_P)>100); minSST_P = minSST_P-273.15; end
maxSST_P = NaN*zeros(size(minSST_P));
maxSSTi_P = NaN*maxSST_P;
for i=1:yL
    fnds = find(~isnan(diff_cbt_P(i,:)),1,'last');
    if length(fnds)>0
        maxSSTi_P(i) = fnds;
        maxSST_P(i) = T(fnds);
    end
end
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
SST_A(SST_A==0) = NaN;
minSST_A = squeeze(min(monmean(SST_A,3,ndays),[],1)');
if (max(minSST_A)>100); minSST_A = minSST_A-273.15; end
maxSST_A = NaN*zeros(size(minSST_A));
maxSSTi_A = NaN*maxSST_A;
for i=1:yL
    fnds = find(~isnan(diff_cbt_A(i,:)),1,'last');
    if length(fnds)>0
        maxSSTi_A(i) = fnds;
        maxSST_A(i) = T(fnds);
    end
end

clvls = [0 10.^[-6:0.25:1] 1e6]
subplot(1,3,1);
contourf(X,Y,diff_cbt_G,clvls,'linestyle','none');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'$\kappa$ (m$^2$s$^{-1}$)','Interpreter','latex');
set(gca,'colorscale','log');
colormap(jet);
caxis([1e-5-1e-7 1e-2]);
hold on;
plot(yt,maxSST_G,':k');
plot(yt,minSST_G,':k');
title('Global');
subplot(1,3,2);
contourf(X,Y,diff_cbt_A,clvls,'linestyle','none');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'$\kappa$ (m$^2$s$^{-1}$)','Interpreter','latex');
set(gca,'colorscale','log');
colormap(jet);
caxis([1e-5-1e-7 1e-2]);
hold on;
plot(yt,maxSST_A,':k');
plot(yt,minSST_A,':k');
title('Atlantic');
subplot(1,3,3);
contourf(X,Y,diff_cbt_P,clvls,'linestyle','none');
xlabel('Latitude ($^\circ$N)');
ylabel('Temperature ($^\circ$C)');
cb = colorbar;
ylabel(cb,'$\kappa$ (m$^2$s$^{-1}$)','Interpreter','latex');
set(gca,'colorscale','log');
colormap(jet);
caxis([1e-5-1e-7 1e-2]);
hold on;
plot(yt,maxSST_P,':k');
plot(yt,minSST_P,':k');
title('Indo-Pacific');

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




