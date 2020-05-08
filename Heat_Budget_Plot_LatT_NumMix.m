% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';


RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[4567]}, ...
% $$$ % MOM025-SIS:
    {'MOM025_kb3seg',[101120],'(a) MOM025 Control'}, ...
% $$$     {'MOM025',[15:19]}, ...
% $$$     {'MOM025_kb1em5',[95:99]}, ...
% $$$     {'MOM025_kb1em6',[30]}, ...
% $$$ % ACCESS-OM2 Gadi runs:
         {'ACCESS-OM2_1deg_jra55_ryf',[31],'(b) ACCESS-OM2-1-KDS50'}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135],'(c) ACCESS-OM2-1-KDS135'}, ...
         {'ACCESS-OM2_025deg_jra55_ryf',[7680],'(d) ACCESS-OM2-025-RG'}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf8485_gmredi',[73]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680],'(e) ACCESS-OM2-025'}, ...
% $$$     {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[80]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[300]}, ...
         {'ACCESS-OM2_01deg_jra55_ryf',[636639],'(f) ACCESS-OM2-01'}, ...
       };

for rr = 1:length(RUNS)
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};
    label = RUNS{rr}{3};

load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
end
% Latitude difference vector for plotting per-degree:
dy = [yu(2)-yu(1); diff(yu)]; % (First-element is done by hand - but dy is equal to second).

regions = {'Global'};
regLets = {'G'};

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

% Global:
dAI_mR_dphi = diff(cat(1,zeros(1,TL+1),ZA_G.AI-ZA_G.AHDR),[],1); % convergence of that transport
ZA_G.Jdia = -ZA_G.N-dAI_mR_dphi; % total diathermal transport
ZA_G.I = -(ZA_G.Jdia+ZA_G.M+ZA_G.KPPNL+ZA_G.F+ZA_G.PI+ZA_G.RED+ZA_G.K33+ZA_G.MDS+ZA_G.SIG); % numerical mixing (both advective and submesoscale)
if (isfield(ZA_G,'NUM_SUBlf'))
    ZA_G.NUM = ZA_G.NUM_SUBlf;
end

% $$$ groups = {'G'};%,'P','A'}%,'SA','SP'};
% $$$ for gi=1:length(groups)
% $$$     eval(['ZA_' groups{gi} '.AIF = 0*ZA_G.F;']); % Heat lost max SST line from AI
% $$$ 
% $$$     % Fix above max SST to calculate budgets between theta or max SST
% $$$     % and -2C:
% $$$     for yi=1:yL
% $$$         %PI:
% $$$         eval(['ZA_' groups{gi} '.PI(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
% $$$               '= ZA_' groups{gi} '.PI(yi,ZA_' groups{gi} '.maxTit(yi));']);
% $$$         %tendency:
% $$$         eval(['ZA_' groups{gi} '.N(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
% $$$               '= ZA_' groups{gi} '.N(yi,ZA_' groups{gi} '.maxTit(yi));']);
% $$$         %Num-mix:
% $$$         eval(['ZA_' groups{gi} '.I(yi,ZA_' groups{gi} '.maxTit(yi):end) ' ...
% $$$               '= ZA_' groups{gi} '.I(yi,ZA_' groups{gi} '.maxTit(yi));']);
% $$$         %AI:
% $$$         eval(['ZA_' groups{gi} '.AI(yi,ZA_' groups{gi} '.maxTiu(yi):end) ' ...
% $$$               '= ZA_' groups{gi} '.AI(yi,ZA_' groups{gi} '.maxTiu(yi));']);
% $$$     end
% $$$     
% $$$     % Correct surface forcing by residual:
% $$$     eval(['ZA_' groups{gi} '.AIF = diff(cat(1,zeros(1,TL+1),ZA_' groups{gi} '.AI-ZA_' ...
% $$$           groups{gi} '.AHDR),[],1)-ZA_' groups{gi} '.PI-ZA_' groups{gi} '.F-ZA_' ...
% $$$           groups{gi} '.M-ZA_' groups{gi} '.I-ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.N;']);
% $$$ 
% $$$     eval(['ZA_' groups{gi} '.Fall = ZA_' groups{gi} '.F + ZA_' groups{gi} '.PI+ZA_' groups{gi} '.AIF;']);
% $$$     eval(['ZA_' groups{gi} '.Mall = ZA_' groups{gi} '.M + ZA_' groups{gi} '.KPPNL+ZA_' groups{gi} '.I;']);
% $$$ end

%% Plot latitude - temperature plane for different basins:
dy = diff(yu);
dy = [dy(1); dy];

% NumMix:
fields = { ...
% $$$           {'I',-1./repmat(dy,[1 TL+1])/1e12,'Numerical Mixing',[-10 10],0.25,'TW/$^\circ$latitude'}, ...
          {'NUM',1./repmat(dy,[1 TL+1])/1e12,'Numerical Mixing',[-30 0],1,'TW/$^\circ$latitude'}, ...
};

cpts = cell(1,length(fields));
for i=1:length(fields)
    cpts{i} = [-1e10 fields{i}{4}(1):fields{i}{5}:fields{i}{4}(2) 1e10];
end
npts = length(cpts{1});
clab = [1 1 1 1 1 1];

% $$$ cmap = redblue(npts-3);
cmap = parula(npts-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

AIsp = 0.25;
latfilt = 1;

%Fluxes only:
% $$$ figure;
set(gcf,'Position',[2125          11        1680         960]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

% 2x3:
% $$$ poss = [0.11     0.5949    0.25      0.3301; ...
% $$$         0.375    0.5949    0.1875    0.3301; ...
% $$$         0.58     0.5949    0.1719    0.3301; ...
% $$$         0.11     0.2300    0.25      0.3301; ...
% $$$         0.375    0.2300    0.1875    0.3301; ...
% $$$         0.58     0.2300    0.1719    0.3301];
poss = [0.11     0.5949    0.21    0.3301; ...
        0.35    0.5949    0.21    0.3301; ...
        0.59     0.5949    0.21    0.3301; ...
        0.11     0.2300    0.21    0.3301; ...
        0.35    0.2300    0.21    0.3301; ...
        0.59     0.2300    0.21    0.3301];
        
letlabs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
for i=1:length(fields)
    subplot(2,3,rr);
    [X,Y] = ndgrid(yt,Te);

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

% $$$         if (i>=1)
% $$$             eval(['VAR = ZA_' regLets{reg} '.AI/1e15;']);        
% $$$             eval(['VAR(ZA_' regLets{reg} '.NaNst==1) = NaN;']);
% $$$             VAR = filter_field(VAR',latfilt,'-t')';
% $$$             [c,h] = contour(X,Y,VAR,[-50:AIsp:-AIsp],'--k');
% $$$             if (clab(i))
% $$$                 clabel(c,h);
% $$$             end
% $$$             [c,h] =  contour(X,Y,VAR,[AIsp:AIsp:50],'-k');
% $$$             if (clab(i))
% $$$                 clabel(c,h);
% $$$             end
% $$$             
% $$$ % $$$             if (i==2 & r == 3)
% $$$ % $$$             % 0-contour on Indo-Pacific mixing:
% $$$ % $$$             VAR(X>39) = 0.25;
% $$$ % $$$             [c,h] = contour(X,Y,VAR,[0 0],'--','color',[0 0.5 0]);
% $$$ % $$$             clabel(c,h);
% $$$ % $$$             end
% $$$         end
        
        ylim([-3 34]);
        caxis(fields{i}{4});
        box on; 
        grid on;
        letno = rr;%3*(i-1)+r;
        xlim([-80 80]);
        text(-79,32.15,label);%[strrep(model,'_','\_')]);%RUNS{letlabs{letno} ' ' fields{i}{3}]);
        set(gca,'xtick',[-90:30:90]);
        if (rr>=4)
            xlabel('Latitude ($^\circ$N)');
        end
        if (rr==1 | rr == 4)
            ylabel('Temperature $\Theta$ ($^\circ$C)');
        end
        if (rr== 3 | rr == 6)
            cb = colorbar;
            ylabel(cb,fields{i}{6});
        end
    end
    set(gca,'Position',poss(rr,:));
colormap(cmap);
end
