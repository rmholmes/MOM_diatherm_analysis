% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';


RUNS = { ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680],'ACCESS-OM2-025'}, ...
       };
doZAremap = 1; % remap too depth.
subCont = 1; % subtract control and plot anomalies.

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
        if (doZAremap)
            ZAtempS(:,:,i) = ZAtemp;
            tempZAS(:,:,i) = tempZA;
            rhoZAS(:,:,i) = rhoZA;
        end
    end
    
    ZAtemp = mean(ZAtempS,3);
    tempZA = mean(tempZAS,3);
    rhoZA = mean(rhoZAS,3);

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

%% Plot latitude - temperature plane for different basins:
dy = diff(yu);
dy = [dy(1); dy];

% NumMix:
fields = { ...
          {'NUM',1./repmat(dy,[1 TL+1])/1e12,'Numerical Mixing',[-25 0],0.5,'TW/$^\circ$latitude'},
          {'M',-1./repmat(dy,[1 TL+1])/1e12,'Vertical Mixing',[-25 0],0.5,'TW/$^\circ$latitude'},
};

cpts = cell(1,length(fields));
for i=1:length(fields)
    cpts{i} = [-1e10 fields{i}{4}(1):fields{i}{5}:fields{i}{4}(2) 1e10];
end
npts = length(cpts{1});
clab = [1 1 1 1 1 1];

cmap = parula(npts-3);
cmap(end,:) = [0.97 0.97 0.8];
cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;

AIsp = 0.25;
latfilt = 1;

set(gcf,'Position',[3          40        1278         963]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

poss = [0.12    0.39      0.35    0.25; ...
        0.5210    0.39      0.35    0.25];
        
labels = {'(a) Numerical Mixing','(b) Vertical Mixing'};
for i=1:length(fields)
    subplot(3,2,2*(rr-1)+i);
    if (doZAremap)    
        X =repmat(yt,[1 TL+1]);
        Y = ZAtemp;
    else
        [X,Y] = ndgrid(yt,Te);
    end

    eval(['VAR = ZA_' regLets{reg} '.' fields{i}{1} '.*fields{i}{2};']);
    if (subCont)
    if (rr==1)
        VARC = VAR;
    else
        VAR = VAR-VARC;
    end
    end
    VAR(VAR==0) = NaN;
    eval(['VAR(ZA_' regLets{reg} '.NaNst==1) = NaN;']);
    VAR = filter_field(VAR',latfilt,'-t')';
    contourf(X,Y,VAR,cpts{i},'linestyle','none');
    hold on;
    col = [0 0 0];
    if (doZAremap)
        [X,Y] = ndgrid(yt,-z);
        [c,h] = contour(X,Y,tempZA,[-2:2:34],'-k');
        clabel(c,h);
        [c,h] = contour(X,Y,tempZA,[5 15 22.5],'-m','linewidth',1);
    else
        eval(['plot(yt,filter_field(ZA_' regLets{reg} '.maxTt,latfilt,''-t''),'':k'');']);
    end
    ylim([-4000 -500]);
    text(-79,-3800,labels{i},'BackgroundColor','w');%[strrep(model,'_','\_')]);%RUNS{letlabs{letno} ' ' fields{i}{3}]);
    caxis([-25 0]);%fields{i}{4});
    box on; 
    grid on;
    letno = rr;%3*(i-1)+r;
    xlim([-80 80]);
    set(gca,'xtick',[-90:30:90]);
    xlabel('Latitude ($^\circ$N)');
    if (i == 1)
        ylabel('Depth (m)');
    else
        set(gca,'yticklabel',[]);
        cb = colorbar;
        ylabel(cb,fields{i}{6});
    end
    set(gca,'color','k');

    set(gca,'Position',poss(i,:));
    
    gca2 = copyobj(gca,gcf);
    axes(gca2);
    colorbar off;
    pos = poss(i,:);
    pos(2) = pos(2)+pos(4);
    pos(4) = pos(4)/4*3;
    set(gca2,'Position',pos);
    set(gca2,'xticklabel',[]);
    xlabel(gca2,'');
    ylabel(gca2,'');
    ylim([-500 0]);
    
end
colormap(cmap);
end

