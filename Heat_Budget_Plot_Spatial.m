% This script makes plots of the spatial structure of the
% diathermal fluxes in the MOM simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
% ACCESS-OM2 Gadi runs:
% $$$ % 1-degree
         {'ACCESS-OM2_1deg_jra55_ryf',[31]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135]}, ...
% 1/4-degree
         {'ACCESS-OM2_025deg_jra55_ryf',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5',[7781]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kbvar',[7781]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_kds75',[7680]}, ...
% 1/10-degree
         {'ACCESS-OM2_01deg_jra55_ryf',[636643]}, ...
         {'ACCESS-OM2_01deg_jra55_ryf_k_smag_iso7p9',[648:655]}, ...
       };
cols = {'b','r','k','m','g'};

rr = 7;
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
        ndays = ndays(1:12);
    end
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
    nyrs = tL/12;
    if (round(nyrs)~=nyrs)
        anavg = 1;
        nyrs = tL;
        months = {[1:1]};
    else
        anavg = 0;
        months = {[1:12]};
    end
    yrs = 1:nyrs;
    ycur = 1;

    TYPE = 'VertInt';
    Tls = 22.5; % Isotherm choice
    labels = {'(a) $22.5^\circ$C'};
    iii = 1;
    name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
    load(name,'FlI','FlM');
    % Choose vertical or numerical mixing:
    FlM = FlI;
    FlM(FlM==0) = NaN;

    try
        obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
        LAND = obj.SST(:,:,1);
    catch
        LAND = zeros(size(FlM(:,:,1)));
    end

    %If MOM01, fix NaN's in grid:
    if (strfind(model,'01'))
        lon = repmat(lon(:,500),[1 yL]);
        latv = nanmean(lat,1);
        lat = repmat(latv,[xL 1]);
    end

    [xL,yL] = size(lon);
    xvec = 1:2:xL;
    yvec = 1:2:yL;
    txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    
% $$$     %Colormap and continents:
    sp = 5;
    clim = [-100 0];
    cCH = 2; % 0 = symmetric redblue
             % 1 = negative definite parula
             % 2 = negative parula with +ve's possible
             % 3 = negative parula with cutoff pink above
    if (cCH==0)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = redblue(npts-3);
        for i=1:(npts-3)
            if (cmap(i,:) == 1.0)
                cmap(i,:) = [0.94 0.94 0.94];
            end
        end
    elseif (cCH>=1)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = parula(npts-3);
        cmap(end,:) = [0.97 0.97 0.8];
        cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
    end
    if (cCH == 2)
        buf = 2;
        clim = [clim(1) buf*sp];
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        cmap(end+1,:) = cmap(end,:); % 1st positive bin
        cmap(end+buf-1,:) = [1 0.7 0.9]; % last pink bin
        for ii = 1:(buf-2)
            cmap(end-buf+1+ii,:) = cmap(end,:)*ii/(buf-1) + ...
                cmap(end-buf+1,:)*(buf-1-ii)/(buf-1);
        end
    end

    tmp = LAND;
    tmp(isnan(LAND)) = clim(1)-sp/2;
    tmp(~isnan(LAND)) = NaN;
    LAND = tmp;
    cmap(2:(end+1),:) = cmap;
    cmap(1,:) = [0 0 0];

    climn = [clim(1)-sp clim(2)];
    figure;
    set(gcf,'Position',[1          36        1920         970]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);

    poss = [0.1300    0.54    0.7403    0.4149; ...
            0.1300    0.0876    0.7403    0.4149];
    i = 1;
    ii = 2;
    subplot(2,1,ii);
% $$$     subplot(3,2,2*(rr-1)+iii);
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

    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z,cpts,'linestyle','none');    
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    cb = colorbar;
    ylabel(cb,'$\mathcal{I}$ (Wm$^{-2}$)');
    ylim(cb,clim);
    hold on;
    if (ii>=3)
        xlabel('Longitude ($^\circ$E)');
    else
        set(gca,'xticklabel',[]);
    end
    if (ii == 1 | ii ==3)
        ylabel('Latitude ($^\circ$N)');
    else
        set(gca,'yticklabel',[]);
    end
        
    set(gca,'xtick',[-270:30:60]);
    set(gca,'ytick',[-75:15:75]);
    ylim([-65 75]);
    text(-277,68,labels{ii},'BackgroundColor','w','Margin',0.5,'FontSize',13);
    colormap(cmap);
