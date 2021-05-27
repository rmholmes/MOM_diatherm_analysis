% This script makes plots depth-longitude or depth-latitude slices 

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
% $$$          {'ACCESS-OM2_01deg_jra55_ryf_k_smag_iso7p9',[648655]}, ...
       };
cols = {'b',[0.3020    0.7451    0.9333],'b','b','b'};
typs = {'-','-','--',':','-.'};

figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

for rr = 1:length(RUNS);
% $$$ rr = 1;
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

    % Load Variable and calculate mean:
    reg = 'EqPM2';
    load([base model sprintf(['_output%03d_varsat_' reg '.mat'],outputs(1))]);

    vars = {'temp','mld','ndif','vdif','vnlc'};

    [xL,zL,tL] = size(temp);
    TL = length(T);

    if (max(max(max(temp)))>100)
        temp = temp-273.15;
    end

% $$$ %% Plot Temp bias against WOA13:
% $$$ % $$$ months = [1:12];
% $$$ % $$$ temp = monmean(temp(:,:,months),3,ndays(months));
% $$$ months = 1:10;
% $$$ temp = mean(temp,3);
% $$$ 
% $$$ % WOA13:
% $$$ WOAname = '/srv/ccrc/data03/z3500785/Data_Products/WOA13/woa13_decav_t00_04.nc';
% $$$ WOAlat = ncread(WOAname,'lat');
% $$$ WOAlon = ncread(WOAname,'lon');
% $$$ WOAdep = ncread(WOAname,'depth');
% $$$ [tmp Eqind] = min(abs(WOAlat));
% $$$ % $$$ [tmp ln140ind] = min(abs(WOAlon+110));
% $$$ 
% $$$ WOAT = squeeze(ncread(WOAname,'t_an',[1 Eqind 1 1],[1440 1 102 1]));
% $$$ % $$$ WOAT = squeeze(ncread(WOAname,'t_an',[ln140ind 1 1 1],[1 720 102 1]));
% $$$ % $$$ 
% $$$ %Shift longitudes:
% $$$ [tmp ind] = min(abs(WOAlon-80));
% $$$ WOAT = cat(1,WOAT(ind+1:end,:),WOAT(1:ind,:));
% $$$ WOAlon = cat(1,WOAlon(ind+1:end)-360,WOAlon(1:ind));
% $$$ [WOAlon,WOAdep] = ndgrid(WOAlon,WOAdep);
% $$$ % $$$ [WOAlat,WOAdep] = ndgrid(WOAlat,WOAdep);
% $$$ % $$$ 
% $$$ % Calculate bias from WOA:
% $$$ Tbias = temp-interp2(WOAlon',-WOAdep',WOAT',Xt,-Zt,'linear');
% $$$ % $$$ Tbias = temp-interp2(WOAlat',-WOAdep',WOAT',Yt,-Zt,'linear');
% $$$ 
% $$$ %Colormap:
% $$$ clim = [-3 3];
% $$$ sp = 0.25;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$     
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ subplot(2,3,rr);
% $$$ contourf(Xt,-Zt,Tbias,cpts,'linestyle','none');
% $$$ % $$$ contourf(Yt,-Zt,Tbias,cpts,'linestyle','none');
% $$$ hold on;
% $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[0:2:35],'-k');
% $$$ % $$$ [c,h] = contour(WOAlat,-WOAdep,WOAT,[0:2:35],'-k');
% $$$ clabel(c,h,[0:2:35]);
% $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[20 20],'-k','linewidth',2);
% $$$ % $$$ [c,h] = contour(WOAlat,-WOAdep,WOAT,[20 20],'-k','linewidth',2);
% $$$ hold on;
% $$$ [c,h] = contour(Xt,-Zt,temp,[20 20],'--k','linewidth',2);
% $$$ % $$$ [c,h] = contour(Yt,-Zt,temp,[20 20],'--k','linewidth',2);
% $$$ ylim([-300 0]);
% $$$ xlim([-220 -80]);
% $$$ % $$$ ylim([-300 0]);
% $$$ % $$$ xlim([-15 15]);
% $$$ if (rr == 3 | rr == 6)
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature Bias ($^\circ$C)');
% $$$ end
% $$$ if (rr <= 3)
% $$$     set(gca,'xticklabel',[]);
% $$$ else
% $$$     xlabel('Longitude ($^\circ$E)');
% $$$ end
% $$$ % $$$ end
% $$$ % $$$ xlabel('Latitude ($^\circ$N)');
% $$$ if (rr == 1 | rr == 4)
% $$$     ylabel('Depth (m)');
% $$$ else
% $$$     set(gca,'yticklabel',[]);
% $$$ end
% $$$ caxis(clim);
% $$$ colormap(cmap);
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ % $$$ title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ' ...
% $$$ % $$$                     ,'AOM'),'deg jra55',''),' may',''),'ryf8485 ','') ' - WOA13 Equatorial T ($^\circ$C)']);
% $$$ text(-219,-285,RUNS{rr}{3},'Backgroundcolor','w');
% $$$ 
% $$$ poss = [0.0693    0.5271    0.2726    0.3753; ...
% $$$         0.3578    0.5271    0.2726    0.3753; ...
% $$$         0.6463    0.5271    0.2726    0.3753; ...
% $$$         0.0693    0.1134    0.2726    0.3753; ...
% $$$         0.3578    0.1134    0.2726    0.3753; ...
% $$$         0.6463    0.1134    0.2726    0.3753];
% $$$ set(gca,'Position',poss(rr,:));
% $$$ end

%%% Plot Diathermal fluxes:

% Depth of isotherms:
Zi = zeros(xL,TL,tL);
for ti=1:tL
    for xi=1:xL
        tvec = squeeze(temp(xi,:,ti));
        zvec = -Zt(xi,:);
        tvec(isnan(tvec)) = -1000;
        tvec(tvec == 0) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(xi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(xi,:,ti)),1,'last');
        Zi(xi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Xi = repmat(Xt(:,1),[1 TL]);

% Choose variable:
var = cumsum(vdif+vnlc,2,'reverse'); % Vertical Mixing Flux
var = ndif; % Numerical mixing

months = {[1:tL]}%:12]};

    %Colormap and continents:
    sp = 2;
    clim = [-30 0];

    cCH = 2; % 0 = symmetric redblue
             % 1 = negative definite parula
             % 2 = negative parula with +ve's possible
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
        onebin = 0; % only allow one pink bin
        clim = [clim(1) buf*sp];
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        if (onebin)
            for ii = 1:buf
                cmap(end+1,:) = cmap(npts-3,:);
            end
            cmap(end,:) = [1 0.7 0.9]; % last pink bin
        else
            cmap(end+1,:) = cmap(end,:); % 1st positive bin
            cmap(end+buf-1,:) = [1 0.7 0.9]; % last pink bin
            for ii = 1:(buf-2)
                cmap(end-buf+1+ii,:) = cmap(end,:)*ii/(buf-1) + ...
                    cmap(end-buf+1,:)*(buf-1-ii)/(buf-1);
            end
        end
    end        

% $$$ rr = 1;
for i=1:length(months)
    subplot(3,4,rr);
    contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
    hold on;
    Tout = nanmonmean(temp(:,:,months{i}),3,ndays(months{i}));
    Tout(Tout==0) = NaN;
    [c,h] = contour(Xt,-Zt,Tout,[0:1:35],'-k');
    clabel(c,h,[0:2:35]);
    ylim([-200 0]);
    xlim([-200 -80]);
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
    xlabel('Longitude ($^\circ$E)');
    ylabel('Depth (m)');
    caxis(clim);
    title(strrep(RUNS{rr}{1},'_',' '))
end

colormap(gca, cmap);
%colormap(cmap);
end

