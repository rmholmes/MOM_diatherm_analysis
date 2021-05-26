% This script makes plots depth-longitude or depth-latitude slices 

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[81:82]}, ...
       };
cols = {'b',[0.3020    0.7451    0.9333],'b','b','b'};
typs = {'-','-','--',':','-.'};

% $$$ figure;
% $$$ set(gcf,'Position',[1923           5        1366         998]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);

rr = 1;
% $$$ for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

% $$$     clearvars -except base RUNS rr outputs model;
    
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
reg = 'GSLS';
load([base model sprintf(['_output%03d_varsat_' reg '.mat'],outputs(1))]);

vars = {'temp','mld','ndif','vdif','vnlc','aiso_bih'};%,'u_sq','v_sq','u','v','w_sq','w','Tdxsq','Tdysq','Tdzsq'};
for i=1:length(vars)
    eval([vars{i} 'all = ' vars{i} ';']);
end

for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_' reg '.mat'],outputs(i))]);
    for i=1:length(vars)
        eval(['sz = size(' vars{i} ');']);
        eval([vars{i} 'all = cat(length(sz),' vars{i} 'all,' vars{i} ');']);
    end
end
for i=1:length(vars)
    eval(['sz = size(' vars{i} ');']);
    eval([vars{i} ' = mean(' vars{i} 'all,length(sz));']);
end

[xL,zL] = size(temp);
TL = length(T);

if (max(max(max(temp)))>100)
    temp = temp-273.15;
end

%%% Plot Diathermal fluxes:

tL = 1;
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

M = cumsum(vdif+vnlc,2,'reverse'); % Vertical Mixing Flux
I = ndif;
sp = 2;
clim = [-150 0];

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
% $$$         buf = 12;
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
% $$$     cmap = flipud(cmap);
    
    % MOM025 kb3seg example:
    labels = {'(a) $\mathcal{I}$, $10^\circ$C isotherm','(b) Numerical Mixing',...
             '(c) Vertical Mixing','(d) Biharmonic viscosity'};
    
    text(-118,61,labels{1},'Backgroundcolor','w','FontSize',15,'margin',0.5);
    hold on;
    plot([-40 -40 -50 -50 -40],[35 62 62 35 35],'-k','linewidth',2);
    
    subplot(2,2,2);
    contourf(Xi,Zi,I,cpts,'linestyle','none');
    hold on;
    [c,h] = contour(Xt,-Zt,temp,[0:1:35],'-k');
    clabel(c,h,[0:2:35]);
    ylim([-1500 0]);
    xlim([35 62]);
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
    xlabel('Latitude ($^\circ$N)');
    ylabel('Depth (m)');
    caxis(clim);
    text(35.1,-1440,labels{2},'Backgroundcolor','w','FontSize',15,'margin',0.5);
    colormap(gca, cmap);
    
    subplot(2,2,3);
    contourf(Xi,Zi,M,cpts,'linestyle','none');
    hold on;
    [c,h] = contour(Xt,-Zt,temp,[0:1:35],'-k');
    clabel(c,h,[0:2:35]);
    ylim([-1500 0]);
    xlim([35 62]);
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
    xlabel('Latitude ($^\circ$N)');
    ylabel('Depth (m)');
    caxis(clim);
    text(35.1,-1440,labels{3},'Backgroundcolor','w','FontSize',15,'margin',0.5);
    colormap(gca, cmap);

    subplot(2,2,4);
    contourf(Xt,-Zt,aiso_bih,[0:1e9:5e10 1e50],'linestyle','none');
    hold on;
    [c,h] = contour(Xt,-Zt,temp,[0:1:35],'-k');
    clabel(c,h,[0:2:35]);
    ylim([-1500 0]);
    xlim([35 62]);
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
    xlabel('Latitude ($^\circ$N)');
    ylabel('Depth (m)');
    text(35.1,-1440,labels{4},'Backgroundcolor','w','FontSize',15,'margin',0.5);

