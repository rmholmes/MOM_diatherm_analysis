% This script makes plots depth-longitude or depth-latitude slices 

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[4567]}, ...
% $$$ % MOM025-SIS:
% $$$     {'MOM025_kb3seg',[101120]}, ...
% $$$     {'MOM025',[15:19]}, ...
% $$$     {'MOM025_kb1em5',[95:99]}, ...
% $$$     {'MOM025_kb1em6',[30]}, ...
% $$$ % ACCESS-OM2 Gadi runs:
         {'ACCESS-OM2_1deg_jra55_ryf',[31]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_sgl',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf',[51]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[80]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[300]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
% $$$          {'ACCESS-OM2_01deg_jra55_ryf',[636639]}, ...
       };


rr = 1;
for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model;
    
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

vars = {'temp','mld','ndif','vdif','vnlc'};%,'u_sq','v_sq','u','v','w_sq','w','Tdxsq','Tdysq','Tdzsq'};
% $$$ for i=1:length(vars)
% $$$     eval(['sz = size(' vars{i} ');']);
% $$$     if (~anavg)
% $$$         sz(end) = 12;
% $$$     end
% $$$     eval([vars{i} 'all = reshape(' vars{i} ',[sz nyrs]);']);
% $$$ end
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf(['_output%03d_varsat_' reg '.mat'],outputs(i))]);
% $$$     for i=1:length(vars)
% $$$         eval(['sz = size(' vars{i} ');']);
% $$$         if (~anavg)
% $$$             sz(end) = 12;
% $$$         end
% $$$         eval([vars{i} 'all = cat(4,' vars{i} 'all,reshape(' vars{i} ',[sz nyrs]));']);
% $$$     end
% $$$ end
% $$$ for i=1:length(vars)
% $$$     eval(['sz = size(' vars{i} 'all);']);
% $$$     if (~anavg)
% $$$         eval([vars{i} ' = mean(' vars{i} 'all,length(sz));']);
% $$$     else
% $$$         eval([vars{i} ' = ' vars{i} 'all;']);
% $$$     end
% $$$     eval(['clear ' vars{i} 'all;']);
% $$$ end

[xL,zL,tL] = size(temp);
TL = length(T);

if (max(max(max(temp)))>100)
    temp = temp-273.15;
end

% $$$ %% Plot Temp bias against WOA13:
% $$$ months = [1:12];
% $$$ temp = monmean(temp(:,:,months),3,ndays(months));
% $$$ 
% $$$ % WOA13:
% $$$ WOAname = '/srv/ccrc/data03/z3500785/WOA13/woa13_decav_t00_04v2.nc';
% $$$ WOAlat = ncread(WOAname,'lat');
% $$$ WOAlon = ncread(WOAname,'lon');
% $$$ WOAdep = ncread(WOAname,'depth');
% $$$ % $$$ [tmp Eqind] = min(abs(WOAlat));
% $$$ [tmp ln140ind] = min(abs(WOAlon+110));
% $$$ 
% $$$ % $$$ WOAT = squeeze(ncread(WOAname,'t_an',[1 Eqind 1 1],[1440 1 102 1]));
% $$$ WOAT = squeeze(ncread(WOAname,'t_an',[ln140ind 1 1 1],[1 720 102 1]));
% $$$ % $$$ 
% $$$ % $$$ %Shift longitudes:
% $$$ % $$$ [tmp ind] = min(abs(WOAlon-80));
% $$$ % $$$ WOAT = cat(1,WOAT(ind+1:end,:),WOAT(1:ind,:));
% $$$ % $$$ WOAlon = cat(1,WOAlon(ind+1:end)-360,WOAlon(1:ind));
% $$$ % $$$ [WOAlon,WOAdep] = ndgrid(WOAlon,WOAdep);
% $$$ [WOAlat,WOAdep] = ndgrid(WOAlat,WOAdep);
% $$$ % $$$ 
% $$$ % Calculate bias from WOA:
% $$$ % $$$ Tbias = temp-interp2(WOAlon',-WOAdep',WOAT',Xt,-Zt,'linear');
% $$$ Tbias = temp-interp2(WOAlat',-WOAdep',WOAT',Yt,-Zt,'linear');
% $$$ 
% $$$ %Colormap:
% $$$ clim = [-5 5];
% $$$ sp = 0.25;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$     
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ subplot(2,3,rr);
% $$$ % $$$ contourf(Xt,-Zt,Tbias,cpts,'linestyle','none');
% $$$ contourf(Yt,-Zt,Tbias,cpts,'linestyle','none');
% $$$ hold on;
% $$$ % $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[0:2:35],'-k');
% $$$ [c,h] = contour(WOAlat,-WOAdep,WOAT,[0:2:35],'-k');
% $$$ clabel(c,h,[0:2:35]);
% $$$ % $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[20 20],'-k','linewidth',2);
% $$$ [c,h] = contour(WOAlat,-WOAdep,WOAT,[20 20],'-k','linewidth',2);
% $$$ hold on;
% $$$ % $$$ [c,h] = contour(Xt,-Zt,temp,[20 20],'--k','linewidth',2);
% $$$ [c,h] = contour(Yt,-Zt,temp,[20 20],'--k','linewidth',2);
% $$$ % $$$ ylim([-250 0]);
% $$$ % $$$ xlim([-220 -80]);
% $$$ ylim([-300 0]);
% $$$ xlim([-15 15]);
% $$$ % $$$ if (rr == 3)
% $$$     cb = colorbar;
% $$$ % $$$ end
% $$$ % $$$ xlabel('Longitude ($^\circ$E)');
% $$$ xlabel('Latitude ($^\circ$N)');
% $$$ ylabel('Depth (m)');
% $$$ caxis(clim);
% $$$ colormap(cmap);
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ' ...
% $$$                     ,'AOM'),'deg jra55',''),' may',''),'ryf8485 ','') ' - WOA13 Equatorial T ($^\circ$C)']);
% $$$ 
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

var = cumsum(vdif+vnlc,2,'reverse'); % Vertical Mixing Flux
var = ndif; % Numerical mixing
% $$$ var = u_sq - u.^2 + v_sq-v.^2; % EKE
% $$$ var = w_sq - w.^2; % Vertical EKE
% $$$ var = Tdxsq+Tdysq; % Horizontal T differences
% $$$ var = Tdzsq;

months = {[1:tL]}%:12]};

% $$$ var = mean(reshape(var,[163 74 12 10]),4);
% $$$ months = {[1:12]}%:12]};
% $$$ months = {[1:12],[3],[7],[11]};
% $$$ monthsu01 = {[1:4],[1],[3],[4]};
% $$$ labels = {'Annual','March','July','November'};

    %Colormap and continents:
    sp = 0.01;
    clim = [0 0.3];
    sp = 0.1e-8;
    clim = [0 3e-8];
    sp = 0.05e-9;
    clim = [0 1e-9];
    sp = 0.5;
    clim = [-20 0];
% $$$     sp = 0.5;
% $$$     clim = [-20 0];
% $$$     sp = 1;
% $$$     clim = [-30 0];
% $$$     sp = 0.01;
% $$$     clim = [0 0.3];
% $$$     sp = 1e-9;
% $$$     clim = [0 3e-8];
% $$$     sp = 0.05e-9;
% $$$     clim = [0 1e-9];
% $$$     sp = 0.2;
% $$$     clim = [0 10];

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
        buf = 12;
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
% $$$     labels = {'(a) Numerical Mixing','(b) Vertical Mixing'};%,'(c) KDS75','(d) KDS100','(e) KDS135'};
% $$$     labels = {strrep(RUNS{rr}{1},'_',' ')}
% $$$     poss = [0.1300    0.4800    0.4154    0.3355; ...
% $$$             0.1300    0.1100    0.4154    0.3355;];    
    % ACCESS-OM2 vertical res:
    labels = {'(a) KDS50','(b) GFDL50','(c) KDS75','(d) KDS100','(e) KDS135','(f) KDS50 - $\Delta t=1.5$hr'};
    poss = [0.0886    0.7160    0.3375    0.2581; ...
            0.4503    0.7160    0.3375    0.2581; ...
            0.0886    0.4163    0.3375    0.2581; ...
            0.4503    0.4163    0.3375    0.2581; ...
            0.0886    0.1142    0.3375    0.2581];
    % ACCESS-OM2 extras:
    labels = {'(a) ACCESS-OM2-1-KDS50','(b) ACCESS-OM2-1-KDS50 with $\Delta t=1.5$hr','(c) ACCESS-OM2-025-RG','(d) $1/4^\circ$','(e) $1/4^\circ$ R'};
            
% $$$     
% $$$     % MOM025 Control dif vars:
% $$$     labels = {'(a) $\overline{u''u''}+\overline{v''v''}$', ...
% $$$               '(b) $\overline{w''w''}$','(c) $|\Delta_x T|^2 + |\Delta_y T|^2$','(d) $|\Delta_z T|^2$'};
% $$$     units = {'$m^2s^{-2}$','$m^2s^{-2}$','$^\circ C^2$','$^\circ C^2$'};
% $$$     poss = [0.05500    0.48100    0.4   0.3355; ...
% $$$            0.53500    0.48100    0.4    0.3355; ...
% $$$            0.05500    0.1100    0.4     0.3355; ...
% $$$            0.53500    0.1100    0.4     0.3355;];    

% $$$ figure;
set(gcf,'Position',[1923           5        1366         998]);
% $$$ set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
% $$$ rr = 2;
for i=1:length(months)
    subplot(3,2,rr);
% $$$     subplot(3,3,9);%rr);
    contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
% $$$     contourf(Xu,-Zu,nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
    hold on;
    Tout = nanmonmean(temp(:,:,months{i}),3,ndays(months{i}));
    Tout(Tout==0) = NaN;
    [c,h] = contour(Xt,-Zt,Tout,[0:1:35],'-k');
    clabel(c,h,[0:2:35]);
    [c,h] = contour(Xt,-Zt,Tout,[22.5 22.5],'-k','linewidth',2);
% $$$ if (strcmp(model,'MOM01'))
% $$$     mnu = monthsu01{i};
% $$$ else
% $$$     mnu = months{i};
% $$$ end
% $$$ ucol = [0.8706    0.4902         0];
% $$$ [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
% $$$                 'color',ucol);
% $$$ [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
% $$$                 'color',ucol);
% $$$ clabel(c,h,'color','w');
    plot(Xu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
    ylim([-200 0]);
% $$$     xlim([-80 -12]);
    xlim([-200 -80]);
    if (rr == 2 | rr == 4 | rr == 5)
% $$$     if (rr >= 1)
        cb = colorbar;
% $$$         ylabel(cb,units{rr});
        ylabel(cb,'Wm$^{-2}$');
    end
    if (rr >=2)
        xlabel('Longitude ($^\circ$E)');
    end
    if (rr==1 | rr == 3 | rr == 5)
        ylabel('Depth (m)');
    end
    if (rr == 2 | rr == 4)
        set(gca,'yticklabel',[]);
    end
    if (rr <= 1)
        set(gca,'xticklabel',[]);
    end
    caxis(clim);
    text(-199,-15,labels{rr},'Backgroundcolor','w','FontSize',15,'margin',0.5);
    set(gca,'Position',poss(rr,:));
% $$$ title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ' ...
% $$$                     ,'AOM'),' jra55',''),' may',''),'ryf8485 ','') ...
% $$$        ' Vertical Mixing']);
end

colormap(gca, cmap);
%colormap(cmap);
end

end


%%% Latitudinal Slices:

% Load Base Variables:
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
outputs = [75:79];
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
ndays = diff(time_snap);

% Load Variable and calculate mean:
lonsl = 110;
load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(1))]);
vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd'};
for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
end
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(i))]);
    for i=1:length(vars)
        eval([vars{i} 'a = ' vars{i} 'a + ' vars{i} ';']);
    end
end
for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a/length(outputs);']);
    eval(['clear ' vars{i} 'a;']);
end
if(max(max(max(temp)))>100)
    temp = temp-273.15;
end

[yL,zL,tL] = size(temp);
TL = length(T);

% Depth of isotherms:
Zi = zeros(yL,TL,tL);
for ti=1:tL
    for yi=1:yL
        tvec = squeeze(temp(yi,:,ti));
        zvec = -Zt(yi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(yi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(yi,:,ti)),1,'last');
        Zi(yi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Yi = repmat(Yt(:,1),[1 TL]);

var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux

%Colormap and continents:
sp = 5;
clim = [-150 0];

cCH = 1; % 0 = symmetric redblue
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
    buf = 6;
    clim = [clim(1) buf*sp];
    cpts = [-1e10 clim(1):sp:clim(2) 1e10];
    cmap(end+1,:) = cmap(end,:); % 1st positive bin
    cmap(end+buf-1,:) = [1 0.7 0.9]; % last pink bin
    for ii = 1:(buf-2)
        cmap(end-buf+1+ii,:) = cmap(end,:)*ii/(buf-1) + ...
            cmap(end-buf+1,:)*(buf-1-ii)/(buf-1);
    end
end        

% $$$ var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
% $$$ var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
% $$$ clim = [-1 1]*1e-5*86400; % FOR WMT
% $$$ sp = 0.1*1e-5*86400;
% $$$ doWMT = 1;

months = {[1:12],[3],[8]};
monthsu01 = {[1:4],[1],[3],[4]};
labels = {'Annual','March','August'};

% $$$ %Save for schematic:
% $$$ Xl = Yi;
% $$$ Yl = nanmonmean(Zi(:,:,months{i}),3,ndays(months{i}));
% $$$ Zl = nanmonmean(vdif(:,:,months{i}),3,ndays(months{i}));
% $$$ XlC = Yt;
% $$$ YlC = -Zt;
% $$$ ZlC = monmean(temp(:,:,months{i}),3,ndays(months{i}));

% $$$ [tmp Eqind] = min(abs(Yu(:,1)));
% $$$ tauweight = abs(taux(Eqind,:))*200;
% $$$ % Wind stress vectors:
% $$$ sp  =5;
% $$$ yvec = Yu(:,1);
figure;
set(gcf,'Position',[1          36        1920         970]);
set(gcf,'defaulttextfontsize',20);
set(gcf,'defaultaxesfontsize',20);
for i=1:length(months)
subplot(3,1,i);
contourf(Yi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
[c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
clabel(c,h,[0:2:35]);
[c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3, ...
                               ndays(months{i})),[21.5 21.5],'-k','linewidth',2);
if (strcmp(model,'MOM01'))
    mnu = monthsu01{i};
else
    mnu = months{i};
end
ucol = [0.8706    0.4902         0];
[c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
                'color',ucol);
[c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
                'color',ucol);
plot(Yu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
ylim([-200 0]);
%ylim([-250 0]);
xlim([-10 10]);
cb = colorbar;
if (doWMT)
    ylabel(cb,'m/day');
else
    ylabel(cb,'Wm$^{-2}$');
end
xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
caxis(clim);
text(-9.6,-188,labels{i},'Backgroundcolor','w','FontSize',20);
text(9.6,-188,[num2str(lonsl) '$^\circ$W'],'Backgroundcolor','w','FontSize',20,'HorizontalAlignment','Right');

% $$$ %Add wind-stress vectors:
% $$$ pos = get(gca,'Position')
% $$$ wsh = axes('Position',[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]);
% $$$ % $$$ for ii=1:sp:length(yvec)
% $$$     plot(0,0,'o','MarkerSize',abs(mean(mean(taux(yvec>=-10 & yvec<=10,months{i}),2),1))*200);
% $$$     hold on;
% $$$ % $$$ end
% $$$ %quiver(yvec,zeros(size(yvec)),mean(taux(1:sp:end,months{i}),2),mean(tauy(1:sp:end,months{i}),2));
% $$$ xlim([-10 10]);
% $$$ ylim([-1 1]);
% $$$ box off;axis off;
% $$$ set(wsh,'Position',[[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]]);

LabelAxes(gca,i,20,0.008,0.95);
end
colormap(cmap);

