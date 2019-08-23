% This script performs a Helmholtz decomposition on the lateral
% heat fluxes following Forget and Ferriera (2019) Nat Geo.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
gname = '/srv/ccrc/data03/z3500785/mom/ocean_grid_mom025.nc';

RUNS = { ...
    {'MOM025_kb3seg',[101110]}, ...
       };

rr = 1;
outputs = RUNS{rr}{2};
model = RUNS{rr}{1};
    
load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
if (~exist('ndays'))
    ndays = diff(time_snap);
    ndays = ndays(1:12);
end
if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;

% Get t-mask:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'SST');
tDry = isnan(SST(:,:,1));
tWet = ~isnan(SST(:,:,1));
clear SST;

% Get x and y masks:
load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
      strrep(num2str(34),'.','p') 'C.mat']);
xflux(xflux==0) = NaN;
yflux(yflux==0) = NaN;
xDry = isnan(xflux);
yDry = isnan(yflux);
xWet = ~isnan(xflux);
yWet = ~isnan(yflux);

% Load lateral fluxes:
TlT = 34;
TlB = -10;
if (TlB<-5)
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlT),'.','p') 'C.mat']);
else
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlT),'.','p') 'C.mat']);
    xfluxT = xflux;yfluxT=yflux;
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlB),'.','p') 'C.mat']);
    xflux = xfluxT-xflux;
    yflux = yfluxT-yflux;
end
    
xflux(xDry)=0;
yflux(yDry)=0;

% Calculate divergence:
divF = zeros(xL,yL);
divF(2:end,2:end) = divF(2:end,2:end) + (xflux(1:(end-1),2:end) - xflux(2:end,2:end) ...
                                     +yflux(2:end,1:(end-1)) - yflux(2:end,2:end));%./area(2:end,2:end);
divF(1,2:end) = divF(1,2:end) + (xflux(end,2:end) - xflux(1,2:end) ...
                             +yflux(1,1:(end-1)) - yflux(1,2:end));%./area(1,2:end);        
divF(tDry) = NaN;

divF(:,end) = 0;

wetI = find(tWet); nn=length(wetI);

% $$$ A=sparse([],[],[],nn,nn,nn*5);
% $$$ 
% $$$ % Less efficient option, naive tripolar points:
% $$$ for ll=1:nn % result index loop
% $$$ 
% $$$     if (mod(ll,100)==0); display([num2str(ll) ' of ' num2str(nn)]); end;
% $$$     
% $$$     [ii,jj] = ind2sub([xL yL],wetI(ll)); % Get full i,j indices
% $$$     
% $$$     A(ll,ll) = -4;
% $$$     
% $$$     % Build up matrix coefficients:
% $$$     % ii+1, jj:
% $$$     if (ii==xL); lind = sub2ind([xL yL],1,jj); else; lind = sub2ind([xL yL],ii+1,jj); end;
% $$$     
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$     
% $$$ 
% $$$     % ii-1, jj:
% $$$     if (ii==1); lind = sub2ind([xL yL],xL,jj); else; lind = sub2ind([xL yL],ii-1,jj); end;
% $$$     
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$     
% $$$ 
% $$$     % ii, jj+1:
% $$$     if (jj==yL)
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     else
% $$$         lind = sub2ind([xL yL],ii,jj+1);
% $$$         ind = find(wetI==lind);
% $$$         if (length(ind)==1) % source wet point -> Add dependence
% $$$             A(ll,ind) = 1; 
% $$$         else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$         end
% $$$     end
% $$$     
% $$$ 
% $$$     % ii, jj-1:
% $$$     lind = sub2ind([xL yL],ii,jj-1);
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$ 
% $$$ end
% $$$ save('HelmholtzA.mat','A');

% $$$ load('HelmholtzA.mat','A');
% $$$ Aold = A;
% $$$ 
% $$$ % Fix tri-polar jj=nj points:
% $$$ [IIinds,JJinds] = ind2sub([xL yL],wetI); % I and J indices of wet points
% $$$ 
% $$$ inds = find(JJinds==yL);
% $$$ for lli=1:length(inds)
% $$$     ll = inds(lli)
% $$$     [ii,jj] = ind2sub([xL yL],wetI(ll)); % Get full i,j indices    
% $$$ 
% $$$     A(ll,:) = 0*A(ll,:);
% $$$     A(ll,ll) = -4;
% $$$     
% $$$     % Build up matrix coefficients:
% $$$     % ii+1, jj:
% $$$     if (ii==xL); lind = sub2ind([xL yL],1,jj); else; lind = sub2ind([xL yL],ii+1,jj); end;
% $$$     
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$     
% $$$     % ii-1, jj:
% $$$     if (ii==1); lind = sub2ind([xL yL],xL,jj); else; lind = sub2ind([xL yL],ii-1,jj); end;
% $$$     
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$     
% $$$     % ii, jj+1:
% $$$     if (jj~=yL); display('ERROR!'); return; end;
% $$$ 
% $$$     lind = sub2ind([xL yL],xL-ii+1,jj); % T_{i,yL+1} = T_{xL-i+1,yL}
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$ 
% $$$     % ii, jj-1:
% $$$     lind = sub2ind([xL yL],ii,jj-1);
% $$$     ind = find(wetI==lind);
% $$$     if (length(ind)==1) % source wet point -> Add dependence
% $$$         A(ll,ind) = 1; 
% $$$     else
% $$$         A(ll,ll) = A(ll,ll)+1;
% $$$     end
% $$$ end
%save('HelmholtzA_tripolarcorrected.mat','A');
load('HelmholtzA_tripolarcorrected.mat','A');

% Solve it:
b = divF(tWet);
P = -A\b;

% unwrap:
Pout = NaN*zeros(xL,yL);
Pout(tWet) = P-mean(P(:));

xfluxD = zeros(size(Pout));
xfluxD(1:(end-1),:) = Pout(2:end,:)-Pout(1:(end-1),:);
xfluxD(end,:) = Pout(1,:)-Pout(end,:);

yfluxD = zeros(size(Pout));
yfluxD(:,1:(end-1)) = Pout(:,2:end)-Pout(:,1:(end-1));

% $$$ MHTD = nansum(yfluxD,1);
% $$$ MHT = nansum(yflux,1);
% $$$ 
% $$$ % MHT:
% $$$ figure;
% $$$ plot(yu,MHT/1e15);
% $$$ hold on;
% $$$ plot(yu,MHTD/1e15,'--r');

divFch = zeros(xL,yL);
divFch(2:end,2:end) = divFch(2:end,2:end) + (xfluxD(1:(end-1),2:end) - xfluxD(2:end,2:end) ...
                                     +yfluxD(2:end,1:(end-1)) - yfluxD(2:end,2:end));%./area(2:end,2:end);
divFch(1,2:end) = divFch(1,2:end) + (xfluxD(end,2:end) - xfluxD(1,2:end) ...
                             +yfluxD(1,1:(end-1)) - yfluxD(1,2:end));%./area(1,2:end);        
divFch(tDry) = NaN;

% Check plotting:
xvec = 1:1:xL;
yvec = 1:1:yL;

% $$$ % Full heat flux:
% $$$ figure;
% $$$ %set(gcf,'Position',[1923           5        1288         998]);
% $$$ set(gcf,'Position',[2165         236        1288         704]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ %subplot(2,1,1);
% $$$ %set(gca,'Position',[0.0916    0.5472    0.7562    0.3956]);
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),divF(xvec,yvec)./area(xvec,yvec));
% $$$ hold on;
% $$$ xvecV = 1:10:xL;
% $$$ yvecV = 1:10:yL;
% $$$ quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),xflux(xvecV,yvecV),yflux(xvecV,yvecV),8);
% $$$ caxis([-200 200]);
% $$$ title(['Full heat flux convergence (Wm$^{-2}$) and $0^\circ$C-reference ' ...
% $$$        'heat flux vectors ' sprintf('%3.0fC - %3.0fC',TlT,TlB)]);
% $$$ set(gca,'color','k');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ colormap('redblue');
% $$$ ylim([-70 75]);

% Divergent part:
figure;
%set(gcf,'Position',[1923           5        1288         998]);
set(gcf,'Position',[2165         236        1288         704]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
%subplot(2,1,2);
%set(gca,'Position',[0.0916    0.0768    0.7562   0.3956]);
cpts = [-1e10 -200:25:200 1e10];
npts = length(cpts);
contourf(lon(xvec,yvec),lat(xvec,yvec),divFch(xvec,yvec)./ ...
         area(xvec,yvec),cpts,'linestyle','none');
hold on;
minv = min(Pout(:));
maxv = max(Pout(:));
n = 40;
contour(lon,lat,Pout,[-1e50 minv:0.1e15:maxv 1e50],'-k');
xvecV = 1:20:xL;
yvecV = 1:20:yL;
quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),xfluxD(xvecV,yvecV),yfluxD(xvecV,yvecV),2);
caxis([-200 200]);
% $$$ title(['Divergent part convergence (Wm$^{-2}$), gradient function ' ...
% $$$        'and divergent vectors ' sprintf('%3.0fC - %3.0fC',TlT,TlB)]);
set(gca,'color','k');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
cmap = colormap('redblue');
cb = colorbar;
ylabel(cb,'Wm$^{-2}$');
ylim([-80 80]);

%%% Add sector values:
bcol = 0.8*[1 1 1];
lcol = 0.4*[1 1 1];%[0 0.5 0];
tcol = [0 0 0];
% Bering Strait:
lt = 66;ln1=-185;ln2=-160;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
BSflux = nansum(yfluxD(xind1:xind2,yind))/1e13
text(-155,65,sprintf('+%3.2fPW',BSflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
% Drake Passage:
ln = -70; lt1=-90;lt2=-50;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
DPflux = nansum(xfluxD(xind,1:yind))/1e13
text(-69,-75,sprintf('+%3.2fPW',DPflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
% South Atlantic:
lt = -30; ln1 = -60; ln2=+20;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
SAflux = nansum(yfluxD(xind1:xind2,yind))/1e13
text(-25,-25,sprintf('+%3.2fPW',SAflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
% Warm Route:
ln = 22; lt1=-90; lt2=-30;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
WRflux = nansum(xfluxD(xind,1:yind))/1e13
text(23,-55,sprintf('+%3.2fPW',WRflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
% 30N Pacific:
lt = 30; ln1 = -240; ln2=-108;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
P30flux = nansum(yfluxD(xind1:xind2,yind))/1e13
text(-170,35,sprintf('+%3.2fPW',P30flux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
% 30N Atlantic:
lt = 30; ln1 = -82.5; ln2=-5;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
A30flux = nansum(yfluxD(xind1:xind2,yind))/1e13
text(-40,35,sprintf('+%3.2fPW',A30flux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
% 30S Pacific:
lt = -30; ln1 = -215; ln2=-70;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Pm30flux = nansum(yfluxD(xind1:xind2,yind))/1e13
text(-170,-25,sprintf('-%3.2fPW',-Pm30flux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
% 30S Indian:
lt = -30; ln1 = 25; ln2=80;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Im30fluxP = nansum(yfluxD(xind1:xind2,yind))/1e13;
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
lt = -30; ln1 = -280; ln2=-240;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Im30fluxM = nansum(yfluxD(xind1:xind2,yind))/1e13;
Im30flux = Im30fluxP+Im30fluxM;
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(50,-25,sprintf('-%3.2fPW',-Im30flux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
% Aus South:
ln = -220; lt1=-90; lt2=-35;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
ASflux = nansum(xfluxD(xind,1:yind))/1e13
text(-219,-55,sprintf('-%3.2fPW',-ASflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
% ITF:
% ITF segments:
% 114.9E:
ln = -245.1; lt1 = -23; lt2 = -8.25;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFflux = nansum(xfluxD(t1,t2:t3))/1e13
plot([ln ln],[lt1 lt2],'--','color',lcol,'linewidth',2);
% 256.9E:
ln = -256.9; lt1 = -0.875; lt2 = 4.121;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFflux = ITFflux+nansum(xfluxD(t1,t2:t3))/1e13
plot([ln ln],[lt1 lt2],'-','color',lcol,'linewidth',2);
% 254.25E:
ln = -254.25; lt1 = -6.362; lt2 = -4.371;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFflux = ITFflux+nansum(xfluxD(t1,t2:t3))/1e13
plot([ln ln],[lt1 lt2],'-','color',lcol,'linewidth',2);
text(-244,-15,sprintf('-%3.2fPW',-ITFflux/100),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);


% Rotational and divergent parts:
xvec = 1:2:xL;
yvec = 1:2:yL;
figure;
set(gcf,'Position',[1923           5        1288         998]);
set(gcf,'defaulttextfontsize',10);
set(gcf,'defaultaxesfontsize',10);
dxu = ncread(gname,'dxu');
dyu = ncread(gname,'dyu');
flxes = {xflux./dyu,yflux./dxu,xfluxD./dyu,yfluxD./dxu,xflux./dyu-xfluxD./dyu,yflux./dxu-yfluxD./dxu};
tlts = {'Zonal flux','Meridional flux','Zonal Divergent flux', ...
        'Meridional Divergent flux','Zonal Rotational flux','Meridional Rotational flux'};
caxs = {[-10 10],[-10 10],[-0.25 0.25],[-0.25 0.25],[-10 10],[-10 10]};
for i=1:6
    subplot(3,2,i);
    flx = flxes{i};
    flx(flx==0) = NaN;
    pcolor(lon(xvec,yvec),lat(xvec,yvec),flx(xvec,yvec)/1e9);
    shading flat;
    hold on;
    caxis(caxs{i});
    title([tlts{i} ' (GW/m)']);
    set(gca,'color','k');
    xlabel('Longitude ($^\circ$E)');
    ylabel('Latitude ($^\circ$N)');
    colormap('redblue');
    colorbar;
    ylim([-70 75]);
end


% Production figure with Rotational and divergent parts:
xvec = 1:2:xL;
yvec = 1:2:yL;
poss = [0.0409    0.5145    0.4257    0.4271;
        0.4956    0.5145    0.4257    0.4271; ...
        0.0409    0.0545    0.4257    0.4271; ...
        0.4956    0.0545    0.4257    0.4271;];
figure;
set(gcf,'Position',[1921           1        1920        1005]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
i = 1;
ttls = {'Vertically-Integrated','20$^\circ$C - $34^\circ$C', ...
        '15$^\circ$C - $20^\circ$C','-2$^\circ$C - $15^\circ$C'};
subplot(2,2,i);
pcolor(lon(xvec,yvec),lat(xvec,yvec),divF(xvec,yvec)./area(xvec,yvec));
shading flat;
hold on;
minv = min(Pout(:));
maxv = max(Pout(:));
n = 30;
contour(lon,lat,Pout,[-1e50 minv:(maxv-minv)/n:maxv 1e50],'--k');
xvecV = 1:20:xL;
yvecV = 1:20:yL;
quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),xfluxD(xvecV,yvecV),yfluxD(xvecV,yvecV),2);
caxis([-300 300]);
set(gca,'color','k');
if (i>=3)
    xlabel('Longitude ($^\circ$E)');
end
if (i==1 | i == 3)
    ylabel('Latitude ($^\circ$N)');
end
if (i == 2 | i == 4)
    cb = colorbar;
    ylabel(cb,'Wm$^{-2}$');
end
colormap('redblue');
ylim([-75 80]);
text(-279,70,ttls{i},'Backgroundcolor','w','margin',0.1);
set(gca,'Position',poss(i,:));

% $$$ %%%% Old attempts:
% $$$ 
% $$$ % $$$ % Neive full matrix way:
% $$$ % $$$ one = ones(xL*yL,1);
% $$$ % $$$ D = spdiags([one -4*one one],[-1 0 1],xL*yL,xL*yL);
% $$$ % $$$ 
% $$$ % $$$ [II,JJ] = ind2sub([xL yL],1:(xL*yL));
% $$$ % $$$ ipjI = II+1;ipjI(II==xL) = 1;
% $$$ % $$$ imjI = II-1;imjI(II==1) = xL;
% $$$ % $$$ D(sub2ind([xL yL],ipjI,JJ))=1;
% $$$ % $$$ D(sub2ind([xL yL],imjI,JJ))=1;
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ jjpI = JJ+1;jjpI(JJ==yL) = 1;
% $$$ % $$$ jjmI = JJ-1;jjmI(JJ==1) = 1;
% $$$ % $$$ D(sub2ind([xL yL],II,))=1;
% $$$ % $$$ D(sub2ind([xL yL],II,))=1;
% $$$ % $$$ 
% $$$ % $$$ D = sub2ind([xL yL],II+1,JJ);
% $$$ % $$$ 
% $$$ % $$$ indsp1 = sub2ind([xL yL],1:xL,1:yL);
% $$$ % $$$ A=spdiags([],[],[],nn,nn,nn*5);
% $$$ 
% $$$ 
% $$$ % $$$ % More efficient attempt:
% $$$ % $$$ wetI = find(tWet); nn=length(wetI);
% $$$ % $$$ [IIinds,JJinds] = ind2sub([xL yL],wetI); % I and J indices of wet points
% $$$ % $$$ 
% $$$ % $$$ iINTinds = [1:nn];
% $$$ % $$$ ixLinds = find(IIinds==xL); % inds of i=xL points
% $$$ % $$$ i1inds = find(IIinds==1); % inds of i=1 points
% $$$ % $$$ jyLinds = find(JJinds==yL); % inds of i=xL points
% $$$ % $$$ j1inds = find(JJinds==1); % inds of i=1 points
% $$$ % $$$ iBDYinds = [ixLinds; i1inds; jyLinds; j1inds]; % inds of any
% $$$ % $$$                                                % boundary points
% $$$ % $$$ iINTinds(iBDYinds) = []; % inds of interior points
% $$$ % $$$ 
% $$$ % $$$ % Construct matrix:
% $$$ % $$$ A=sparse([],[],[],nn,nn,nn*5);
% $$$ % $$$ for ii=1:nn; A(ii,ii) = -4; end; % This step takes ages...
% $$$ % $$$ 
% $$$ % $$$ % Do interior points:
% $$$ % $$$ % ii+1, j:
% $$$ % $$$ iip1 = sub2ind([xL yL],IIinds(iINTinds)+1,JJinds(iINTinds));
% $$$ % $$$ iswet = ismember(iip1,wetI);
% $$$ % $$$ iswetinds = find(iswet);
% $$$ % $$$ iswetindsI = 0*iswetinds;
% $$$ % $$$ %for i=1:length(iswetinds); display(i); iswetindsI(i) = find(wetI==iip1(iswetinds(i))); ...
% $$$ % $$$ %end;
% $$$ % $$$ 
% $$$ % $$$ inds = sub2ind([nn nn],iINTinds(iswet),wetI(iip1(iswet)));
% $$$ % $$$ A(inds) = 1;
% $$$ % $$$ inds = sub2ind([nn nn],iINTinds(~iswet),iINTinds(~iswet));
% $$$ % $$$ A(inds) = A(inds)+1;
% $$$ % $$$ 
% $$$ % $$$ %
% $$$ % $$$ 
% $$$ % $$$ iim1 = sub2ind([xL yL],IIinds(iINTinds)-1,JJinds(iINTinds));
% $$$ % $$$ ijp1 = sub2ind([xL yL],IIinds(iINTinds),JJinds(iINTinds)+1);
% $$$ % $$$ ijm1 = sub2ind([xL yL],IIinds(iINTinds),JJinds(iINTinds)-1);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ %check for domain edge points where no exchange is possible:
% $$$ % $$$ tmp1=tWet; tmp1(:)=1; tmp2=exch_T_N(tmp1);
% $$$ % $$$ for iF=1:tWet.nFaces;
% $$$ % $$$ tmp3=tWet{iF}; tmp4=tmp2{iF}; 
% $$$ % $$$ tmp4=tmp4(2:end-1,1:end-2)+tmp4(2:end-1,3:end)+tmp4(1:end-2,2:end-1)+tmp4(3:end,2:end-1);
% $$$ % $$$ if ~isempty(find(isnan(tmp4)&~isnan(tmp3))); fprintf('warning: mask was modified\n'); end;
% $$$ % $$$ tmp3(isnan(tmp4))=NaN; tWet{iF}=tmp3; 
% $$$ % $$$ end;
% $$$ % $$$ %why does this matters... see set-ups with open boundary conditions...
% $$$ % $$$ 
% $$$ % $$$ %put 0 first guess if needed and switch land mask:
% $$$ % $$$ fld(find(isnan(fld)))=0; fld=fld.*tWet;
% $$$ 
% $$$ %define mapping from global array to (no nan points) global vector
% $$$ kk=find(tWet); nn=length(kk);
% $$$ KK=tWet; KK(:)=0; KK(kk)=kk; %global array indices
% $$$ LL=tWet; LL(:)=0; LL(kk)=[1:nn]; %global vector indices
% $$$ 
% $$$ %note:	diffsmooth2D_extrap_inv.m uses a Dirichlet boundary condition
% $$$ %       so that it uses a nan/1 mask in KK/LL // here we use a 
% $$$ %       Neumann boundary condition so we use a 0/1 mask (see below also)
% $$$ 
% $$$ %form matrix problem:
% $$$ A=sparse([],[],[],nn,nn,nn*5);
% $$$ for ii=1:3; for jj=1:3;
% $$$ 
% $$$ %1) seed points (FLDones) and neighborhood of influence (FLDkkFROM)
% $$$     FLDones=divF; FLDones(:)=0;
% $$$     FLDones(ii:3:end,jj:3:end)=1; 
% $$$     FLDones(KK==0)=0;
% $$$ 
% $$$     FLDkkFROMtmp=divF; FLDkkFROMtmp(:)=0;
% $$$     FLDkkFROMtmp(ii:3:end,jj:3:end)=KK(ii:3:end,jj:3:end);
% $$$     FLDkkFROMtmp(find(isnan(divF)))=0;
% $$$     
% $$$     FLDkkFROM=exch_T_N(FLDkkFROMtmp); FLDkkFROM(isnan(FLDkkFROM))=0; 
% $$$     tmp1=FLDkkFROM; tmp2=zeros(size(tmp1)-2);
% $$$     for ii2=1:3; for jj2=1:3; tmp2=tmp2+tmp1(ii2:end-3+ii2,jj2:end-3+jj2); end; end;
% $$$     FLDkkFROM=tmp2; 
% $$$ 
% $$$ %2) compute effect of each point on neighboring target point:
% $$$     [tmpU,tmpV]=calc_T_grad(FLDones,0);
% $$$ %unlike calc_T_grad, we work in grid point index, so we need to omit grid factors
% $$$     tmpU=tmpU.*mygrid.DXC; tmpV=tmpV.*mygrid.DYC;
% $$$ %and accordingly we use no grid scaling factor in div.
% $$$     [dFLDdt]=calc_UV_conv(tmpU,tmpV,{});
% $$$ 
% $$$ %note:	diffsmooth2D_extrap_inv.m uses a Dirichlet boundary condition
% $$$ %	so that it needs to apply mskFreeze to dFLDdt // here we use a 
% $$$ %	Neumann boundary condition so we do not mask dFLDdt (see above also)
% $$$ 
% $$$ %3) include seed contributions in matrix:
% $$$     FLDkkFROM=convert2array(FLDkkFROM);
% $$$ %3.1) for wet points
% $$$     dFLDdtWet=convert2array(dFLDdt.*tWet);
% $$$     tmp1=find(dFLDdtWet~=0&~isnan(dFLDdtWet)); 
% $$$     dFLDdtWet=dFLDdtWet(tmp1); FLDkkFROMtmp=FLDkkFROM(tmp1); FLDkkTOtmp=KK(tmp1);
% $$$     A=A+sparse(LL(FLDkkTOtmp),LL(FLDkkFROMtmp),dFLDdtWet,nn,nn);
% $$$ %3.2) for dry points (-> this part reflects the neumann boundary condition)
% $$$     dFLDdtDry=convert2array(dFLDdt.*mskDry);
% $$$     tmp1=find(dFLDdtDry~=0&~isnan(dFLDdtDry));
% $$$     dFLDdtDry=dFLDdtDry(tmp1); FLDkkFROMtmp=FLDkkFROM(tmp1);
% $$$     A=A+sparse(LL(FLDkkFROMtmp),LL(FLDkkFROMtmp),dFLDdtDry,nn,nn);
% $$$ 
% $$$ end; end;
% $$$ 
% $$$ %to check results: 
% $$$ %figure; spy(A);
% $$$ 
% $$$ %4) solve for potential:
% $$$ yy=convert2array(fld); yy=yy(find(KK~=0));
% $$$ xx=A\yy;
% $$$ yyFROMxx=A*xx; 
% $$$ 
% $$$ 
% $$$ 
% $$$ kk = find(tWet); nn = length(kk);
% $$$ KK = 0*tWet; KK(kk) = kk; 
% $$$ 
% $$$ 
% $$$ % $$$ dxt = ncread(gname,'dxt');dxt(isnan(dxt))=0;
% $$$ % $$$ dyt = ncread(gname,'dyt');dyt(isnan(dyt))=0;
% $$$ % $$$ dxu = ncread(gname,'dxu');dxu(isnan(dxu))=0;
% $$$ % $$$ dyu = ncread(gname,'dyu');dyu(isnan(dyu))=0;
% $$$ % $$$ 
% $$$ % $$$ % Derive some variables (approximate!):
% $$$ % $$$ dxte = cat(1,(dxt(end,:)+dxt(1,:))/2,avg(dxt,1));
% $$$ 
% $$$ % Make second derivative matrix:
% $$$ 
% $$$ % How to do boundary conditions? Force flux (i.e. all derivatives) to
% $$$ % be zero within all continents?
% $$$ 
% $$$ % dF/dx(i) = (F(i+1/2)-F(i-1/2))/dx(i)
% $$$ % d2F/dx2(i-1/2) = ((F(i+1/2)-F(i-1/2))/dx(i)-(F(i-1/2)-F(i-3/2))/dx(i-1))/dx(i-1/2)
% $$$ % d2F/dx2(i) = ((F(i+1)-F(i))/dx(i+1/2)-(F(i)-F(i-1))/dx(i-1/2))/dx(i)
% $$$ 
% $$$ % Multiply by mask?
% $$$ 
% $$$ %Construct Dirichlett Laplacian Matrix:
% $$$ one = ones(xL*yL,1);
% $$$ 
% $$$ coldx = spdiags([1/dxu
% $$$                 
% $$$ 
% $$$ %load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ 
% $$$ 
