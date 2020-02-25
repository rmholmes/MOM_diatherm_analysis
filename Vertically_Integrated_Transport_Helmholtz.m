% This script analyses the vertically-integrated volume and heat
% transports, also performing a Helmholtz decomposition on the lateral
% heat fluxes following Forget and Ferriera (2019) Nat Geo.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
gname = '/srv/ccrc/data03/z3500785/mom/ocean_grid_mom025.nc';

RUNS = { ...
% $$$     {'MOM025_kb3seg',[101110]}, ...
    {'ACCESS-OM2_1deg_ryf',[31]}, ...
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
qxflux(qxflux==0) = NaN;
qyflux(qyflux==0) = NaN;
xDry = isnan(qxflux);
yDry = isnan(qyflux);
xWet = ~isnan(qxflux);
yWet = ~isnan(qyflux);

% Load lateral fluxes:
TlT = 34;
TlB = -10;
if (TlB<-5)
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlT),'.','p') 'C.mat']);
else
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlT),'.','p') 'C.mat']);
    qxfluxT = qxflux;qyfluxT=qyflux;
    load([base model sprintf('_output%03d_XYtrans_T',outputs(1)) ...
          strrep(num2str(TlB),'.','p') 'C.mat']);
    qxflux = qxfluxT-qxflux;
    qyflux = qyfluxT-qyflux;
end
    
qxflux(xDry)=0;
qyflux(yDry)=0;
mxflux(xDry)=0;
myflux(yDry)=0;

% Calculate divergence:
qdivF = zeros(xL,yL);
qdivF(2:end,2:end) = qdivF(2:end,2:end) + (qxflux(1:(end-1),2:end) - qxflux(2:end,2:end) ...
                                     +qyflux(2:end,1:(end-1)) - qyflux(2:end,2:end));%./area(2:end,2:end);
qdivF(1,2:end) = qdivF(1,2:end) + (qxflux(end,2:end) - qxflux(1,2:end) ...
                             +qyflux(1,1:(end-1)) - qyflux(1,2:end));%./area(1,2:end);        
qdivF(tDry) = NaN;

qdivF(:,end) = 0;

% Calculate divergence:
mdivF = zeros(xL,yL);
mdivF(2:end,2:end) = mdivF(2:end,2:end) + (mxflux(1:(end-1),2:end) - mxflux(2:end,2:end) ...
                                     +myflux(2:end,1:(end-1)) - myflux(2:end,2:end));%./area(2:end,2:end);
mdivF(1,2:end) = mdivF(1,2:end) + (mxflux(end,2:end) - mxflux(1,2:end) ...
                             +myflux(1,1:(end-1)) - myflux(1,2:end));%./area(1,2:end);        
mdivF(tDry) = NaN;

mdivF(:,end) = 0;

wetI = find(tWet); nn=length(wetI);

A=sparse([],[],[],nn,nn,nn*5);

% Less efficient option, naive tripolar points:
for ll=1:nn % result index loop

    if (mod(ll,100)==0); display([num2str(ll) ' of ' num2str(nn)]); end;
    
    [ii,jj] = ind2sub([xL yL],wetI(ll)); % Get full i,j indices
    
    A(ll,ll) = -4;
    
    % Build up matrix coefficients:
    % ii+1, jj:
    if (ii==xL); lind = sub2ind([xL yL],1,jj); else; lind = sub2ind([xL yL],ii+1,jj); end;
    
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end
    

    % ii-1, jj:
    if (ii==1); lind = sub2ind([xL yL],xL,jj); else; lind = sub2ind([xL yL],ii-1,jj); end;
    
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end
    

    % ii, jj+1:
    if (jj==yL)
        A(ll,ll) = A(ll,ll)+1;
    else
        lind = sub2ind([xL yL],ii,jj+1);
        ind = find(wetI==lind);
        if (length(ind)==1) % source wet point -> Add dependence
            A(ll,ind) = 1; 
        else
        A(ll,ll) = A(ll,ll)+1;
        end
    end
    

    % ii, jj-1:
    lind = sub2ind([xL yL],ii,jj-1);
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end

end
save('HelmholtzA_AOM2_1deg.mat','A');

load('HelmholtzA_AOM2_1deg.mat','A');
Aold = A;

% Fix tri-polar jj=nj points:
[IIinds,JJinds] = ind2sub([xL yL],wetI); % I and J indices of wet points

inds = find(JJinds==yL);
for lli=1:length(inds)
    ll = inds(lli)
    [ii,jj] = ind2sub([xL yL],wetI(ll)); % Get full i,j indices    

    A(ll,:) = 0*A(ll,:);
    A(ll,ll) = -4;
    
    % Build up matrix coefficients:
    % ii+1, jj:
    if (ii==xL); lind = sub2ind([xL yL],1,jj); else; lind = sub2ind([xL yL],ii+1,jj); end;
    
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end
    
    % ii-1, jj:
    if (ii==1); lind = sub2ind([xL yL],xL,jj); else; lind = sub2ind([xL yL],ii-1,jj); end;
    
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end
    
    % ii, jj+1:
    if (jj~=yL); display('ERROR!'); return; end;

    lind = sub2ind([xL yL],xL-ii+1,jj); % T_{i,yL+1} = T_{xL-i+1,yL}
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end

    % ii, jj-1:
    lind = sub2ind([xL yL],ii,jj-1);
    ind = find(wetI==lind);
    if (length(ind)==1) % source wet point -> Add dependence
        A(ll,ind) = 1; 
    else
        A(ll,ll) = A(ll,ll)+1;
    end
end
save('HelmholtzA_tripolarcorrected_AOM2_1deg.mat','A');
load('HelmholtzA_tripolarcorrected_AOM2_1deg.mat','A');

% Solve it:
b = qdivF(tWet);
P = -A\b;
bm = mdivF(tWet);
Pm = -A\bm;

% unwrap:
Pout = NaN*zeros(xL,yL);
Pout(tWet) = P-mean(P(:));

qxfluxD = zeros(size(Pout));
qxfluxD(1:(end-1),:) = Pout(2:end,:)-Pout(1:(end-1),:);
qxfluxD(end,:) = Pout(1,:)-Pout(end,:);

qyfluxD = zeros(size(Pout));
qyfluxD(:,1:(end-1)) = Pout(:,2:end)-Pout(:,1:(end-1));

% unwrap mass:
Pmout = NaN*zeros(xL,yL);
Pmout(tWet) = Pm-mean(Pm(:));

mxfluxD = zeros(size(Pmout));
mxfluxD(1:(end-1),:) = Pmout(2:end,:)-Pmout(1:(end-1),:);
mxfluxD(end,:) = Pmout(1,:)-Pmout(end,:);

myfluxD = zeros(size(Pmout));
myfluxD(:,1:(end-1)) = Pmout(:,2:end)-Pmout(:,1:(end-1));

% $$$ MHTD = nansum(qyfluxD,1);
% $$$ MHT = nansum(qyflux,1);
% $$$ 
% $$$ % MHT:
% $$$ figure;
% $$$ plot(yu,MHT/1e15);
% $$$ hold on;
% $$$ plot(yu,MHTD/1e15,'--r');

qdivFch = zeros(xL,yL);
qdivFch(2:end,2:end) = qdivFch(2:end,2:end) + (qxfluxD(1:(end-1),2:end) - qxfluxD(2:end,2:end) ...
                                     +qyfluxD(2:end,1:(end-1)) - qyfluxD(2:end,2:end));%./area(2:end,2:end);
qdivFch(1,2:end) = qdivFch(1,2:end) + (qxfluxD(end,2:end) - qxfluxD(1,2:end) ...
                             +qyfluxD(1,1:(end-1)) - qyfluxD(1,2:end));%./area(1,2:end);        
qdivFch(tDry) = NaN;


mdivFch = zeros(xL,yL);
mdivFch(2:end,2:end) = mdivFch(2:end,2:end) + (mxfluxD(1:(end-1),2:end) - mxfluxD(2:end,2:end) ...
                                     +myfluxD(2:end,1:(end-1)) - myfluxD(2:end,2:end));%./area(2:end,2:end);
mdivFch(1,2:end) = mdivFch(1,2:end) + (mxfluxD(end,2:end) - mxfluxD(1,2:end) ...
                             +myfluxD(1,1:(end-1)) - myfluxD(1,2:end));%./area(1,2:end);        
mdivFch(tDry) = NaN;

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
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),qdivF(xvec,yvec)./area(xvec,yvec));
% $$$ hold on;
% $$$ xvecV = 1:10:xL;
% $$$ yvecV = 1:10:yL;
% $$$ quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),qxflux(xvecV,yvecV),qyflux(xvecV,yvecV),8);
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
contourf(lon(xvec,yvec),lat(xvec,yvec),qdivFch(xvec,yvec)./ ...
         area(xvec,yvec),cpts,'linestyle','none');
% $$$ contourf(lon(xvec,yvec),lat(xvec,yvec),qdivF(xvec,yvec)./ ...
% $$$          area(xvec,yvec),cpts,'linestyle','none');
hold on;
minv = min(Pout(:));
maxv = max(Pout(:));
n = 40;
contour(lon,lat,Pout,[-1e50 minv:0.1e15:maxv 1e50],'-k');
xvecV = 1:5:xL;yvecV = 1:5:yL;
% $$$ xvecV = 1:20:xL;yvecV = 1:20:yL;
quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),qxfluxD(xvecV,yvecV),qyfluxD(xvecV,yvecV),2);
%quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),qxflux(xvecV,yvecV),qyflux(xvecV,yvecV),2);
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

%%% Add sector values (AOM2 1-degree):
bcol = 0.8*[1 1 1];
lcol = 0.4*[1 1 1];%[0 0.5 0];
tcol = [0 0 0];
tcol1 = 'r';tcol2 = 'b';tcol3 = 'm';
Tref2 = 20;
% Bering Strait:
lt = 66;ln1=-185;ln2=-160;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
BSfluxD = nansum(qyfluxD(xind1:xind2,yind));
BSflux = nansum(qyflux(xind1:xind2,yind));
BSfluxmD = nansum(myfluxD(xind1:xind2,yind));
BSfluxm = nansum(myflux(xind1:xind2,yind));
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(-155-21,65,sprintf('%3.2fSv',BSfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-155-21,65-5,sprintf('%3.2fSv',BSfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-155,65,sprintf('%3.2fPW',BSflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-155,65-5,sprintf('%3.2fPW',BSfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-155+21,65,sprintf('%3.2fPW',(BSflux-rho0*Cp*Tref2*BSfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-155+21,65-5,sprintf('%3.2fPW',(BSfluxD-rho0*Cp*Tref2*BSfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% Drake Passage:
ln = -70; lt1=-90;lt2=-50;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
DPfluxD = nansum(qxfluxD(xind,1:yind));
DPflux = nansum(qxflux(xind,1:yind));
DPfluxmD = nansum(mxfluxD(xind,1:yind));
DPfluxm = nansum(mxflux(xind,1:yind));
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
text(-69-21,-75,sprintf('%3.2fSv',DPfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-69-21,-75-5,sprintf('%3.2fSv',DPfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-69,-75,sprintf('%3.2fPW',DPflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-69,-75-5,sprintf('%3.2fPW',DPfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-69+21,-75,sprintf('%3.2fPW',(DPflux-rho0*Cp*Tref2*DPfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-69+21,-75-5,sprintf('%3.2fPW',(DPfluxD-rho0*Cp*Tref2*DPfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% South Atlantic:
lt = -30; ln1 = -60; ln2=+20;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
SAfluxD = nansum(qyfluxD(xind1:xind2,yind));
SAflux = nansum(qyflux(xind1:xind2,yind));
SAfluxmD = nansum(myfluxD(xind1:xind2,yind));
SAfluxm = nansum(myflux(xind1:xind2,yind));
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(-25-21,-25,sprintf('%3.2fSv',SAfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-25-21,-25-5,sprintf('%3.2fSv',SAfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-25,-25,sprintf('%3.2fPW',SAflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-25,-25-5,sprintf('%3.2fPW',SAfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-25+21,-25,sprintf('%3.2fPW',(SAflux-rho0*Cp*Tref2*SAfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-25+21,-25-5,sprintf('%3.2fPW',(SAfluxD-rho0*Cp*Tref2*SAfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% Warm Route:
ln = 22; lt1=-90; lt2=-30;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
WRfluxD = nansum(qxfluxD(xind,1:yind));
WRflux = nansum(qxflux(xind,1:yind));
WRfluxmD = nansum(mxfluxD(xind,1:yind));
WRfluxm = nansum(mxflux(xind,1:yind));
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
text(23-21,-55,sprintf('%3.2fSv',WRfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(23-21,-55-5,sprintf('%3.2fSv',WRfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(23,-55,sprintf('%3.2fPW',WRflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(23,-55-5,sprintf('%3.2fPW',WRfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(23+21,-55,sprintf('%3.2fPW',(WRflux-rho0*Cp*Tref2*WRfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(23+21,-55-5,sprintf('%3.2fPW',(WRfluxD-rho0*Cp*Tref2*WRfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% 30N Pacific:
lt = 30; ln1 = -240; ln2=-108;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
P30fluxD = nansum(qyfluxD(xind1:xind2,yind));
P30flux = nansum(qyflux(xind1:xind2,yind));
P30fluxmD = nansum(myfluxD(xind1:xind2,yind));
P30fluxm = nansum(myflux(xind1:xind2,yind));
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(-170-21,35,sprintf('%3.2fSv',P30fluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-170-21,35-5,sprintf('%3.2fSv',P30fluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-170,35,sprintf('%3.2fPW',P30flux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-170,35-5,sprintf('%3.2fPW',P30fluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-170+21,35,sprintf('%3.2fPW',(P30flux-rho0*Cp*Tref2*P30fluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-170+21,35-5,sprintf('%3.2fPW',(P30fluxD-rho0*Cp*Tref2*P30fluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% 30N Atlantic:
lt = 30; ln1 = -82.5; ln2=-5;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
A30fluxD = nansum(qyfluxD(xind1:xind2,yind));
A30flux = nansum(qyflux(xind1:xind2,yind));
A30fluxmD = nansum(myfluxD(xind1:xind2,yind));
A30fluxm = nansum(myflux(xind1:xind2,yind));
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(-40-21,35,sprintf('%3.2fSv',A30fluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-40-21,35-5,sprintf('%3.2fSv',A30fluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-40,35,sprintf('%3.2fPW',A30flux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-40,35-5,sprintf('%3.2fPW',A30fluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-40+21,35,sprintf('%3.2fPW',(A30flux-rho0*Cp*Tref2*A30fluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-40+21,35-5,sprintf('%3.2fPW',(A30fluxD-rho0*Cp*Tref2*A30fluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% 30S Pacific:
lt = -30; ln1 = -215; ln2=-70;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Pm30fluxD = nansum(qyfluxD(xind1:xind2,yind));
Pm30flux = nansum(qyflux(xind1:xind2,yind));
Pm30fluxmD = nansum(myfluxD(xind1:xind2,yind));
Pm30fluxm = nansum(myflux(xind1:xind2,yind));
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(-170-21,-25,sprintf('%3.2fSv',Pm30fluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-170-21,-25-5,sprintf('%3.2fSv',Pm30fluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-170,-25,sprintf('%3.2fPW',Pm30flux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-170,-25-5,sprintf('%3.2fPW',Pm30fluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-170+21,-25,sprintf('%3.2fPW',(Pm30flux-rho0*Cp*Tref2*Pm30fluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-170+21,-25-5,sprintf('%3.2fPW',(Pm30fluxD-rho0*Cp*Tref2*Pm30fluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% 30S Indian:
lt = -30; ln1 = 25; ln2=80;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Im30fluxPD = nansum(qyfluxD(xind1:xind2,yind));;
Im30fluxP = nansum(qyflux(xind1:xind2,yind));;
Im30fluxPmD = nansum(myfluxD(xind1:xind2,yind));;
Im30fluxPm = nansum(myflux(xind1:xind2,yind));;
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
lt = -30; ln1 = -280; ln2=-240;
[tmp yind] = min(abs(yu-lt));[tmp xind1] = min(abs(xu-ln1));[tmp xind2] = min(abs(xu-ln2));
Im30fluxMD = nansum(qyfluxD(xind1:xind2,yind));;
Im30fluxM = nansum(qyflux(xind1:xind2,yind));;
Im30fluxMmD = nansum(myfluxD(xind1:xind2,yind));;
Im30fluxMm = nansum(myflux(xind1:xind2,yind));;
Im30fluxD = Im30fluxPD+Im30fluxMD;
Im30flux = Im30fluxP+Im30fluxM;
Im30fluxmD = Im30fluxPmD+Im30fluxMmD;
Im30fluxm = Im30fluxPm+Im30fluxMm;
plot([ln1 ln2],[lt lt],'--','color',lcol,'linewidth',2);
text(45-21,-25,sprintf('%3.2fSv',Im30fluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(45-21,-25-5,sprintf('%3.2fSv',Im30fluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(45,-25,sprintf('%3.2fPW',Im30flux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(45,-25-5,sprintf('%3.2fPW',Im30fluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(45+21,-25,sprintf('%3.2fPW',(Im30flux-rho0*Cp*Tref2*Im30fluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(45+21,-25-5,sprintf('%3.2fPW',(Im30fluxD-rho0*Cp*Tref2*Im30fluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% Aus South:
ln = -220; lt1=-90; lt2=-35;
[tmp xind] = min(abs(xu-ln));[tmp yind] = min(abs(yu-lt2));
ASfluxD = nansum(qxfluxD(xind,1:yind));
ASflux = nansum(qxflux(xind,1:yind));
ASfluxmD = nansum(mxfluxD(xind,1:yind));
ASfluxm = nansum(mxflux(xind,1:yind));
plot([ln ln],[-75 lt2],'--','color',lcol,'linewidth',2);
text(-219-21,-55,sprintf('%3.2fSv',ASfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-219-21,-55-5,sprintf('%3.2fSv',ASfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-219,-55,sprintf('%3.2fPW',ASflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-219,-55-5,sprintf('%3.2fPW',ASfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-219+21,-55,sprintf('%3.2fPW',(ASflux-rho0*Cp*Tref2*ASfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-219+21,-55-5,sprintf('%3.2fPW',(ASfluxD-rho0*Cp*Tref2*ASfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);
% ITF:
% ITF segments:
% 114.9E:
ln = -245.1; lt1 = -23; lt2 = -8.25;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFfluxD = nansum(qxfluxD(t1,t2:t3));
ITFflux = nansum(qxflux(t1,t2:t3));
ITFfluxmD = nansum(mxfluxD(t1,t2:t3));
ITFfluxm = nansum(mxflux(t1,t2:t3));
plot([ln ln],[lt1 lt2],'--','color',lcol,'linewidth',2);
% 256.9E:
ln = -256.9; lt1 = -0.875; lt2 = 4.121;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFfluxD = ITFfluxD+nansum(qxfluxD(t1,t2:t3));
ITFflux = ITFflux+nansum(qxflux(t1,t2:t3));
ITFfluxmD = ITFfluxmD+nansum(mxfluxD(t1,t2:t3));
ITFfluxm = ITFfluxm+nansum(mxflux(t1,t2:t3));
plot([ln ln],[lt1 lt2],'-','color',lcol,'linewidth',2);
% 254.25E:
ln = -254.25; lt1 = -6.362; lt2 = -4.371;
[~, t1] = min(abs(xu-ln));[~, t2] = min(abs(yu-lt1));[~, t3] = min(abs(yu-lt2));
ITFfluxD = ITFfluxD+nansum(qxfluxD(t1,t2:t3));
ITFflux = ITFflux+nansum(qxflux(t1,t2:t3));
ITFfluxmD = ITFfluxmD+nansum(mxfluxD(t1,t2:t3));
ITFfluxm = ITFfluxm+nansum(mxflux(t1,t2:t3));
plot([ln ln],[lt1 lt2],'-','color',lcol,'linewidth',2);
text(-244-21,-15,sprintf('%3.2fSv',ITFfluxm/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-244-21,-15-5,sprintf('%3.2fSv',ITFfluxmD/1e6),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-244,-15,sprintf('%3.2fPW',ITFflux/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-244,-15-5,sprintf('%3.2fPW',ITFfluxD/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-244+21,-15,sprintf('%3.2fPW',(ITFflux-rho0*Cp*Tref2*ITFfluxm)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-244+21,-15-5,sprintf('%3.2fPW',(ITFfluxD-rho0*Cp*Tref2*ITFfluxmD)/1e15),'Backgroundcolor',bcol,'margin',0.01,'color',tcol3);

text(-255-21,50,'Vol','Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-255-21,50-5,'Div. Vol','Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-255,50,'Heat 0C','Backgroundcolor',bcol,'margin',0.01,'color',tcol);
text(-255,50-5,'Div. Heat 0C','Backgroundcolor',bcol,'margin',0.01,'color',tcol1);
text(-255+30,50,'Heat 20C','Backgroundcolor',bcol,'margin',0.01,'color',tcol2);
text(-255+30,50-5,'Div. Heat 20C','Backgroundcolor',bcol,'margin',0.01,'color',tcol3);

title(['ACCESS-OM2 1-degree 10yr-average heat flux divergence ' ...
       '(color), vector potential (contour) and divergent heat ' ...
      'flux vectors']);


% Plot section transports as a function of the reference
% temperature:
Tref = -2:1:34;
nT = length(Tref);
secs = {'BS','DP','SA','WR','P30','A30','Pm30','Im30','AS','ITF'};
nsec = length(secs);

full = zeros(nsec,nT);
div = full;
rot = full;

for si = 1:nsec
    eval(['qflx = ' secs{si} 'flux;']);
    eval(['mflx = ' secs{si} 'fluxm;']);
% $$$     eval(['qflxD = ' secs{si} 'fluxD;']);
% $$$     eval(['mflxD = ' secs{si} 'fluxmD;']);

    full(si,:) = qflx - rho0*Cp*Tref*mflx;
% $$$     div(si,:) = qflxD + rho0*Cp*Tref*mflxD;
% $$$     rot(si,:) = full(si,:) - div(si,:);
end

cols = {'k','r','b','m',[0.5 0.5 0.5],'k','r','b','m',[0.5 0.5 0.5]};
stys = {'-','-','-','-','-','--','--','--','--','--'};
% Plot:
figure;
subplot(2,1,1);
for si=1:nsec
    plot(Tref,full(si,:)/1e15,stys{si},'color',cols{si});
    hold on;
end
legend(secs);
xlabel('Reference Temperature ($^\circ$C)');
ylabel('Heat Transport (PW)');




% Rotational and divergent parts:
xvec = 1:2:xL;
yvec = 1:2:yL;
figure;
set(gcf,'Position',[1923           5        1288         998]);
set(gcf,'defaulttextfontsize',10);
set(gcf,'defaultaxesfontsize',10);
dxu = ncread(gname,'dxu');
dyu = ncread(gname,'dyu');
flxes = {qxflux./dyu,qyflux./dxu,qxfluxD./dyu,qyfluxD./dxu,qxflux./dyu-qxfluxD./dyu,qyflux./dxu-qyfluxD./dxu};
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
pcolor(lon(xvec,yvec),lat(xvec,yvec),qdivF(xvec,yvec)./area(xvec,yvec));
shading flat;
hold on;
minv = min(Pout(:));
maxv = max(Pout(:));
n = 30;
contour(lon,lat,Pout,[-1e50 minv:(maxv-minv)/n:maxv 1e50],'--k');
xvecV = 1:20:xL;
yvecV = 1:20:yL;
quiver(lon(xvecV,yvecV),lat(xvecV,yvecV),qxfluxD(xvecV,yvecV),qyfluxD(xvecV,yvecV),2);
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
% $$$     FLDones=qdivF; FLDones(:)=0;
% $$$     FLDones(ii:3:end,jj:3:end)=1; 
% $$$     FLDones(KK==0)=0;
% $$$ 
% $$$     FLDkkFROMtmp=qdivF; FLDkkFROMtmp(:)=0;
% $$$     FLDkkFROMtmp(ii:3:end,jj:3:end)=KK(ii:3:end,jj:3:end);
% $$$     FLDkkFROMtmp(find(isnan(qdivF)))=0;
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
% $$$ %and accordingly we use no grid scaling factor in qdiv.
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
