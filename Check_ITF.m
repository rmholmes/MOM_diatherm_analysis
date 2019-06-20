
%Check ITF calculations

baseL = '/short/e14/rmh561/mom/archive/';
model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/'];
rstbaseD = baseD;
outD = [baseD 'mat_data/'];

output = 96;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output)];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
else
    fname = [base 'ocean.nc'];
end
gname = [base 'ocean_grid.nc'];
wname = [base 'ocean_wmass.nc'];

% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');

% Vertical grid  -----------------------------------------
z = ncread(fname,'st_ocean');zL = length(z);

% Time  -----------------------------------------
time = ncread(fname,'time');

tL = length(time);

% Temperature grid  -----------------------------------------
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

% Mask information:  -----------------------------------------
mask_t = ones(size(area)); %Pacific Mask: 1 = Pacific water, 0 = Elsewhere
mask_u = ones(size(lonu));

mask_u_Ny = 0*mask_u; %North y-trans mask
mask_u_Nx = 0*mask_u; %North x-trans mask
mask_u_Sy = 0*mask_u; %South y-trans mask
mask_u_Sx = 0*mask_u; %South x-trans mask
mask_u_Wy = 0*mask_u; %West y-trans mask
mask_u_Wx = 0*mask_u; %West x-trans mask

% ITF segments:
% 114.9E:
[~, t1] = min(abs(lonv_u+245.1));[~, t2] = min(abs(latv_u+23));[~, t3] = min(abs(latv_u+8.25));
mask_u_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;mask_u(1:t1,t2:t3) = 0;
% 256.9E:
[~, t1] = min(abs(lonv_u+256.9));[~, t2] = min(abs(latv_u+0.875));[~, t3] = min(abs(latv_u-4.121));
mask_u_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;mask_u(1:t1,t2:t3) = 0;
% 254.25E:
[~, t1] = min(abs(lonv_u+254.25));[~, t2] = min(abs(latv_u+6.362));[~, t3] = min(abs(latv_u+4.371));
mask_u_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;mask_u(1:t1,t2:t3) = 0;

% SOUTH Pacific:
% 46S:
[~, t1] = min(abs(latv_u+46));[~, t2] = min(abs(lonv_u+213.5));[~, t3] = min(abs(lonv_u+73));
mask_u_Sy(t2:t3,t1) = 1;mask_t(:,1:t1) = 0;mask_u(:,1:t1) = 0;

% 213.5E:
[~, t1] = min(abs(lonv_u+213.5));[~, t2] = min(abs(latv_u+46));[~, t3] = min(abs(latv_u+37));
mask_u_Sx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;mask_u(1:t1,t2:t3) = 0;

% NORTH Pacific:
[~, t1] = min(abs(latv_u-66));[~, t2] = min(abs(lonv_u+172));[~, t3] = min(abs(lonv_u+166));
mask_u_Ny(t2:t3,t1) = 1;
[~, t1a] = min(abs(latv_u-66.5));
mask_t(:,t1a:end) = 0;mask_u(:,t1a:end) = 0;mask_t(t2:end,t1:end) = 0;mask_u(t2:end,t1:end) = 0;

% $$$ % Pacific Mask:
% $$$ [~, t1] = min(abs(lonv_u+68));mask_t(t1:end,:) = 0;mask_u(t1:end,:) = 0;
% $$$ [~, t1] = min(abs(lonv_u+98));[~, t2] = min(abs(latv_u-18));mask_t(t1:end,t2:end) = 0;
% $$$ mask_u(t1:end,t2:end) = 0;[~, t1] = min(abs(lonv_u+89));[~, t2] = min(abs(latv_u-15));
% $$$ mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;[~, t1] = min(abs(lonv_u+84));
% $$$ [~, t2] = min(abs(latv_u-10));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
% $$$ [~, t1] = min(abs(lonv_u+82.5));[~, t2] = min(abs(latv_u-9.5));mask_t(t1:end,t2:end) = 0;
% $$$ mask_u(t1:end,t2:end) = 0;[~, t2] = min(abs(lonv_u+80.5));[~, t3] = min(abs(latv_u-9));
% $$$ mask_t(t1:t2,t3:end) = 0;mask_u(t1:t2,t3:end) = 0;[~, t1] = min(abs(lonv_u+78.25));
% $$$ [~, t2] = min(abs(latv_u-9));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
% $$$ [~, t1] = min(abs(lonv_u+77.5));[~, t2] = min(abs(latv_u-7.25));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
% $$$ [~, t1] = min(abs(lonv_u+217));[~, t2] = min(abs(latv_u+22.29));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ [~, t1] = min(abs(lonv_u+261));mask_t(1:t1,:) = 0;mask_u(1:t1,:) = 0;
% $$$ [~, t1] = min(abs(lonv_u+256.75));[~, t2] = min(abs(latv_u+1));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ [~, t1] = min(abs(lonv_u+252));[~, t2] = min(abs(latv_u+6.5));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ [~, t1] = min(abs(lonv_u+248));[~, t2] = min(abs(latv_u+7.25));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ [~, t1] = min(abs(lonv_u+260));[~, t2] = min(abs(latv_u-8));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ [~, t1] = min(abs(lonv_u+258.5));[~, t2] = min(abs(latv_u-6.5));mask_t(1:t1,1:t2) = 0;mask_u(1:t1,1:t2) = 0;
% $$$ mask_t(81,525) = 0;mask_t(799,535)=0;
% $$$ 
% $$$ txtrans = sum(ncread(wname,'tx_trans_nrho',[1 1 1 1],[xL yL TL 1]),3);
SST = ncread(fname,'temp',[1 1 1 1],[xL yL 1 1]);
mask_t(isnan(SST)) = 0;
mask_u(txtrans==0) = 0;

outputs = 96;

ln1 = 50;
ln2 = 150;
lt1 = 400;
lt2 = 550;



mask = mask_u_Wx(ln1:ln2,lt1:lt2);
lonus = lonu(ln1:ln2,lt1:lt2);
latus = latu(ln1:ln2,lt1:lt2);
lonusm = lonus(find(mask));
latusm = latus(find(mask));
[xL,yL] = size(mask);
txtrans = zeros(xL,yL,TL);
hxtrans = zeros(xL,yL,TL);

for output = outputs
    base = [baseD sprintf('output%03d/',output)];
    wname = sprintf([base 'ocean_wmass.nc'],output)
    
    for i=1:12
        i
        txtrans = txtrans+ncread(wname,'tx_trans_nrho',[ln1 lt1 1 i],[ln2-ln1+1 ...
                            lt2-lt1+1 TL 1])*1e9/rho0;
        hxtrans = hxtrans+ncread(wname,'temp_xflux_adv_on_nrho',[ln1 lt1 1 i],[ln2-ln1+1 ...
                            lt2-lt1+1 TL 1]);
    end
end
txtrans = txtrans/12;
hxtrans = hxtrans/12;

txtrans_sum = zeros(TL,1);
hxtrans_sum = zeros(TL,1);

for i=1:TL
    tx = txtrans(:,:,i);
    txtrans_sum(i) = nansum(tx(find(mask==1)));
    hx = hxtrans(:,:,i);
    hxtrans_sum(i) = nansum(hx(find(mask==1)));
end

PSI = cat(1,0,cumsum(txtrans_sum));
A = cat(1,0,cumsum(hxtrans_sum));
AE = rho0*Cp*PSI.*Te;
AI = A-AE;

figure;
subplot(1,2,1);
plot(PSI/1e6,Te,'-k');
xlabel('ITF $\Psi(T)$ (Sv)');
ylabel('Temperature ($^\circ$C)');
grid on;

subplot(1,2,2);
plot(A/1e15,Te,'-k');
hold on;
plot(AE/1e15,Te,'-b');
plot(AI/1e15,Te,'-r');
xlabel('ITF heat transport (PW)');
legend('A','AE','AI');
grid on;

