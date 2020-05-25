% This script calculates global average heat content as a function
% of depth from ACCESS-OM2 runs.

baseL = '/g/data/e14/mv7494/access-om2/archive/';

% MOM-SIS025:
model = 'ACCESS-OM2_025deg_jra55_iaf';
baseD = [baseL '025deg_jra55_iaf/']; %Data Directory.
ICdir = [baseL '025deg_jra55_iaf/'];

% $$$ outD = [baseD 'mat_data/'];
outD = ['/g/data/e14/rmh561/access-om2/archive/025deg_jra55_iaf/mat_data/'];
rstbaseD = baseD;

post = 'ocean/'; % For ACCESS-OM2 output coulpled;
% $$$ post = ''; % For MOM-SIS.

restart = output-1;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];
if (output==0)
    baser = ICdir;
else
    baser = [rstbaseD sprintf('restart%03d/',restart) post];
end
fname = [base 'ocean.nc'];
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
tname = [base 'time_stamp.out'];
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
xt = ncread(gname,'xt_ocean');xu = ncread(gname,'xu_ocean');
yt = ncread(gname,'yt_ocean');yu = ncread(gname,'yu_ocean');

% Vertical grid  -----------------------------------------
z = ncread(fname,'st_ocean');zL = length(z);
if (output ==0)
    % Initial dzt accounting for partial bottom cells:
    ht = ncread(gname,'ht');ze = ncread(fname,'st_edges_ocean');
    kmt = ncread(gname,'kmt');
    dztI = repmat(permute(diff(ze),[3 2 1]),[size(ht) 1]);
    for ii = 1:xL
        for jj = 1:yL
            if (kmt(ii,jj)>1)
                dztI(ii,jj,kmt(ii,jj)) = ht(ii,jj) - ze(kmt(ii,jj));
            end
        end
    end    
end

% 3D mask ------------------------------------------------
mask = ncread(fname,'temp',[1 1 1 1],[xL yL zL 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

if (~strcmp(region,'Global'))
    [tmaskREG,umaskREG,~] = Heat_Budget_Mask(region,gname,wname,outD,model);
else
    tmaskREG = ones(xL,yL);
    umaskREG = ones(xL,yL);
end

% Time  -----------------------------------------
time = ncread(fname,'time');
ndays = ncread(fname,'average_DT');
tL = length(time);
tLoc = length(ncread(fname,'time'));

rnameT = [baser 'ocean_temp_salt.res.nc'];
rnameZ = [baser 'ocean_thickness.res.nc'];
if (exist(rnameT) & exist(rnameZ))
    found_rst = 1;rstti = 1;
else
    found_rst = 0;rstti = tL;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
end

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

V = zeros(zL,tL);
H = zeros(zL,tL);

%Do other times for Vsnap and Hsnap:
for ti=1:tL
    for zi=1:zL
        sprintf('Calculating V and H months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
        
        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        temp(~mask(:,:,zi)) = NaN;
        if (max(max(temp))>120);temp = temp-273.15;end;
        Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Vol(isnan(Vol)) = 0;

        V(zi,ti) = nansum(nansum(Vol));
        H(zi,ti) = rho0*Cp*nansum(nansum(temp.*Vol));
    end
end

save([outD model sprintf('_output%03d',output) '_VHofz.mat'],'V','H','z','time','-v7.3');
