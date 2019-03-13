% This script calculates vertically-integrated and time-averaged
% lateral transports from the MOM025 control simulation

% $$$ baseL = '/short/e14/rmh561/mom/archive/';
baseL = '/g/data/e14/rmh561/mom/';
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/srv/ccrc/data03/z3500785/';

% MOM-SIS025-WOMBAT:
% $$$ model = 'MOM025';
% $$$ baseD = [baseL 'MOM_wombat/']; %Data Directory.
% $$$ % MOM-SIS025:
model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
% ACCESS-OM2:
% $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ baseD = [baseL '1deg_jra55_ryf8485_kds50_may/']; %Data Directory.
% $$$ ICdir = '/g/data1/ua8/MOM/initial_conditions/WOA/10_KDS50/';
% MOM-SIS01:
% $$$ model = 'MOM01';
% $$$ baseD = [baseL 'MOM01_HeatDiag/']; %Data Directory.

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

% $$$ post = 'ocean/'; % For ACCESS-OM2 output coulpled;
post = ''; % For MOM-SIS.

haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
haveGM = 0; % 1 = GM is on, 0 = off;
haveSUB = 1; % 1 = submeso is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;
haveMIX = 1; % 1 = Do mixing components (vdiffuse_diff_cbt_*), 0 = don't. 
haveHND = 1; % 1 = Do numerical mixing via heat budget.

% scaling constant on the transports:
if (strcmp(model(1),'A')) %ACCESS-OM2, transport in kg/s
    tsc = 1;
else % MOM-SIS, transport in 1e9 kg/s
    tsc = 1e9;
end

for output = 86:90;
% $$$ output=91;
restart = output-1;

region = 'Global';

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];
if (output==0)
    baser = ICdir;
else
    baser = [rstbaseD sprintf('restart%03d/',restart) post];
end
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
else
    fname = [base 'ocean.nc'];
end
fname2 = [base 'ocean_month.nc'];
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
tname = [base 'time_stamp.out'];
if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
else
    found_rst = 0;rstti = 12;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
end

% Time  -----------------------------------------
time = ncread(fname,'time');

dys = [31 28 31 30 31 30 31 31 30 31 30 31];
C = textread(tname, '%s','delimiter', '\n');
C = strsplit(C{1});
rtime = [str2num(C{1}) str2num(C{2}) str2num(C{3}) str2num(C{4}) str2num(C{5}) str2num(C{6})];
time_snap = [(rtime(1)-1)*365+sum(dys(1:(rtime(2)-1)))+(rtime(3)-1)+rtime(4)/24+rtime(5)/24/60+rtime(6)/24/60/60;
             ncread(sname,'time')];

tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

for ti=1:tL
    for Ti=1:TL
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        if (haveSUB)
            txtrans = txtrans + ncread(wname,'tx_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        end
        if (haveGM)
            txtrans = txtrans + ncread(wname,'tx_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        end
        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveSUB)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            qytrans = qytrans + ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end 
        if (haveGM)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
    end
end


save([outD model sprintf('_output%03d',output) '_LatTrans.mat'],'qx','qy','tx','ty','rho0','Cp','Tmax','-v7.3');
