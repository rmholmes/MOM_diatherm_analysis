% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseL = '/short/e14/rmh561/mom/archive/';
% $$$ baseL = '/g/data/e14/rmh561/';
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/srv/ccrc/data03/z3500785/';

% MOM-SIS025-WOMBAT:
% $$$ model = 'MOM025';
% $$$ baseD = [baseL 'MOM_wombat/']; %Data Directory.
% $$$ % MOM-SIS025:
% $$$ model = 'MOM025_kb3seg';
% $$$ baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
% ACCESS-OM2:
% $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ baseD = [baseL '1deg_jra55_ryf8485_kds50_may/']; %Data Directory.
% $$$ ICdir = '/g/data1/ua8/MOM/initial_conditions/WOA/10_KDS50/';
% MOM-SIS01:
model = 'MOM01';
baseD = [baseL 'MOM01_HeatDiag/']; %Data Directory.

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

% $$$ post = 'ocean/'; % For ACCESS-OM2 output coulpled;
post = ''; % For MOM-SIS.

haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
haveGM = 0; % 1 = GM is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;
haveMIX = 1; % 1 = Do mixing components (vdiffuse_diff_cbt_*), 0 = don't. 
haveHND = 1; % 1 = Do numerical mixing via heat budget.

% scaling constant on the transports:
if (strcmp(model(1),'A')) %ACCESS-OM2, transport in kg/s
    tsc = 1;
else % MOM-SIS, transport in 1e9 kg/s
    tsc = 1e9;
end

% $$$ for output = 86:90;
output=4;
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
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');

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
    [maskREG,~,~,~,~,~,~] = Heat_Budget_Mask(region,gname,fname, ...
                                             wname,outD,model);
else
    maskREG = ones(xL,yL);
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

%latitude vector for heat function:
latv = max(lat,[],1);
late = [-90 (latv(2:end)+latv(1:(end-1)))/2 90];

save([outD model sprintf('_output%03d',output) '_BaseVars.mat'], ...
     'T','Te','TL','dT','Cp','rho0','time','time_snap','tL', ...
     'z','zL','lon','lat','area','xL','yL','latv','late', ...
     'lonu','latu','-v7.3');

%% Calculate dVdt, dHdt and save back into wmass file:
%% Also, calculate TENMON

% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
try
    id = netcdf.inqVarID(ncid,'dVdt');
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutral');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    dVdtID = netcdf.defVar(ncid,'dVdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dVdtID,'long_name','Change in time of volume within temperature bin');
    netcdf.putAtt(ncid,dVdtID,'units','Sv (10^9 kg/s)');
    netcdf.putAtt(ncid,dVdtID,'_FillValue',single(-1e20));
    dHdtID = netcdf.defVar(ncid,'dHdt','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,dHdtID,'long_name','Change in time of heat content within temperature bin');
    netcdf.putAtt(ncid,dHdtID,'units','Watts');
    netcdf.putAtt(ncid,dHdtID,'_FillValue',single(-1e20));
    netcdf.endDef(ncid);
else
    dVdtID = netcdf.inqVarID(ncid,'dVdt');
    dHdtID = netcdf.inqVarID(ncid,'dHdt');
end

Vsnap = zeros(xL,yL,TL);
Hsnap = zeros(xL,yL,TL);

%Do IC for Vsnap and Hsnap:
for zi = 1:zL
    sprintf('Calculating Vsnap and Hsnap IC, depth %02d of %02d',zi,zL)
    %Temperature snapshot:
    tempsnap = ncread(rnameT,'temp',[1 1 zi rstti],[xL yL 1 1]);
    tempsnap(~mask(:,:,zi)) = NaN;
    if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
    
    if (found_rst)
        if (output == 0) % Initial dzt:
            Volsnap = dztI(:,:,zi).*area;
        else
            Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
        end
    else
        Volsnap = ncread(rnameT,'dzt',[1 1 zi rstti],[xL yL 1 1]).*area;
    end
    Volsnap(isnan(Volsnap)) = 0;
    
    for Ti=1:TL
        %Accumulate sums:
        inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
        Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
        Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
        Hlay(isnan(Hlay)) = 0;
        Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
    end
end

VsnapM = Vsnap;
HsnapM = Hsnap;
Vsnap = zeros(xL,yL,TL);
Hsnap = zeros(xL,yL,TL);

TENMON = zeros(TL+1,tL);

%Do other times for Vsnap and Hsnap:
for ti=1:tL
    for zi=1:zL
        sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        temp(~mask(:,:,zi)) = NaN;
        if (max(max(temp))>120);temp = temp-273.15;end;

        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap(~mask(:,:,zi)) = NaN;
        if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap(isnan(Volsnap)) = 0;
        
        TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL ...
                            yL 1 1]);

        for Ti=1:TL
            %Accumulate sums:
            inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
            Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
            Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
            Hlay(isnan(Hlay)) = 0;
            Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
            inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
            TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));
        end
        inds = find(temp>=Te(TL+1));
        TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));
    end
    
    % Integrate to get to T'>T:
    TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
    
    for Ti=1:TL
    netcdf.putVar(ncid,dVdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Vsnap(:,:,Ti)-VsnapM(:,:,Ti)) ...
                        /(time_snap(ti+1)-time_snap(ti))/86400*rho0/1e9);
    netcdf.putVar(ncid,dHdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Hsnap(:,:,Ti)-HsnapM(:,:,Ti)) ...
                        /(time_snap(ti+1)-time_snap(ti))/86400);
    end
    VsnapM = Vsnap;
    HsnapM = Hsnap;
    Vsnap = zeros(xL,yL,TL);
    Hsnap = zeros(xL,yL,TL);
end
netcdf.close(ncid);

% $$$ % Calculate non-snap zonal-average annual-average volumes:
% $$$ V = zeros(yL,TL,tL);
% $$$ H = zeros(yL,TL,tL);
% $$$ for ti=1:tL
% $$$     for zi=1:zL
% $$$         sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
% $$$ 
% $$$         temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         temp(~mask(:,:,zi)) = NaN;
% $$$         if (max(max(temp))>120);temp = temp-273.15;end;
% $$$         Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
% $$$         Vol(isnan(Vol)) = 0;
% $$$         
% $$$         for Ti=1:TL
% $$$             %Accumulate sums:
% $$$             inds = temp>=Te(Ti) & temp<Te(Ti+1);
% $$$             V(:,Ti,ti) = V(:,Ti,ti) + nansum(Vol.*inds,1)';
% $$$             Hlay = Vol.*temp.*inds*rho0*Cp;
% $$$             Hlay(isnan(Hlay)) = 0;
% $$$             H(:,Ti,ti) = H(:,Ti,ti) + nansum(Hlay,1)';
% $$$         end
% $$$     end
% $$$ end

%% Calculate numerical mixing by residual from heat budget and save
%% back into netcdf file:

% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
try
    id = netcdf.inqVarID(ncid,'temp_numdiff_heat_on_nrho');
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutralrho_edges');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    ndifID = netcdf.defVar(ncid,'temp_numdiff_heat_on_nrho','NC_FLOAT',[xid yid zid tid]);
    netcdf.putAtt(ncid,ndifID,'long_name',['Diffusion of heat due ' ...
                        'to numerical mixing estimated from heat ' ...
                        'fluxes binned to neutral density']);
    netcdf.putAtt(ncid,ndifID,'units','Watts/m^2');
    netcdf.putAtt(ncid,ndifID,'_FillValue',single(-1e20));
    netcdf.endDef(ncid);
else
    ndifID = netcdf.inqVarID(ncid,'temp_numdiff_heat_on_nrho');
end

for ti=1:tL
    Ti = TL;
    sprintf('Calculating numdif time %03d of %03d, Temp %03d of %03d',ti,tL,Ti,TL)
    JS = ncread(wname,'mass_pmepr_on_nrho',[1 1 Ti ti],[xL yL 1 1])./area/rho0;
    dVdt = double(ncread(wname,'dVdt',[1 1 Ti ti],[xL yL 1 1]))*1e9/rho0./area;
    dHdt = double(ncread(wname,'dHdt',[1 1 Ti ti],[xL yL 1 1]))./area;
    dift = ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
           ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    if (haveRedi)
        dift = dift + ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
    end
    if (haveMDS)
        dift = dift + ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end    
    
    txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
              ncread(wname,'tx_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
    tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
              ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
    if (haveGM)
        txtrans = txtrans + ncread(wname,'tx_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
    end
    qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])+ ...
              ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])+ ...
              ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    if (haveGM)
        qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);% + ...
% $$$                   ncread(wname,'temp_xflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);% + ...
% $$$                   ncread(wname,'temp_yflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end
    JI = zeros(xL,yL);
    JI(2:end,2:end) = (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                         +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
    JI(1,2:end) = (txtrans(end,2:end) - txtrans(1,2:end) ...
                                 +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
    QI = zeros(xL,yL);
    QI(2:end,2:end) = (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                         +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
    QI(1,2:end) = (qxtrans(end,2:end) - qxtrans(1,2:end) ...
                                 +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        

    ndif = dHdt - (dVdt - JI - JS)*rho0*Cp*Te(Ti) - dift - QI;
    
    netcdf.putVar(ncid,ndifID,[0 0 Ti ti-1],[xL yL 1 1],zeros(xL,yL));
    netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    for Ti=TL-1:-1:1
        sprintf(['Calculating numdif time %03d of %03d, Temp %03d ' ...
                 'of %03d'],ti,tL,Ti,TL)
        JS = JS + ncread(wname,'mass_pmepr_on_nrho',[1 1 Ti ti],[xL yL 1 1])./area/rho0;
        dVdt = dVdt + double(ncread(wname,'dVdt',[1 1 Ti ti],[xL ...
                            yL 1 1]))*1e9/rho0./area;
        dHdt = dHdt + double(ncread(wname,'dHdt',[1 1 Ti ti],[xL ...
                            yL 1 1]))./area;
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
                  ncread(wname,'tx_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
                  ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        if (haveGM)
            txtrans = txtrans + ncread(wname,'tx_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
            tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        end
        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])+ ...
                  ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])+ ...
                  ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveGM)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);% + ...
% $$$                   ncread(wname,'temp_xflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);% + ...
% $$$                   ncread(wname,'temp_yflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
            
        JI(2:end,2:end) = JI(2:end,2:end)+(txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end) = JI(1,2:end)+(txtrans(end,2:end) - txtrans(1,2:end) ...
            +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        

        QI(2:end,2:end) = QI(2:end,2:end)+(qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end) = QI(1,2:end)+(qxtrans(end,2:end) - qxtrans(1,2:end) ...
            +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        

        
        dift = dift + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
               ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveRedi)
            dift = dift + ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                          ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveMDS)
            dift = dift + ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
             
        ndif = dHdt - (dVdt - JI - JS)*rho0*Cp*Te(Ti) - dift - QI;
        netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    end
end
netcdf.close(ncid);

% $$$ %% Calculate Monthly and Annual binned volume integrated budget:
% $$$ 
% $$$ GWBmon.SWH    = zeros(TL+1,tL);GWBmon.VDS    = zeros(TL+1,tL);
% $$$ GWBmon.RMX    = zeros(TL+1,tL);GWBmon.PME    = zeros(TL+1,tL);
% $$$ GWBmon.FRZ    = zeros(TL+1,tL);GWBmon.ETS    = zeros(TL+1,tL);
% $$$ GWBmon.SUB    = zeros(TL+1,tL);GWBmon.VDF    = zeros(TL+1,tL);
% $$$ GWBmon.KNL    = zeros(TL+1,tL);GWBmon.ADV    = zeros(TL+1,tL);
% $$$ GWBmon.TEN    = zeros(TL+1,tL);GWBmon.SFW    = zeros(TL+1,tL);
% $$$ 
% $$$ GWBann = GWBmon;
% $$$ 
% $$$ for zi=1:zL
% $$$     tempAN = zeros(xL,yL);
% $$$     SWHAN = zeros(xL,yL);VDSAN = zeros(xL,yL);RMXAN = zeros(xL,yL);
% $$$     SUBAN = zeros(xL,yL);VDFAN = zeros(xL,yL);KNLAN = zeros(xL,yL);
% $$$     ADVAN = zeros(xL,yL);TENAN = zeros(xL,yL);SFWAN = zeros(xL,yL);
% $$$     PMEAN = zeros(xL,yL);FRZAN = zeros(xL,yL);ETSAN = zeros(xL,yL);
% $$$ 
% $$$     for ti=1:tL
% $$$         sprintf('Calculating MON/AN binned time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
% $$$ 
% $$$         temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
% $$$         temp(~mask(:,:,zi)) = NaN;
% $$$         if (max(max(temp))>120);temp = temp-273.15;end;
% $$$ 
% $$$         tempAN = tempAN + temp;
% $$$         
% $$$         TEN = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);
% $$$         ADV = area.*ncread(hname,'temp_advection',[1 1 zi ti],[xL yL 1 1]);
% $$$         SUB = area.*ncread(hname,'temp_submeso',[1 1 zi ti],[xL yL 1 1]);
% $$$         RMX = area.*ncread(hname,'temp_rivermix',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDS = area.*ncread(hname,'temp_vdiffuse_sbc',[1 1 zi ti],[xL yL 1 1]);
% $$$         SWH = area.*ncread(hname,'sw_heat',[1 1 zi ti],[xL yL 1 1]);
% $$$         VDF = area.*ncread(hname,'temp_vdiffuse_diff_cbt',[1 1 zi ti],[xL yL 1 1]);
% $$$         KNL = area.*ncread(hname,'temp_nonlocal_KPP',[1 1 zi ti],[xL yL 1 1]);
% $$$ 
% $$$         if (zi == 1)
% $$$             FRZ = area.*ncread(hname,'frazil_2d',[1 1 ti],[xL yL 1]);
% $$$             ETS = area.*ncread(hname,'temp_eta_smooth',[1 1 ti],[xL yL 1]);
% $$$             PME = area.*ncread(hname,'sfc_hflux_pme',[1 1 ti],[xL yL 1]);
% $$$         end
% $$$         
% $$$         TENAN = TENAN + TEN;
% $$$         ADVAN = ADVAN + ADV;
% $$$         SUBAN = SUBAN + SUB;
% $$$         RMXAN = RMXAN + RMX;
% $$$         VDSAN = VDSAN + VDS;
% $$$         SWHAN = SWHAN + SWH;
% $$$         VDFAN = VDFAN + VDF;
% $$$         KNLAN = KNLAN + KNL;
% $$$         
% $$$         if (zi == 1)
% $$$             PMEAN = PMEAN + PME;
% $$$             FRZAN = FRZAN + FRZ;
% $$$             ETSAN = ETSAN + ETS;
% $$$         end
% $$$         
% $$$         for Ti=1:TL
% $$$             %Accumulate sums:
% $$$             inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
% $$$             GWBmon.TEN(Ti,ti) = GWBmon.TEN(Ti,ti)+nansum(TEN(inds));
% $$$             GWBmon.ADV(Ti,ti) = GWBmon.ADV(Ti,ti)+nansum(ADV(inds));
% $$$             GWBmon.SUB(Ti,ti) = GWBmon.SUB(Ti,ti)+nansum(SUB(inds));
% $$$             GWBmon.RMX(Ti,ti) = GWBmon.RMX(Ti,ti)+nansum(RMX(inds));
% $$$             GWBmon.VDS(Ti,ti) = GWBmon.VDS(Ti,ti)+nansum(VDS(inds));
% $$$             GWBmon.SWH(Ti,ti) = GWBmon.SWH(Ti,ti)+nansum(SWH(inds));
% $$$             GWBmon.VDF(Ti,ti) = GWBmon.VDF(Ti,ti)+nansum(VDF(inds));
% $$$             GWBmon.KNL(Ti,ti) = GWBmon.KNL(Ti,ti)+nansum(KNL(inds));
% $$$             
% $$$             if (zi == 1)
% $$$                 GWBmon.FRZ(Ti,ti) = GWBmon.FRZ(Ti,ti)+nansum(FRZ(inds));
% $$$                 GWBmon.ETS(Ti,ti) = GWBmon.ETS(Ti,ti)+nansum(ETS(inds));
% $$$                 GWBmon.PME(Ti,ti) = GWBmon.PME(Ti,ti)+nansum(PME(inds));
% $$$             end
% $$$         end
% $$$         inds = find(temp>=Te(TL+1));
% $$$         GWBmon.TEN(TL+1,ti) = GWBmon.TEN(TL+1,ti)+nansum(TEN(inds));
% $$$         GWBmon.ADV(TL+1,ti) = GWBmon.ADV(TL+1,ti)+nansum(ADV(inds));
% $$$         GWBmon.SUB(TL+1,ti) = GWBmon.SUB(TL+1,ti)+nansum(SUB(inds));
% $$$         GWBmon.RMX(TL+1,ti) = GWBmon.RMX(TL+1,ti)+nansum(RMX(inds));
% $$$         GWBmon.VDS(TL+1,ti) = GWBmon.VDS(TL+1,ti)+nansum(VDS(inds));
% $$$         GWBmon.SWH(TL+1,ti) = GWBmon.SWH(TL+1,ti)+nansum(SWH(inds));
% $$$         GWBmon.VDF(TL+1,ti) = GWBmon.VDF(TL+1,ti)+nansum(VDF(inds));
% $$$         GWBmon.KNL(TL+1,ti) = GWBmon.KNL(TL+1,ti)+nansum(KNL(inds));
% $$$ 
% $$$         if (zi == 1)
% $$$             GWBmon.FRZ(TL+1,ti) = GWBmon.FRZ(TL+1,ti)+nansum(FRZ(inds));
% $$$             GWBmon.ETS(TL+1,ti) = GWBmon.ETS(TL+1,ti)+nansum(ETS(inds));
% $$$             GWBmon.PME(TL+1,ti) = GWBmon.PME(TL+1,ti)+nansum(PME(inds));
% $$$         end
% $$$     end
% $$$     
% $$$     ti = 1;
% $$$     for Ti=1:TL
% $$$         %Accumulate sums:
% $$$         inds = find(tempAN/tL>=Te(Ti) & tempAN/tL<Te(Ti+1));
% $$$         GWBann.TEN(Ti,ti) = GWBann.TEN(Ti,ti)+nansum(TENAN(inds)/tL);
% $$$         GWBann.ADV(Ti,ti) = GWBann.ADV(Ti,ti)+nansum(ADVAN(inds)/tL);
% $$$         GWBann.SUB(Ti,ti) = GWBann.SUB(Ti,ti)+nansum(SUBAN(inds)/tL);
% $$$         GWBann.RMX(Ti,ti) = GWBann.RMX(Ti,ti)+nansum(RMXAN(inds)/tL);
% $$$         GWBann.VDS(Ti,ti) = GWBann.VDS(Ti,ti)+nansum(VDSAN(inds)/tL);
% $$$         GWBann.SWH(Ti,ti) = GWBann.SWH(Ti,ti)+nansum(SWHAN(inds)/tL);
% $$$         GWBann.VDF(Ti,ti) = GWBann.VDF(Ti,ti)+nansum(VDFAN(inds)/tL);
% $$$         GWBann.KNL(Ti,ti) = GWBann.KNL(Ti,ti)+nansum(KNLAN(inds)/tL);
% $$$     
% $$$         if (zi == 1)
% $$$             GWBann.FRZ(Ti,ti) = GWBann.FRZ(Ti,ti)+nansum(FRZAN(inds)/tL);
% $$$             GWBann.ETS(Ti,ti) = GWBann.ETS(Ti,ti)+nansum(ETSAN(inds)/tL);
% $$$             GWBann.PME(Ti,ti) = GWBann.PME(Ti,ti)+nansum(PMEAN(inds)/tL);
% $$$         end
% $$$     end
% $$$     inds = find(temp>=Te(TL+1));
% $$$     GWBann.TEN(TL+1,ti) = GWBann.TEN(TL+1,ti)+nansum(TENAN(inds)/tL);
% $$$     GWBann.ADV(TL+1,ti) = GWBann.ADV(TL+1,ti)+nansum(ADVAN(inds)/tL);
% $$$     GWBann.SUB(TL+1,ti) = GWBann.SUB(TL+1,ti)+nansum(SUBAN(inds)/tL);
% $$$     GWBann.RMX(TL+1,ti) = GWBann.RMX(TL+1,ti)+nansum(RMXAN(inds)/tL);
% $$$     GWBann.VDS(TL+1,ti) = GWBann.VDS(TL+1,ti)+nansum(VDSAN(inds)/tL);
% $$$     GWBann.SWH(TL+1,ti) = GWBann.SWH(TL+1,ti)+nansum(SWHAN(inds)/tL);
% $$$     GWBann.VDF(TL+1,ti) = GWBann.VDF(TL+1,ti)+nansum(VDFAN(inds)/tL);
% $$$     GWBann.KNL(TL+1,ti) = GWBann.KNL(TL+1,ti)+nansum(KNLAN(inds)/tL);
% $$$ 
% $$$     if (zi == 1)
% $$$         GWBann.FRZ(TL+1,ti) = GWBann.FRZ(TL+1,ti)+nansum(FRZAN(inds)/tL);
% $$$         GWBann.ETS(TL+1,ti) = GWBann.ETS(TL+1,ti)+nansum(ETSAN(inds)/tL);
% $$$         GWBann.PME(TL+1,ti) = GWBann.PME(TL+1,ti)+nansum(PMEAN(inds)/tL);
% $$$     end
% $$$ end
% $$$ 
% $$$ % Integrate to get to T'>T:
% $$$ names = fieldnames(GWBmon);
% $$$ for i=1:length(names)
% $$$     for ti=1:tL
% $$$         eval(['GWBmon.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBmon.' ...
% $$$               names{i} '(:,ti))));']);
% $$$     end
% $$$     ti = 1;
% $$$     eval(['GWBann.' names{i} '(:,ti) = flipud(cumsum(flipud(GWBann.' ...
% $$$           names{i} '(:,ti))));']);
% $$$ end
% $$$ 
% $$$ save([outD model sprintf('_output%03d',output) '_GlobalHBud_MonAnBin.mat'],'GWBmon','GWBann','-v7.3');
% $$$ 

%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------

GWB.dVdt   = zeros(TL+1,tL); % Sv 
try
    GWB.TENMON = TENMON; % Tendency from monthly averages (from previous code block)
catch
    GWB.TENMON = GWB.dVdt;
end
GWB.dHdt   = zeros(TL+1,tL); % Total W
GWB.SWH    = zeros(TL+1,tL); % W due to SW redistribution
GWB.VDS    = zeros(TL+1,tL); % W due to vdiffuse_sbc.
GWB.RMX    = zeros(TL+1,tL); % W due to rivermix.
GWB.PME    = zeros(TL+1,tL); % W due to P-E.
GWB.FRZ    = zeros(TL+1,tL); % W due to frazil.
GWB.ETS    = zeros(TL+1,tL); % W due to eta_smoothing.
GWB.SUB    = zeros(TL+1,tL); % W due to submesoscale.
GWB.VDF    = zeros(TL+1,tL); % W due to vdiffusion
GWB.KNL    = zeros(TL+1,tL); % W due to KPP non-local
if (haveMIX)
    GWB.VDFkppiw = zeros(TL+1,tL); % W due to KPP Internal Wave
    GWB.VDFkppish = zeros(TL+1,tL); % W due to KPP Shear
    GWB.VDFkppicon = zeros(TL+1,tL); % W due to KPP Convection
    GWB.VDFkppbl = zeros(TL+1,tL); % W due to KPP Boundary Layer
    GWB.VDFkppdd = zeros(TL+1,tL); % W due to KPP Double-Diffusion
    GWB.VDFwave = zeros(TL+1,tL); % W due to Simmons Internal Wave
end    
if (haveRedi)
    GWB.K33    = zeros(TL+1,tL); % W due to K33
    GWB.RED    = zeros(TL+1,tL); % W due to Redi diffusion
end
if (haveGM)
    GWB.NGM    = zeros(TL+1,tL); % W due to GM
end
if (haveMDS)
    GWB.MDS    = zeros(TL+1,tL); % W due to mixdownslope
end
GWB.ADV    = zeros(TL+1,tL); % W due to advection
GWB.TEN    = zeros(TL+1,tL); % W due to tendency
GWB.SFW    = zeros(TL+1,tL); % surface volume flux into ocean (m3s-1)

if (haveHND)
    GWB.NUM    = zeros(TL+1,tL); % W due to numerical mixing from heat budget
    % Get NUM:
    for ti=1:tL
        for Ti = (TL+1):-1:1
            sprintf('Calculating NUMH heat budget time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
            GWB.NUM(Ti,ti) = nansum(nansum(area.*ncread(wname, ...
                                                        'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]),1),2);
        end
    end
end

for ti=1:tL
    ii = TL;
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveMIX)
    GWB.VDFkppiw(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end        
    if (haveRedi)
    GWB.K33(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = nansum(nansum(maskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = nansum(nansum(maskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SUB(ii,ti) = GWB.SUB(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.PME(ii,ti) = GWB.PME(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RMX(ii,ti) = GWB.RMX(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDS(ii,ti) = GWB.VDS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SWH(ii,ti) = GWB.SWH(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.KNL(ii,ti) = GWB.KNL(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveMIX)
    GWB.VDFkppiw(ii,ti)   = GWB.VDFkppiw(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppish(ii,ti)  = GWB.VDFkppish(ii+1,ti)  + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppicon(ii,ti) = GWB.VDFkppicon(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppbl(ii,ti)   = GWB.VDFkppbl(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFkppdd(ii,ti)   = GWB.VDFkppdd(ii+1,ti)   + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDFwave(ii,ti)    = GWB.VDFwave(ii+1,ti)    + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);        
    end
    if (haveRedi)
    GWB.K33(ii,ti) = GWB.K33(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = GWB.RED(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = GWB.NGM(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveMDS)
    GWB.MDS(ii,ti) = GWB.MDS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.FRZ(ii,ti) = GWB.FRZ(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ETS(ii,ti) = GWB.ETS(ii+1,ti) + nansum(nansum(maskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.SFW(ii,ti) = GWB.SFW(ii+1,ti) + nansum(nansum(maskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
end
end
save([outD model sprintf('_output%03d',output) '_' region 'HBud.mat'],'GWB','-v7.3');

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
Tls = [0:2.5:27.5];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlF = zeros(xL,yL,tL); % surface forcing
FlP = zeros(xL,yL,tL); % P-E+R
% $$$ FlA = zeros(xL,yL,tL); % advection + submeso + GM
if (haveMIX)
    FlMdif  = zeros(xL,yL,tL); % vdiffuse
    FlMkppiw  = zeros(xL,yL,tL);
    FlMkppish  = zeros(xL,yL,tL);
    FlMkppicon  = zeros(xL,yL,tL);
    FlMkppbl  = zeros(xL,yL,tL);
    FlMkppdd  = zeros(xL,yL,tL);
    FlMwave  = zeros(xL,yL,tL);
end
if (haveRedi)
    FlK = zeros(xL,yL,tL); % K33
    FlR = zeros(xL,yL,tL); % Redi
end
% $$$ if (haveGM)
% $$$     FlG = zeros(xL,yL,tL); % GM
% $$$ end
% $$$ FlT = zeros(xL,yL,tL); % tendency
FlSP = zeros(xL,yL,tL); % solar penetration
FlI = zeros(xL,yL,tL); % numerical mixing
 
while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
% $$$         FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveRedi)
            FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
% $$$         FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$             ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
% $$$         if (haveGM)
% $$$             FlG(:,:,ti) = FlG(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
% $$$             FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
% $$$         end
        FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveMDS)
            FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveMIX)
            FlMdif(:,:,ti) = FlMdif(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppiw(:,:,ti) = FlMkppiw(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppish(:,:,ti) = FlMkppish(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppicon(:,:,ti) = FlMkppicon(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppbl(:,:,ti) = FlMkppbl(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMkppdd(:,:,ti) = FlMkppdd(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            FlMwave(:,:,ti) = FlMwave(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end            
        FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
            ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlSP(:,:,ti) = FlSP(:,:,ti)+ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlI(:,:,ti) = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlSP','FlF','FlP','FlI','Tl','-v7.3');
% $$$         save(name,'FlM','FlSP','FlF','FlT','FlA','FlP','FlI','FlIH','Tl','-v7.3');
        if (haveMIX)
            save(name,'FlMdif','FlMkppiw','FlMkppish', ...
                 'FlMkppicon','FlMkppbl','FlMkppdd','FlMwave','-append');
        end
        if (haveRedi)
            save(name,'FlK','FlR','-append');
        end
% $$$         if (haveGM)
% $$$             save(name,'FlG','-append');
% $$$         end
        
        Nremain = Nremain-1;
    end

    Ti = Ti-1;
end

% $$$ %% Calculate WMT due to different (resolved) terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ Tls = [0:2.5:27.5]-0.25;
% $$$ 
% $$$ for ii = 1:length(Tls)
% $$$     Tl = Tls(ii);
% $$$ 
% $$$     WMTM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
% $$$     WMTF = NaN*zeros(xL,yL,tL); % surface forcing
% $$$     WMTI = NaN*zeros(xL,yL,tL); % numerical
% $$$     WMTP = NaN*zeros(xL,yL,tL); % P-E+R
% $$$     WMTSP = NaN*zeros(xL,yL,tL); % solar penetration
% $$$     if (haveRedi)
% $$$         WMTK = NaN*zeros(xL,yL,tL); % K33
% $$$         WMTR = NaN*zeros(xL,yL,tL); % Redi
% $$$     end
% $$$     T = ncread(wname,'neutral');
% $$$     Te = ncread(wname,'neutralrho_edges');
% $$$     [tmp Ti] = min(abs(T-Tl));
% $$$     
% $$$     for ti=1:tL
% $$$         sprintf('Calculating WMT time %03d of %03d, temp %03d of %03d',ti,tL,ii,length(Tls))
% $$$         WMTP(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$                                       ncread(wname,'temp_rivermix_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$         WMTM(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$                                       ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$         WMTF(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$                                       ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$                                       ncread(wname,'frazil_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
% $$$                                       ncread(wname,'temp_eta_smooth_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$         WMTSP(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'sw_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$ % $$$         WMTI(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$         if (haveRedi)
% $$$             WMTK(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]));
% $$$             WMTR(:,:,ti) = 1/rho0/Cp/dT*(ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]));
% $$$         end
% $$$         if (haveMDS)
% $$$             WMTM(:,:,ti) = WMTM(:,:,ti) + 1/rho0/Cp/dT* ...
% $$$                 ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
% $$$         end
% $$$     end
% $$$     save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTM','WMTSP','WMTP','WMTF','WMTI','Tl','-v7.3');
% $$$     if (haveRedi)
% $$$         save([outD model sprintf('_output%03d',output) '_WMT_T' strrep(num2str(Tl),'.','p') 'C.mat'],'WMTK','WMTR','-append');
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ %% Save isotherm depths -------------------------------------------------------------------------------------
% $$$ Tls = [22 22.5 23];
% $$$ 
% $$$ for ii = 1:length(Tls)
% $$$     Tl = Tls(ii);
% $$$ 
% $$$     ziso = NaN*zeros(xL,yL,tL);
% $$$ 
% $$$     for ti=1:tL
% $$$         sprintf('Calculating isotherm depth time %03d of %03d, temp %03d of %03d',ti,tL,ii,length(Tls))
% $$$         temp = reshape(ncread(fname,'temp',[1 1 1 ti],[xL yL zL 1]),[xL*yL zL 1]);
% $$$         tmax = max(temp,[],2);
% $$$         inds = tmax >= Tl;
% $$$         temp = temp(inds,:);
% $$$         temp(isnan(temp)) = -1000;
% $$$         temp = temp -0.001*repmat(1:zL,[length(temp(:,1)) 1]);
% $$$         zisot = zeros(length(temp(:,1)),1);
% $$$         for jj=1:length(temp(:,1))
% $$$             zisot(jj) = interp1(temp(jj,:),z,Tl,'linear');
% $$$         end
% $$$         tmp = NaN*zeros(xL,yL);
% $$$         tmp(inds) = zisot;
% $$$         ziso(:,:,ti) = tmp;
% $$$     end
% $$$     save([outD model sprintf('_output%03d',output) '_Ziso_T' strrep(num2str(Tl),'.','p') 'C.mat'],'ziso');
% $$$ end
% $$$ 

%% Save surface heat flux, wind stress, SST, meridional heat flux:
try
    shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tL]);
catch
    shflux = ncread(fname2,'net_sfc_heating',[1 1 1],[xL yL tL]);
end    
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
% $$$ taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
% $$$ tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);

save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST');%,'taux','tauy');

% $$$ % Do meridional heat flux (note: Inferred from transports! Not
% $$$ % fully accurate):
% $$$ mhflux = zeros(yL,tL);
% $$$ for ti = 1:tL
% $$$     for Ti = 1:TL
% $$$         sprintf('Calculating meridional heat flux time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
% $$$         tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0 + ...
% $$$                   ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
% $$$         if (haveGM)
% $$$             tytrans = tytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
% $$$         end
% $$$         mhflux(:,ti) = mhflux(:,ti) + rho0*Cp*T(Ti)*nansum(tytrans,1)';
% $$$     end
% $$$ end
% $$$ save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'mhflux','-append');

%% Zonally-averaged fluxes -------------------------------------------------------------
ZA.F = zeros(yL,TL+1,tL); % Wdeg-1 due to F
ZA.P = zeros(yL,TL+1,tL); % Wdeg-1 due to P
ZA.M = zeros(yL,TL+1,tL); % Wdeg-1 due to M
ZA.dVdt = zeros(yL,TL+1,tL); % Wdeg-1 dVdt
ZA.dHdt = zeros(yL,TL+1,tL); % Wdeg-1 dHdt
ZA.SWH = zeros(yL,TL+1,tL); % Wdeg-1 due to SW redistribution
if (haveMIX)
    ZA.Mkppiw = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mkppish = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mwave = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Mkppbl = zeros(yL,TL+1,tL); %Wdeg-1 due to kppiw;
    ZA.Moth   = zeros(yL,TL+1,tL); %Wdeg-1 due to kppicon+kppdd+KPPnloc;
end
if (haveRedi)
    ZA.K33    = zeros(yL,TL+1,tL); % W due to K33
    ZA.RED    = zeros(yL,TL+1,tL); % W due to Redi diffusion
end
ZA.JS = zeros(yL,TL+1,tL); % m3s-1deg-1 due to surface volume flux
ZA.PSI = zeros(yL,TL+1,tL); % m3s-1 northward transport
ZA.AHD = zeros(yL,TL+1,tL); % W A direct using heat fluxes

yto = diff(ncread(gname,'yt_ocean'));
yto = [yto(1); (yto(2:end)+yto(1:(end-1)))/2; yto(end)];
yuo = diff([-90; ncread(gname,'yu_ocean')]);

% ignoring tri-polar for now.
for ti=1:tL
    for ii = 1:TL
    sprintf('Calculating zonally-averaged water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    ZA.F(:,ii+1,ti) = ZA.F(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(maskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(maskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
    ZA.P(:,ii+1,ti) = ZA.P(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(maskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
    ZA.M(:,ii+1,ti) = ZA.M(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(maskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
    if (haveMDS)
        ZA.M(:,ii+1,ti) = ZA.M(:,ii+1,ti) + nansum(maskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)'./yuo;
    end
    ZA.dVdt(:,ii+1,ti) = ZA.dVdt(:,ii,ti) + nansum(maskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1)'./yuo;
    ZA.dHdt(:,ii+1,ti) = ZA.dHdt(:,ii,ti) + nansum(maskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1)'./yuo;
    ZA.PSI(:,ii+1,ti) =  ZA.PSI(:,ii,ti) + nansum(maskREG.*ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]) + ...
                              maskREG.*ncread(wname,'ty_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
    if (haveGM) 
        ZA.PSI(:,ii+1,ti) = ZA.PSI(:,ii+1,ti) + nansum(maskREG.*ncread(wname,'ty_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
    end
    if (haveHND)
        ZA.AHD(:,ii+1,ti) = ZA.AHD(:,ii,ti) + nansum(maskREG.*ncread(wname,'temp_yflux_adv_on_nrho',[1 1 ii ti],[xL yL 1 1])+ ...
                                 maskREG.*ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';

        if (haveGM)
            ZA.AHD(:,ii+1,ti) = ZA.AHD(:,ii+1,ti) + nansum(maskREG.*ncread(wname,'temp_yflux_gm_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
        end
    end

    ZA.SWH(:,ii+1,ti) = ZA.SWH(:,ii,ti) + nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)'./yuo;
    ZA.F(:,ii+1,ti) = ZA.F(:,ii+1,ti) + nansum(maskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)'./yuo;
    if (haveMIX)
        ZA.Mkppiw(:,ii+1,ti) = ZA.Mkppiw(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
        ZA.Mkppish(:,ii+1,ti) = ZA.Mkppish(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
        ZA.Mwave(:,ii+1,ti) = ZA.Mwave(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
        ZA.Mkppbl(:,ii+1,ti) = ZA.Mkppbl(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
        ZA.Moth(:,ii+1,ti) = ZA.Moth(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(maskREG.*area.*ncread(wname, ...
                                                'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
    end
    if (haveRedi)
        ZA.K33(:,ii+1,ti) = ZA.K33(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
        ZA.RED(:,ii+1,ti) = ZA.RED(:,ii,ti) + (nansum(maskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1))'./yuo;
    end        
    ZA.JS(:,ii+1,ti) = ZA.JS(:,ii,ti) + (nansum(maskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)/rho0)'./yuo;
    end
    save([outD model sprintf('_output%03d',output) '_' region '_ZAHBud.mat'],'ZA','yuo','yto','-v7.3');
end

% $$$ %% Annual average/max, zonal-average/max depth of isotherms:
% $$$ 
% $$$ tempZA = zeros(yL,zL);
% $$$ saltZA = zeros(yL,zL);
% $$$ tempZM = -100*zeros(yL,zL);
% $$$ tempZMA = -100*zeros(yL,zL);
% $$$ 
% $$$ for ti=1:tL
% $$$     sprintf('Calculating zonally-averaged temperature %03d of %03d',ti,tL)
% $$$     temp = ncread(fname,'temp',[1 1 1 ti],[xL yL zL 1]);
% $$$     salt = ncread(fname,'salt',[1 1 1 ti],[xL yL zL 1]);
% $$$     areaNAN = repmat(area,[1 1 zL]).*(~isnan(temp));
% $$$     
% $$$     tempZ = squeeze(nansum(areaNAN.*temp,1))./squeeze(nansum(areaNAN,1));
% $$$     tempZA = tempZA + tempZ;
% $$$ 
% $$$     saltZ = squeeze(nansum(areaNAN.*salt,1))./squeeze(nansum(areaNAN,1));
% $$$     saltZA = saltZA + saltZ;
% $$$ 
% $$$     maxt = squeeze(max(temp,[],1));
% $$$     tempZM = max(tempZM,maxt);
% $$$     tempZMA = max(tempZMA,tempZ);
% $$$ end
% $$$ tempZA = tempZA/tL;
% $$$ saltZA = saltZA/tL;
% $$$ 
% $$$ % Depth of isotherms:
% $$$ ZAtemp = zeros(yL,TL+1);
% $$$ ZMtemp = zeros(yL,TL+1);
% $$$ ZMAtemp = zeros(yL,TL+1);
% $$$ 'Calculating zonally-averaged isotherms'
% $$$ for yi=1:yL
% $$$     tvec = squeeze(tempZA(yi,:));
% $$$     zvec = -z;
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZAtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZAtemp(yi,:)),1,'last');
% $$$     ZAtemp(yi,(ind+1):end) = max(zvec);
% $$$ 
% $$$     tvec = squeeze(tempZM(yi,:));
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZMtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZMtemp(yi,:)),1,'last');
% $$$     ZMtemp(yi,(ind+1):end) = max(zvec);
% $$$ 
% $$$     tvec = squeeze(tempZMA(yi,:));
% $$$     tvec(isnan(tvec)) = -1000;
% $$$     tvec = tvec - 0.01*(1:zL);
% $$$     ZMAtemp(yi,:) = interp1(tvec,zvec,Te,'linear');
% $$$     ind = find(~isnan(ZMAtemp(yi,:)),1,'last');
% $$$     ZMAtemp(yi,(ind+1):end) = max(zvec);
% $$$ end
% $$$ 
% $$$ save([outD model sprintf('_output%03d',output) '_ZAHBud.mat'],'tempZA','saltZA','ZAtemp','tempZM','ZMtemp','tempZMA','ZMAtemp','z','Te','latv','-append');

% $$$ end

% $$$ %% Swap in non-NaN'd lon/lat:
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',8)]);
% $$$ region = 'Global';
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ model = 'MOM025_btide';
% $$$ for output = [20:21]
% $$$     save([base model sprintf('_output%03d_BaseVars.mat',output)], ...
% $$$          'lon','lat','lonu','latu','area','-append');
% $$$ end
% $$$ 
% $$$ 

