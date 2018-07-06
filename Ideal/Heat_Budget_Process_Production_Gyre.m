 % This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseL = '/short/e14/rmh561/mom/archive/';

model = 'MOM_Gyre';
baseD = [baseL 'MOM_Gyre/']; %Data Directory.

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

tsc=1e9;

for output = 1;
restart = output-1;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output)];
basem1 = [baseD sprintf('output%03d/',output-1)];
baser = [rstbaseD sprintf('restart%03d/',restart)];
hname = [base 'ocean_heat.nc'];
fname = [base 'ocean.nc'];
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

% 3D mask ------------------------------------------------
mask = ncread(fname,'temp',[1 1 1 rstti],[xL yL zL 1]);
mask(~isnan(mask)) = 1; mask(isnan(mask)) = 0;
mask = mask == 1;

% Time  -----------------------------------------
try
    time = ncread(fname,'time');
catch
    time = ncread(fname,'Time');
end

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
        Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
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
    
    netcdf.putVar(ncid,dVdtID,[0 0 0 ti-1],[xL yL TL 1],(Vsnap-VsnapM) ...
                  /(time_snap(ti+1)-time_snap(ti))/86400*rho0/1e9);
    netcdf.putVar(ncid,dHdtID,[0 0 0 ti-1],[xL yL TL 1],(Hsnap-HsnapM) ...
                  /(time_snap(ti+1)-time_snap(ti))/86400);
    VsnapM = Vsnap;
    HsnapM = Hsnap;
    Vsnap = zeros(xL,yL,TL);
    Hsnap = zeros(xL,yL,TL);
end
netcdf.close(ncid);

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
    dVdt = double(ncread(wname,'dVdt',[1 1 Ti ti],[xL yL 1 1]))*1e9/rho0./area;
    dHdt = double(ncread(wname,'dHdt',[1 1 Ti ti],[xL yL 1 1]))./area;
    dift = ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    
    txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
    tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
    qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);

    JI = zeros(xL,yL);
    JI(2:end,2:end) = (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                         +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
    JI(1,2:end) = (- txtrans(1,2:end) ...
                                 +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
    JI(2:end,1) = (txtrans(1:(end-1),1) - txtrans(2:end,1) ...
                       - tytrans(2:end,1))./area(2:end,1);        
    JI(1,1) = (-txtrans(1,1)-tytrans(1,1))./area(1,1);
    
    QI = zeros(xL,yL);
    QI(2:end,2:end) = (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                         +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
    QI(1,2:end) = (- qxtrans(1,2:end) ...
                                 +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        
    QI(2:end,1) = (qxtrans(1:(end-1),1) - qxtrans(2:end,1) ...
                       - qytrans(2:end,1))./area(2:end,1);        
    QI(1,1) = (-qxtrans(1,1)-qytrans(1,1))./area(1,1);

    ndif = dHdt - (dVdt - JI)*rho0*Cp*Te(Ti) - dift - QI;
    
    netcdf.putVar(ncid,ndifID,[0 0 Ti ti-1],[xL yL 1 1],zeros(xL,yL));
    netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    for Ti=TL-1:-1:1
        sprintf(['Calculating numdif time %03d of %03d, Temp %03d ' ...
                 'of %03d'],ti,tL,Ti,TL)
        dVdt = dVdt + double(ncread(wname,'dVdt',[1 1 Ti ti],[xL ...
                            yL 1 1]))*1e9/rho0./area;
        dHdt = dHdt + double(ncread(wname,'dHdt',[1 1 Ti ti],[xL ...
                            yL 1 1]))./area;
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            
        JI(2:end,2:end) = JI(2:end,2:end)+(txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end) = JI(1,2:end) + (- txtrans(1,2:end) ...
                                 +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
        JI(2:end,1) = JI(2:end,1) + (txtrans(1:(end-1),1) - txtrans(2:end,1) ...
                       - tytrans(2:end,1))./area(2:end,1);        
        JI(1,1) = JI(1,1) + (-txtrans(1,1)-tytrans(1,1))./area(1,1);

        QI(2:end,2:end) = QI(2:end,2:end)+(qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end) = QI(1,2:end) + (- qxtrans(1,2:end) ...
                       +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        
        QI(2:end,1) = QI(2:end,1) + (qxtrans(1:(end-1),1) - qxtrans(2:end,1) ...
                       - qytrans(2:end,1))./area(2:end,1);        
        QI(1,1) = QI(1,1) + (-qxtrans(1,1)-qytrans(1,1))./area(1,1);
        
        dift = dift + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
             
        ndif = dHdt - (dVdt - JI)*rho0*Cp*Te(Ti) - dift - QI;
        netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    end
end
netcdf.close(ncid);

%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------

GWB.TENMON = TENMON; % Tendency from monthly averages (from
                     % previous code block
GWB.dVdt   = zeros(TL+1,tL); % Sv 
GWB.dHdt   = zeros(TL+1,tL); % Total W
GWB.VDF    = zeros(TL+1,tL); % W due to vdiffusion
GWB.ADV    = zeros(TL+1,tL); % W due to advection
GWB.TEN    = zeros(TL+1,tL); % W due to tendency
GWB.NUM    = zeros(TL+1,tL); % W due to numerical mixing from heat budget

% Get NUM:
for ti=1:tL
    for Ti = (TL+1):-1:1
        sprintf('Calculating NUMH heat budget time %03d of %03d, temp %03d of %03d',ti,tL,Ti,TL)
        GWB.NUM(Ti,ti) = nansum(nansum(area.*ncread(wname, ...
                                                     'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]),1),2);
    end
end

for ti=1:tL
    ii = TL;
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = nansum(nansum(ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
end
end
save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'],'GWB','-v7.3');

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
Tls = [0:2:24];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlA = zeros(xL,yL,tL); % advection + submeso + GM
FlT = zeros(xL,yL,tL); % tendency
FlI = zeros(xL,yL,tL);

while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
        FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        FlI(:,:,ti) = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlT','FlA','FlI','Tl','-v7.3');
        Nremain = Nremain-1;
    end
    Ti = Ti-1;    
end

end

