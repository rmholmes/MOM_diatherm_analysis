% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseL = '/scratch/e14/rmh561/mom/archive/';
baseL = '/srv/ccrc/data03/z3500785/mom/';

% MOM-SIS025:
model = 'MOM_Gyre_Run014';
baseD = [baseL 'MOM_Gyre/Run004/']; %Data Directory.
% $$$ ICdir = ['/scratch/e14/rmh561/mom/input/gyre1/'];

outD = [baseD 'mat_data/'];
rstbaseD = baseD;

post = ''; % For MOM-SIS.

% term options:
haveSF = 1; % 1= have surface forcing (vdiffuse_sbc)
haveRedi = 0; % 1 = Redi diffusion is on, 0 = off
haveKPP = 0; % 1 = KPP is on (need nonlocal term)
haveGM = 0; % 1 = GM is on, 0 = off;
haveLT  = 0; % 1 = Laplacian lateral tracer diffusion on. 
haveSUB = 0; % 1 = submeso is on, 0 = off;
haveMDS = 0; % 1 = MDS is on, 0 = off;
haveSIG = 0; % 1 = SIG is on, 0 = off;
haveMIX = 0; % 1 = Do mixing components (vdiffuse_diff_cbt_*), 0 = don't. 

% Processing options:
doBASE     = 0; % 1 = save BaseVars.mat file
dodVdtdHdt = 0; % 1 = calculate dVdt/dHdt and save into .nc file
doVHza     = 0; % 1 = save zonally-integrated V and H fields from
                % average time slots in a .mat file.
doNUMDIF   = 0; % 1 = calculate tempdiff x,y,T,t and save into .nc file
doSGMviac  = 0; % 1 = calculate SUB/GM influence via binned
                % convergence (otherwise uses lateral flux). The
                % better option is to use the lateral flux -> This
                % calculates the net numerical mixing (i.e. it
                % includes the numerical mixing associated with the GM
                % and SUB schemes).

doGWB      = 1; % 1 = calculate global online budget terms
doXY       = 1; % 1 = calculate spatial fluxes-on-an-isotherm
doWMT      = 0; % 1 = calculate WMT volume fluxes spatial structure
doSURF     = 0; % 1 = calculate surface flux field and SST
doZA       = 0; % 1 = calculate zonal average budget

doHND      = 1; % 1 = calculate global online numdif
doTENMON   = 0; % 1 = do monthly eulerian tendency binning
doMONANN   = 0; % 1 = calculate monthly and annually binned eulerian global budget
doXYall    = 0; % 1 = do all XY calcs (most not used)
doXYtran   = 0; % 1 = calculate vertically-integrated heat
                % transports below given isotherm/s.
tsc = 1e9;

% $$$ for output = 0:5;
restart = output-1;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];
if (output==0)
    baser = ICdir;
else
    baser = [rstbaseD sprintf('restart%03d/',restart) post];
end
hname = [base 'ocean_heat.nc'];
fname_month = [base 'ocean_month.nc'];
fname = [base 'ocean.nc'];
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
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

tmaskREG = ones(xL,yL);
umaskREG = ones(xL,yL);

% Time  -----------------------------------------
time = ncread(wname,'time');
ndays = ncread(wname,'average_DT');
tL = length(time);
tLoc = length(ncread(fname,'Time'));

if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
else
    found_rst = 0;rstti = tL;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
end

Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

if (doBASE)
save([outD model sprintf('_output%03d',output) '_BaseVars.mat'], ...
     'T','Te','TL','dT','Cp','rho0','time','ndays','tL', ...
     'z','zL','lon','lat','area','xL','yL','yt','yu', 'xt','xu','lonu','latu','-v7.3');
end

%% Calculate dVdt, dHdt and save back into wmass file:
%% Also, calculate TENMON
if (dodVdtdHdt)
% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
try
    id = netcdf.inqVarID(ncid,'dHdt');
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutral');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    dVdtID = netcdf.defVar(ncid,'dVdt','NC_DOUBLE',[xid yid zid tid]);
    netcdf.putAtt(ncid,dVdtID,'long_name','Change in time of volume within temperature bin');
    netcdf.putAtt(ncid,dVdtID,'units','Sv (10^9 kg/s)');
% $$$     netcdf.putAtt(ncid,dVdtID,'_FillValue',single(-1e20));
    dHdtID = netcdf.defVar(ncid,'dHdt','NC_DOUBLE',[xid yid zid tid]);
    netcdf.putAtt(ncid,dHdtID,'long_name','Change in time of heat content within temperature bin');
    netcdf.putAtt(ncid,dHdtID,'units','Watts');
% $$$     netcdf.putAtt(ncid,dHdtID,'_FillValue',single(-1e20));
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

if (doTENMON)
    TENMON = zeros(TL+1,tL);
end

%Do other times for Vsnap and Hsnap:
for ti=1:tL
    for zi=1:zL
        sprintf('Calculating Vsnap and Hsnap later months time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)
        
        if (doTENMON)
            temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
            temp(~mask(:,:,zi)) = NaN;
            if (max(max(temp))>120);temp = temp-273.15;end;
        end

        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap(~mask(:,:,zi)) = NaN;
        if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap(isnan(Volsnap)) = 0;

        if (doTENMON)
            TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL ...
                                yL 1 1]);
        end
        
        for Ti=1:TL
            %Accumulate sums:
            inds = tempsnap>=Te(Ti) & tempsnap<Te(Ti+1);
            Vsnap(:,:,Ti) = Vsnap(:,:,Ti) + Volsnap.*inds;
            Hlay = Volsnap.*tempsnap.*inds*rho0*Cp;
            Hlay(isnan(Hlay)) = 0;
            Hsnap(:,:,Ti) = Hsnap(:,:,Ti) + Hlay;
            if (doTENMON)
                inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
                TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));
            end
        end
        if (doTENMON)
            inds = find(temp>=Te(TL+1));
            TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));
        end
    end
    
    if (doTENMON)
        % Integrate to get to T'>T:
        TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
    end
    
    for Ti=1:TL
    netcdf.putVar(ncid,dVdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Vsnap(:,:,Ti)-VsnapM(:,:,Ti)) ...
                        /ndays(ti)/86400*rho0/1e9);
    netcdf.putVar(ncid,dHdtID,[0 0 Ti-1 ti-1],[xL yL 1 1],(Hsnap(:,:,Ti)-HsnapM(:,:,Ti)) ...
                        /ndays(ti)/86400);
    end
    VsnapM = Vsnap;
    HsnapM = Hsnap;
    Vsnap = zeros(xL,yL,TL);
    Hsnap = zeros(xL,yL,TL);
end
netcdf.close(ncid);

end

if (doVHza)
    % Calculate non-snap zonal-average annual-average volumes:
    V = zeros(yL,TL,tL);
    H = zeros(yL,TL,tL);
    for ti=1:tL
        for zi=1:zL
            sprintf('Calculating V and H time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

            temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
            temp(~mask(:,:,zi)) = NaN;
            if (max(max(temp))>120);temp = temp-273.15;end;
            Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
            Vol(isnan(Vol)) = 0;
        
            for Ti=1:TL
                %Accumulate sums:
                inds = temp>=Te(Ti) & temp<Te(Ti+1);
                V(:,Ti,ti) = V(:,Ti,ti) + nansum(Vol.*inds,1)';
                Hlay = Vol.*temp.*inds*rho0*Cp;
                Hlay(isnan(Hlay)) = 0;
                H(:,Ti,ti) = H(:,Ti,ti) + nansum(Hlay,1)';
            end
        end
    end
    save([outD model sprintf('_output%03d',output) '_VHza.mat'],'V','H','-v7.3');
end


%% Calculate numerical mixing by residual from heat budget and save
%% back into netcdf file:
if (doNUMDIF)
% Create variables:
ncid = netcdf.open(wname,'NC_WRITE');
varname = 'temp_numdiff_heat_on_nrho';
try
    id = netcdf.inqVarID(ncid,varname);
    not_there = 0;
catch
    not_there = 1;
end
%If variable not there, add it:
if (not_there)
    xid = netcdf.inqDimID(ncid,'grid_xt_ocean');yid = netcdf.inqDimID(ncid,'grid_yt_ocean');zid = netcdf.inqDimID(ncid,'neutralrho_edges');tid = netcdf.inqDimID(ncid,'time');
    netcdf.reDef(ncid);
    ndifID = netcdf.defVar(ncid,varname,'NC_DOUBLE',[xid yid zid tid]);
    netcdf.putAtt(ncid,ndifID,'long_name',['Diffusion of heat due ' ...
                        'to numerical mixing estimated from heat ' ...
                        'fluxes binned to neutral density']);
    netcdf.putAtt(ncid,ndifID,'units','Watts/m^2');
    %    netcdf.putAtt(ncid,ndifID,'_FillValue',single(-1e20));
    netcdf.endDef(ncid);
else
    ndifID = netcdf.inqVarID(ncid,varname);
end


for ti=1:tL
    
    JI = zeros(xL,yL);
    QI = zeros(xL,yL);
    dVdt = zeros(xL,yL);
    dHdt = zeros(xL,yL);
    dift = zeros(xL,yL);
    netcdf.putVar(ncid,ndifID,[0 0 TL ti-1],[xL yL 1 1],zeros(xL,yL));

    for Ti=TL:-1:1
        sprintf('Calculating numdif time %03d of %03d, Temp %03d of %03d',ti,tL,Ti,TL)
        dVdt = dVdt + double(ncread(wname,'dVdt',[1 1 Ti ti],[xL yL 1 1]))*1e9/rho0./area;
        dHdt = dHdt + double(ncread(wname,'dHdt',[1 1 Ti ti],[xL yL 1 1]))./area;
        dift = dift + ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveSF)
            dift = dift + ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveKPP)
            dift = dift + ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveRedi)
            dift = dift + ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1])+...
                   ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveMDS)
            dift = dift + ncread(wname,'mixdownslope_temp_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end    
        if (haveSIG)
            dift = dift + ncread(wname,'temp_sigma_diff_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end    
        if (haveLT)
            dift = dift + ncread(wname,'temp_h_diffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
    
        txtrans = ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        tytrans = ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*tsc/rho0;
        % NOTE: THE tx_trans_submeso contribution SHOULD NOT be
        % added if it is implemented via skew-diffusion.
        % NOTE: THE tx_trans_gm contribution SHOULD NOT be
        % added if it is implemented via skew-diffusion.

        qxtrans = ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        qytrans = ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1]);

        % Submesoscale and GM: two options:
        if (haveSUB)
            if (doSGMviac)
                dift = dift + ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            else
                qxtrans = qxtrans + ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                qytrans = qytrans + ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
        end
        if (haveGM)
            if (doSGMviac)
                dift = dift + ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
            else
                qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
        end

        JI(2:end,2:end) = JI(2:end,2:end) + (txtrans(1:(end-1),2:end) - txtrans(2:end,2:end) ...
                                             +tytrans(2:end,1:(end-1)) - tytrans(2:end,2:end))./area(2:end,2:end);
        JI(1,2:end) = JI(1,2:end) + (- txtrans(1,2:end) ...
                                     +tytrans(1,1:(end-1)) - tytrans(1,2:end))./area(1,2:end);        
        JI(2:end,1) = JI(2:end,1) + (txtrans(1:(end-1),1) - txtrans(2:end,1) ...
                       - tytrans(2:end,1))./area(2:end,1);        
        JI(1,1) = JI(1,1) + (-txtrans(1,1)-tytrans(1,1))./area(1,1);
        
        QI(2:end,2:end) = QI(2:end,2:end) + (qxtrans(1:(end-1),2:end) - qxtrans(2:end,2:end) ...
                                             +qytrans(2:end,1:(end-1)) - qytrans(2:end,2:end))./area(2:end,2:end);
        QI(1,2:end) = QI(1,2:end) + (- qxtrans(1,2:end) ...
                                     +qytrans(1,1:(end-1)) - qytrans(1,2:end))./area(1,2:end);        
        QI(2:end,1) = QI(2:end,1) + (qxtrans(1:(end-1),1) - qxtrans(2:end,1) ...
                       - qytrans(2:end,1))./area(2:end,1);        
        QI(1,1) = QI(1,1) + (-qxtrans(1,1)-qytrans(1,1))./area(1,1);

        ndif = dHdt - (dVdt - JI)*rho0*Cp*Te(Ti) - dift - QI;
    
        netcdf.putVar(ncid,ndifID,[0 0 Ti-1 ti-1],[xL yL 1 1],ndif);
    end
end
netcdf.close(ncid);

end


%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------
if (doGWB)
    
GWB.dVdt   = zeros(TL+1,tL); % Sv 
try
    GWB.TENMON = TENMON; % Tendency from monthly averages (from previous code block)
catch
    GWB.TENMON = GWB.dVdt;
end
GWB.dHdt   = zeros(TL+1,tL); % Total W
GWB.VDS    = zeros(TL+1,tL); % W due to vdiffuse_sbc.
GWB.VDF    = zeros(TL+1,tL); % W due to vdiffusion
if (haveSUB)
    GWB.SUB    = zeros(TL+1,tL); % W due to submesoscale.
end
if (haveRedi)
    GWB.K33    = zeros(TL+1,tL); % W due to K33
    GWB.RED    = zeros(TL+1,tL); % W due to Redi diffusion
end
if (haveGM)
    GWB.NGM    = zeros(TL+1,tL); % W due to GM
end
if (haveLT)
    GWB.LTD    = zeros(TL+1,tL); % W due to Lateral Tracer diffusion
end
GWB.ADV    = zeros(TL+1,tL); % W due to advection
GWB.TEN    = zeros(TL+1,tL); % W due to tendency

if (doHND)
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
    GWB.dVdt(ii,ti) = nansum(nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = nansum(nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveSF)
        GWB.VDS(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.VDF(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveKPP)
        GWB.VDF(ii,ti) = GWB.VDF(ii,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveSUB)
    GWB.SUB(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveRedi)
    GWB.K33(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveLT)
    GWB.LTD(ii,ti) = nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_h_diffuse_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
for ii=TL-1:-1:1
    sprintf('Calculating global water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    GWB.dVdt(ii,ti) = GWB.dVdt(ii+1,ti) + nansum(nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1),2);
    GWB.dHdt(ii,ti) = GWB.dHdt(ii+1,ti) + nansum(nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.TEN(ii,ti) = GWB.TEN(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.ADV(ii,ti) = GWB.ADV(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveSF)
        GWB.VDS(ii,ti) = GWB.VDS(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    GWB.VDF(ii,ti) = GWB.VDF(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveKPP)
        GWB.VDF(ii,ti) = GWB.VDF(ii,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveSUB)
    GWB.SUB(ii,ti) = GWB.SUB(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveRedi)
    GWB.K33(ii,ti) = GWB.K33(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    GWB.RED(ii,ti) = GWB.RED(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    GWB.NGM(ii,ti) = GWB.NGM(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveLT)
    GWB.LTD(ii,ti) = GWB.LTD(ii+1,ti) + nansum(nansum(tmaskREG.*area.*ncread(wname,'temp_h_diffuse_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
end
end
save([outD model sprintf('_output%03d',output) '_HBud.mat'],'GWB','-v7.3');

end

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
if (doXY)
Tls = [0:2:26];
Nremain = length(Tls);
Ti = TL;

FlM = zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
FlI = zeros(xL,yL,tL); % numerical mixing
if (haveSUB)
    FlSUB = zeros(xL,yL,tL);
end
if (haveGM)
    FlGM  = zeros(xL,yL,tL);
end
if (doXYall)
    FlF = zeros(xL,yL,tL); % surface forcing
    if (haveRedi)
        FlK = zeros(xL,yL,tL); % K33
        FlR = zeros(xL,yL,tL); % Redi
    end
    if (haveLT)
        FlL = zeros(xL,yL,tL); % Lateral diffusion
    end
end
 
while (Nremain > 0 & Ti >= 1)
    Tl = Te(Ti);

    for ti=1:tL
        sprintf(['Calculating water-mass heat budget time %03d of ' ...
                 '%03d, temp %2.2f, going down to %2.2f'],ti,tL,Te(Ti),min(Tls))
        FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveKPP)
           FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        FlI(:,:,ti) = ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        if (haveSUB)
            FlSUB(:,:,ti) = FlSUB(:,:,ti)+ncread(wname,'temp_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (haveGM)
            FlGM(:,:,ti) = FlGM(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
        end
        if (doXYall)
            if (haveRedi)
                FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
                FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 Ti ti],[xL yL 1 1]);
            end
            if (haveSF)
                FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
            if (haveLT)
                FlL(:,:,ti) = FlL(:,:,ti)+ncread(wname,'temp_h_diffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1]);
            end
        end
    end

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'FlM','FlI','Tl','-v7.3');
        if (haveSUB)
            save(name,'FlSUB','-append');
        end
        if (haveGM)
            save(name,'FlGM','-append');
        end
        if (doXYall)
            if (haveSF)
                save(name,'FlF','-append');
            end
            if (haveRedi)
                save(name,'FlK','FlR','-append');
            end
            if (haveLT)
                save(name,'FlL','-append');
            end
        end
        
        Nremain = Nremain-1;
    end

    Ti = Ti-1;
end

end

%% Save surface heat flux, wind stress, SST, meridional heat flux:
if (doSURF)
try
    shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tLoc]);
catch
    shflux = ncread(fname_month,'net_sfc_heating',[1 1 1],[xL yL tLoc]);
end    
SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tLoc]));
% $$$ taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
% $$$ tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);

save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST');%,'taux','tauy');
end

%% Save vertically-integrated heat transports:
if (doXYtran)
Tls = [10 12.5 15 20 34];
Nremain = length(Tls);
Ti = 1;

qxflux = zeros(xL,yL);
qyflux = zeros(xL,yL);

mxflux = zeros(xL,yL);
myflux = zeros(xL,yL);

while (Nremain > 0 & Ti <= TL)
    Tl = Te(Ti+1);

    qxtrans = zeros(xL,yL);
    qytrans = zeros(xL,yL);
    mxtrans = zeros(xL,yL);
    mytrans = zeros(xL,yL);
    for ti=1:tL
        sprintf(['Calculating vertically-integrated heat and volume fluxes time %03d of ' ...
                 '%03d, temp %2.2f, going up to %2.2f'],ti,tL,Te(Ti),max(Tls))

        qxtrans = qxtrans+ncread(wname,'temp_xflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
        qytrans = qytrans+ncread(wname,'temp_yflux_adv_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
        mxtrans = mxtrans+ncread(wname,'tx_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;
        mytrans = mytrans+ncread(wname,'ty_trans_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;

        % Submesoscale and GM: two options:
        if (haveSUB)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qytrans = qytrans + ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            mxtrans = mxtrans + ncread(wname,'tx_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;
            mytrans = mytrans + ncread(wname,'ty_trans_nrho_submeso',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;
        end
        if (haveGM)
            qxtrans = qxtrans + ncread(wname,'temp_xflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qytrans = qytrans + ncread(wname,'temp_yflux_gm_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qxtrans = qxtrans + ncread(wname,'temp_xflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            qytrans = qytrans + ncread(wname,'temp_yflux_ndiffuse_on_nrho',[1 1 Ti ti],[xL yL 1 1])*ndays(ti);
            mxtrans = mxtrans + ncread(wname,'tx_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;
            mytrans = mytrans + ncread(wname,'ty_trans_nrho_gm',[1 1 Ti ti],[xL yL 1 1])*ndays(ti)*tsc/rho0;
        end
    end
    qxflux = qxflux + qxtrans/sum(ndays);
    qyflux = qyflux + qytrans/sum(ndays);
    mxflux = mxflux + mxtrans/sum(ndays);
    myflux = myflux + mytrans/sum(ndays);

    % Save heat flux terms:
    [sp,ind] = min(abs(Tls-Tl));
    if (abs(sp) <= dT/4)
        name = [outD model sprintf('_output%03d',output) '_XYtrans_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']

        save(name,'mxflux','myflux','qxflux','qyflux','Tl','-v7.3');
        Nremain = Nremain-1;
    end

    Ti = Ti+1;
end
end

%% Zonally-averaged fluxes -------------------------------------------------------------
if (doZA)
ZA.F = zeros(yL,TL+1,tL); % Wdeg-1 due to F
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
    ZA.AHDR   = zeros(yL,TL+1,tL); % Meridional heat flux due to
                                   % Redi diffusion
end
if (haveLT)
    ZA.LTD    = zeros(yL,TL+1,tL); % W due to Lateral diffusion
    ZA.AHDLT   = zeros(yL,TL+1,tL); % Meridional heat flux due to
                                   % lateral diffusion
end
ZA.JS = zeros(yL,TL+1,tL); % m3s-1deg-1 due to surface volume flux
ZA.PSI = zeros(yL,TL+1,tL); % m3s-1 northward transport
ZA.AHD = zeros(yL,TL+1,tL); % W A direct using heat fluxes
if (haveSUB)
    ZA.SUB    = zeros(yL,TL+1,tL);
    ZA.PSISUB = zeros(yL,TL+1,tL);
    ZA.AHDSUB = zeros(yL,TL+1,tL);
end
if (haveGM)
    ZA.GM    = zeros(yL,TL+1,tL);
    ZA.PSIGM = zeros(yL,TL+1,tL);
    ZA.AHDGM = zeros(yL,TL+1,tL);
end
if (doHND)
    ZA.NUM   = zeros(yL,TL+1,tL);
end
if (haveMDS)
    ZA.MDS = zeros(yL,TL+1,tL); % W
end
if (haveSIG)
    ZA.SIG = zeros(yL,TL+1,tL); % W
end
% ignoring tri-polar for now.
for ti=1:tL
    if (doHND)
        for ii = 1:(TL+1)
            ZA.NUM(:,ii,ti) = nansum(tmaskREG.*area.*ncread(wname,'temp_numdiff_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
        end
    end

    for ii = 1:TL
    sprintf('Calculating zonally-averaged water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    ZA.F(:,ii+1,ti) = ZA.F(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    ZA.P(:,ii+1,ti) = ZA.P(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                    nansum(tmaskREG.*area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    ZA.M(:,ii+1,ti) = ZA.M(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.KPPNL(:,ii+1,ti) = ZA.KPPNL(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    if (haveMDS)
        ZA.MDS(:,ii+1,ti) = ZA.MDS(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'mixdownslope_temp_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    if (haveSIG)
        ZA.SIG(:,ii+1,ti) = ZA.SIG(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'temp_sigma_diff_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    ZA.dVdt(:,ii+1,ti) = ZA.dVdt(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'dVdt',[1 1 ii ti],[xL yL 1 1])*1e9/rho0,1)';
    ZA.dHdt(:,ii+1,ti) = ZA.dHdt(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'dHdt',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.PSI(:,ii+1,ti)  = ZA.PSI(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]))'*tsc/rho0;
    ZA.AHD(:,ii+1,ti)  = ZA.AHD(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_adv_on_nrho',[1 1 ii ti],[xL yL 1 1]))';
    if (haveSUB)
        ZA.SUB(:,ii+1,ti) = ZA.SUB(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
        ZA.PSISUB(:,ii+1,ti) = ZA.PSISUB(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho_submeso',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
        ZA.AHDSUB(:,ii+1,ti) = ZA.AHDSUB(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end
    if (haveGM) 
        ZA.GM(:,ii+1,ti) = ZA.GM(:,ii,ti) + nansum(tmaskREG.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1)';
        ZA.PSIGM(:,ii+1,ti) = ZA.PSIGM(:,ii,ti) + nansum(umaskREG.*ncread(wname,'ty_trans_nrho_gm',[1 1 ii ti],[xL yL 1 1]),1)'*tsc/rho0;
        ZA.AHDGM(:,ii+1,ti) = ZA.AHDGM(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_gm_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end

    ZA.SWH(:,ii+1,ti) = ZA.SWH(:,ii,ti) + nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    ZA.F(:,ii+1,ti) = ZA.F(:,ii+1,ti) + nansum(tmaskREG.*area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    if (haveMIX)
        ZA.Mkppiw(:,ii+1,ti) = ZA.Mkppiw(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppiw_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mkppish(:,ii+1,ti) = ZA.Mkppish(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppish_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mwave(:,ii+1,ti) = ZA.Mwave(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_wave_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Mkppbl(:,ii+1,ti) = ZA.Mkppbl(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppbl_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.Moth(:,ii+1,ti) = ZA.Moth(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppicon_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_diff_cbt_kppdd_on_nrho',[1 1 ii ti],[xL yL 1 1]),1) + ...
                            nansum(tmaskREG.*area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
    end
    if (haveRedi)
        ZA.K33(:,ii+1,ti) = ZA.K33(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.RED(:,ii+1,ti) = ZA.RED(:,ii,ti) + (nansum(tmaskREG.*area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1))';
        ZA.AHDR(:,ii+1,ti) = ZA.AHDR(:,ii,ti) + nansum(umaskREG.*ncread(wname,'temp_yflux_ndiffuse_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)';
    end        
    ZA.JS(:,ii+1,ti) = ZA.JS(:,ii,ti) + (nansum(tmaskREG.*ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1]),1)/rho0)';
    end
    save([outD model sprintf('_output%03d',output) '_ZAHBud.mat'],'ZA','yu','yt','-v7.3');
end
end % End doZA

% $$$ end

