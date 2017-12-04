% This script processes the heat budget and associated variables in
% MOM025 or MOM01 simulations and save's into .mat files

baseD = '/short/e14/rmh561/mom/archive/MOM_wombat/'; %Data Directory.
model = 'MOM025';
rstbaseD = baseD;%'/short/e14/rmh561/mom/archive/MOM_HeatDiag/'; %Data
outD = '/short/e14/rmh561/mom/archive/MOM_wombat/mat_data/'; %Data

haveRedi = 1; % 1 = Redi diffusion is on, 0 = off
haveGM = 1; % 1 = GM is on, 0 = off;

output=1978
restart = output-1;

% file-names -----------------------------------------
base = [baseD sprintf('output%03d/',output)];
basem1 = [baseD sprintf('output%03d/',output-1)];
baser = [rstbaseD sprintf('restart%03d/',restart)];
hname = [base 'ocean_heat.nc'];
if (strfind(baseD,'01'))
    fname = [base 'ocean_month.nc'];
else
    fname = [base 'ocean.nc'];
end
gname = [base 'ocean_grid.nc'];
sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
if (exist(baser))
    found_rst = 1;rstti = 1;
    rnameT = [baser 'ocean_temp_salt.res.nc'];
    rnameZ = [baser 'ocean_thickness.res.nc'];
    rnametime = [baser 'coupler.res'];
else
    found_rst = 0;rstti = 12;
    rnameT = [basem1 'ocean_snap.nc'];
    rnameZ = [basem1 'ocean_snap.nc'];
    rnametime = [basem1 'ocean_snap.nc'];
end    
         
% Horizontal Grid  -----------------------------------------
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonv_t = ncread(gname,'xt_ocean');lonv_u = ncread(gname,'xu_ocean');
latv_t = ncread(gname,'yt_ocean');latv_u = ncread(gname,'yu_ocean');

% Vertical grid  -----------------------------------------
z = ncread(hname,'st_ocean');zL = length(z);

% Time  -----------------------------------------
time = ncread(hname,'time');

if (found_rst)
    dys = [31 28 31 30 31 30 31 31 30 31 30 31];
    C = textread(rnametime, '%s','delimiter', '\n');
    C = strsplit(C{3});
    rtime = [str2num(C{1}) str2num(C{2}) str2num(C{3}) str2num(C{4}) str2num(C{5}) str2num(C{6})];
    time_snap = [(rtime(1)-1)*365+sum(dys(1:(rtime(2)-1)))+(rtime(3)-1)+rtime(4)/24+rtime(5)/24/60+rtime(6)/24/60/60;
                 ncread(sname,'time')];
else
    time_snapl = ncread(rnametime,'time');
    time_snap = [time_snapl(end); ncread(sname,'time')];
end

time_snap = mod(time_snap,365);
if (time_snap(end) == 0) time_snap(end) = 365;end

tL = length(time);

Cp = 3992.1; % J kg-1 degC-1
rho0 = 1035; % kgm-3

T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

%latitude vector for heat function:
latv = max(lat,[],1);
late = [-90 (latv(2:end)+latv(1:(end-1)))/2 90];

% $$$ save([outD model sprintf('_output%03d',output) '_BaseVars.mat'], ...
% $$$      'T','Te','TL','dT','Cp','rho0','time','time_snap','tL', ...
% $$$      'z','zL','lon','lat','area','xL','yL','latv','late', ...
% $$$      'lonu','latu');

%% Calculate volume integrated budget from online T-binned values -----------------------------------------------------------------------------------------------------------
V      = zeros(TL+1,tL); % Volume of water (m3) above temperature T
H      = zeros(TL+1,tL); % Heat content (J) above temperature T
Vsnap  = zeros(TL+1,tL+1); % Volume from snapshots (m3)
Hsnap  = zeros(TL+1,tL+1); % Heat content from snapshots (J)
Temp   = zeros(zL,tL); % Temperature as a function of depth (degC)
SWH    = zeros(TL+1,tL); % W due to SW redistribution
VDS    = zeros(TL+1,tL); % W due to vdiffuse_sbc.
RMX    = zeros(TL+1,tL); % W due to rivermix.
PME    = zeros(TL+1,tL); % W due to P-E.
FRZ    = zeros(TL+1,tL); % W due to frazil.
ETS    = zeros(TL+1,tL); % W due to eta_smoothing.
SUB    = zeros(TL+1,tL); % W due to submesoscale.
VDF    = zeros(TL+1,tL); % W due to vdiffusion
KNL    = zeros(TL+1,tL); % W due to KPP non-local
if (haveRedi)
    K33    = zeros(TL+1,tL); % W due to K33
    RED    = zeros(TL+1,tL); % W due to Redi diffusion
end
if (haveGM)
    NGM    = zeros(TL+1,tL); % W due to GM
end
ADV    = zeros(TL+1,tL); % W due to advection
TEN    = zeros(TL+1,tL); % W due to tendency
SFW    = zeros(TL+1,tL); % surface volume flux into ocean (m3s-1)
TENMON = zeros(TL+1,tL); % W due to tendency from Offline Monthly

%Do IC for Vsnap and Hsnap:
for zi = 1:zL
    sprintf('Doing snapshot IC, depth %02d of %02d',zi,zL)
    %Temperature snapshot:
    tempsnap = ncread(rnameT,'temp',[1 1 zi rstti],[xL yL 1 1]);
    tempsnap(tempsnap==0) = NaN; %This is included because the
                                 %restarts don't have any NaNs in
                                 %temp, just lots of 0s. 
    if (found_rst)
        Volsnap = ncread(rnameZ,'rho_dzt',[1 1 zi rstti],[xL yL 1 1]).*area/rho0;
    else
        Volsnap = ncread(rnameT,'dzt',[1 1 zi rstti],[xL yL 1 1]).*area;
    end
    %Accumulate sums:
    for Ti=1:TL
        inds = find(tempsnap>=Te(Ti) & tempsnap<Te(Ti+1));
        Vsnap(Ti,1) = Vsnap(Ti,1)+nansum(Volsnap(inds));
        Hsnap(Ti,1) = Hsnap(Ti,1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
    end
    inds = find(tempsnap>=Te(TL+1));
    Vsnap(TL+1,1) = Vsnap(TL+1,1)+nansum(Volsnap(inds));
    Hsnap(TL+1,1) = Hsnap(TL+1,1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
end
%Integrate to get to T'>T:
Vsnap(:,1) = flipud(cumsum(flipud(Vsnap(:,1))));
Hsnap(:,1) = flipud(cumsum(flipud(Hsnap(:,1))));

%Do Eulerian budget calculations:
for ti=1:tL
    for zi=1:zL
        sprintf('Doing Eul Bud. time %03d of %03d, depth %02d of %02d',ti,tL,zi,zL)

        temp = ncread(fname,'temp',[1 1 zi ti],[xL yL 1 1]);
        tempsnap = ncread(sname,'temp',[1 1 zi ti],[xL yL 1 1]);
        Vol = ncread(fname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;
        Volsnap = ncread(sname,'dzt',[1 1 zi ti],[xL yL 1 1]).*area;

        %Calculate T(z):
        areaNaN = area;
        areaNaN(isnan(temp)) = NaN;
        Temp(zi,ti) = squeeze(nansum(nansum(temp.*area,1),2)./nansum(nansum(areaNaN,1),2));

        %Tendency from Monthly snapshots:
        TENf = area.*ncread(hname,'temp_tendency',[1 1 zi ti],[xL yL 1 1]);

        %Accumulate sums:
        for Ti=1:TL
            inds = find(tempsnap>=Te(Ti) & tempsnap<Te(Ti+1));
            Vsnap(Ti,ti+1) = Vsnap(Ti,ti+1)+nansum(Volsnap(inds));
            Hsnap(Ti,ti+1) = Hsnap(Ti,ti+1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
            inds = find(temp>=Te(Ti) & temp<Te(Ti+1));
            V(Ti,ti) = V(Ti,ti)+nansum(Vol(inds));
            H(Ti,ti) = H(Ti,ti)+nansum(Vol(inds).*temp(inds)*rho0*Cp);
            TENMON(Ti,ti) = TENMON(Ti,ti)+nansum(TENf(inds));            
        end
        inds = find(tempsnap>=Te(TL+1));
        Vsnap(TL+1,ti+1) = Vsnap(TL+1,ti+1)+nansum(Volsnap(inds));
        Hsnap(TL+1,ti+1) = Hsnap(TL+1,ti+1)+nansum(Volsnap(inds).*tempsnap(inds)*rho0*Cp);
        inds = find(temp>=Te(TL+1));
        V(TL+1,ti) = V(TL+1,ti)+nansum(Vol(inds));
        H(TL+1,ti) = H(TL+1,ti)+nansum(Vol(inds).*temp(inds)*rho0*Cp);
        TENMON(TL+1,ti) = TENMON(TL+1,ti)+nansum(TENf(inds));            
    end
    %Integrate to get to T'>T:
    Vsnap(:,ti+1) = flipud(cumsum(flipud(Vsnap(:,ti+1))));
    V(:,ti) = flipud(cumsum(flipud(V(:,ti))));
    Hsnap(:,ti+1) = flipud(cumsum(flipud(Hsnap(:,ti+1))));
    H(:,ti) = flipud(cumsum(flipud(H(:,ti))));
    TENMON(:,ti) = flipud(cumsum(flipud(TENMON(:,ti))));
end

save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'], ...
     'Vsnap','Hsnap','V','H','Temp','TENMON');

for ti=1:tL
    ii = TL;
    sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    TEN(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ADV(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SUB(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    PME(ii,ti) = nansum(nansum(area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RMX(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDS(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SWH(ii,ti) = nansum(nansum(area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDF(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    KNL(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    K33(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RED(ii,ti) = nansum(nansum(area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    NGM(ii,ti) = nansum(nansum(area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    FRZ(ii,ti) = nansum(nansum(area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ETS(ii,ti) = nansum(nansum(area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SFW(ii,ti) = nansum(nansum(ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
for ii=TL-1:-1:1
    sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,ii,TL)
    TEN(ii,ti) = TEN(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ADV(ii,ti) = ADV(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SUB(ii,ti) = SUB(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    PME(ii,ti) = PME(ii+1,ti) + nansum(nansum(area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RMX(ii,ti) = RMX(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDS(ii,ti) = VDS(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SWH(ii,ti) = SWH(ii+1,ti) + nansum(nansum(area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    VDF(ii,ti) = VDF(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    KNL(ii,ti) = KNL(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    if (haveRedi)
    K33(ii,ti) = K33(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    RED(ii,ti) = RED(ii+1,ti) + nansum(nansum(area.*ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    if (haveGM)
    NGM(ii,ti) = NGM(ii+1,ti) + nansum(nansum(area.*ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]),1),2);
    end
    FRZ(ii,ti) = FRZ(ii+1,ti) + nansum(nansum(area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    ETS(ii,ti) = ETS(ii+1,ti) + nansum(nansum(area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]),1),2);
    SFW(ii,ti) = SFW(ii+1,ti) + nansum(nansum(ncread(wname,'mass_pmepr_on_nrho',[1 1 ii ti],[xL yL 1 1])/rho0,1),2);
end
end

save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'],'SWH','VDS','RMX','PME','FRZ', ...
     'ETS','SUB','VDF','KNL','ADV','TEN','SFW','-append');
if (haveRedi)
save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'],'K33','RED','-append');
end
if (haveGM)
    save([outD model sprintf('_output%03d',output) '_GlobalHBud.mat'],'NGM','-append');
end

%% Vertical Integrate down to level from online T-binned values -----------------------------------------------------------------------------------------------------------
Tls = [5 10 15 17.5 20 22.5 25 27.5];

for ii = 1:length(Tls)
    Tl = Tls(ii);

% $$$     FlM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
% $$$     FlF = NaN*zeros(xL,yL,tL); % surface forcing
% $$$     FlP = NaN*zeros(xL,yL,tL); % P-E+R
% $$$     FlA = NaN*zeros(xL,yL,tL); % advection + submeso + GM
% $$$     if (haveRedi)
% $$$         FlK = NaN*zeros(xL,yL,tL); % K33
% $$$         FlR = NaN*zeros(xL,yL,tL); % Redi
% $$$     end
    if (haveGM)
        FlG = NaN*zeros(xL,yL,tL); % GM
    end
% $$$     FlT = NaN*zeros(xL,yL,tL); % tendency
% $$$     FlSP = NaN*zeros(xL,yL,tL); % solar penetration
    T = ncread(wname,'neutral');
    Te = ncread(wname,'neutralrho_edges');
    [tmp Ti] = min(abs(Te-Tl));

    for ti=1:tL
        ii = TL;
        sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,TL-ii+1,TL-Ti+1)
% $$$     FlT(:,:,ti) = ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$         if (haveRedi)
% $$$             FlK(:,:,ti) = ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$             FlR(:,:,ti) = ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]);
% $$$         end
% $$$     FlA(:,:,ti) = ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
        if (haveGM)
            FlG(:,:,ti) = ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]);
% $$$             FlA(:,:,ti) = FlA(:,:,ti) + FlG(:,:,ti);
        end
% $$$     FlP(:,:,ti) = ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$         FlM(:,:,ti) = ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$             ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlF(:,:,ti) = ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$         FlSP(:,:,ti) = ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]);
    
        for ii=TL-1:-1:Ti
            sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,TL-ii+1,TL-Ti+1)
% $$$     FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$             if (haveRedi)
% $$$                 FlK(:,:,ti) = FlK(:,:,ti)+ncread(wname,'temp_vdiffuse_k33_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$                 FlR(:,:,ti) = FlR(:,:,ti)+ncread(wname,'neutral_diffusion_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]);
% $$$             end
% $$$     FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
            if (haveGM)
                FlG(:,:,ti) = FlG(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]);
% $$$                 FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'neutral_gm_on_nrho_temp',[1 1 ii ti],[xL yL 1 1]);
            end
% $$$     FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$             FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                 ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$             FlSP(:,:,ti) = FlSP(:,:,ti)+ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]);
        end
    end

% $$$     save([outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tl),'.','p') 'C.mat'],'FlM','FlSP','Tl');
% $$$     if (haveRedi)
% $$$         save([outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tl),'.','p') 'C.mat'],'FlK','FlR','-append');
% $$$     end
    if (haveGM)
        save([outD model sprintf('_output%03d',output) '_VertInt_T' strrep(num2str(Tl),'.','p') 'C.mat'],'FlG','-append');
    end
end

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
% $$$ %% Save surface heat flux, wind stress and SST:
% $$$ shflux = ncread(fname,'net_sfc_heating',[1 1 1],[xL yL tL]);
% $$$ SST = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
% $$$ taux = ncread(fname,'tau_x',[1 1 1],[xL yL tL]);
% $$$ tauy = ncread(fname,'tau_y',[1 1 1],[xL yL tL]);
% $$$ 
% $$$ save([outD model sprintf('_output%03d',output) '_SurfaceVars.mat'],'shflux','SST','taux','tauy');
% $$$ 
% $$$ %% Latitude-depth overturning:
% $$$ PSI = zeros(yL,TL,tL);
% $$$ 
% $$$ % Ignoring tri-polar (wrong >60N).
% $$$ for ti=1:tL
% $$$     ti
% $$$     PSI(:,:,ti) = squeeze(nansum(ncread(wname,'ty_trans_nrho',[1 1 1 ti],[xL yL TL 1]),1));
% $$$ end
% $$$ save([outD model sprintf('_output%03d',output) '_Tpsi.mat'],'PSI');

% $$$     
% $$$     for TI = 1:TL
% $$$         tytrans = ncread(wname,'ty_trans_nrho',[1 1 TI ti],[xL yL 1 1]);
% $$$         for yi=1:yL
% $$$             inds = lat >= late(yi) & lat <= 
% $$$             PSI(yi,TI,ti) = sum(tytrans(


% $$$ %% Heat Function From tendencies -------------------------------------------------
% $$$ HFSWH    = zeros(yL,TL+1,tL); % W due to SW redistribution
% $$$ HFVDS    = zeros(yL,TL+1,tL); % W due to vdiffuse_sbc.
% $$$ HFRMX    = zeros(yL,TL+1,tL); % W due to rivermix.
% $$$ HFPME    = zeros(yL,TL+1,tL); % W due to P-E.
% $$$ HFFRZ    = zeros(yL,TL+1,tL); % W due to frazil.
% $$$ HFETS    = zeros(yL,TL+1,tL); % W due to eta_smoothing.
% $$$ HFSUB    = zeros(yL,TL+1,tL); % W due to submesoscale.
% $$$ HFVDF    = zeros(yL,TL+1,tL); % W due to vdiffusion
% $$$ HFKNL    = zeros(yL,TL+1,tL); % W due to KPP non-local
% $$$ HFADV    = zeros(yL,TL+1,tL); % W due to advection
% $$$ HFTEN    = zeros(yL,TL+1,tL); % W due to tendency
% $$$ 
% $$$ tph = zeros(yL,1);
% $$$ for yi=1:yL
% $$$     tph(yi) = all(lat(:,yi) == latv(yi));
% $$$ end
% $$$ tph = ~tph;
% $$$ tplast = find(tph==1,1,'first');
% $$$ 
% $$$ for ti=1:tL
% $$$     for ii=1:TL
% $$$         TEN = area.*ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);TEN(isnan(TEN)) = 0;
% $$$         ADV = area.*ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1]);ADV(isnan(ADV)) = 0;
% $$$         SUB = area.*ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);SUB(isnan(SUB)) = 0;
% $$$         PME = area.*ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1]);PME(isnan(PME)) = 0;
% $$$         RMX = area.*ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);RMX(isnan(RMX)) = 0;
% $$$         VDS = area.*ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1]);VDS(isnan(VDS)) = 0;
% $$$         SWH = area.*ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1]);SWH(isnan(SWH)) = 0;
% $$$         VDF = area.*ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1]);VDF(isnan(VDF)) = 0;
% $$$         KNL = area.*ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);KNL(isnan(KNL)) = 0;
% $$$         FRZ = area.*ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1]);FRZ(isnan(FRZ)) = 0;
% $$$         ETS = area.*ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);ETS(isnan(ETS)) = 0;
% $$$         
% $$$         % Tripolar region:
% $$$         for yi=yL:-1:(tplast-1)
% $$$             sprintf('Calculating heat function lat %03d of %03d, time %03d of %03d, temp %03d of %03d',yL-yi+1,yL,ti,tL,ii,TL)
% $$$             inds = lat >= latv(yi);
% $$$             HFTEN(yi,ii+1,ti) = HFTEN(yi,ii,ti) + sum(TEN(inds));
% $$$             HFADV(yi,ii+1,ti) = HFADV(yi,ii,ti) + sum(ADV(inds));
% $$$             HFSUB(yi,ii+1,ti) = HFSUB(yi,ii,ti) + sum(SUB(inds));
% $$$             HFPME(yi,ii+1,ti) = HFPME(yi,ii,ti) + sum(PME(inds));
% $$$             HFRMX(yi,ii+1,ti) = HFRMX(yi,ii,ti) + sum(RMX(inds));
% $$$             HFVDS(yi,ii+1,ti) = HFVDS(yi,ii,ti) + sum(VDS(inds));
% $$$             HFSWH(yi,ii+1,ti) = HFSWH(yi,ii,ti) + sum(SWH(inds));
% $$$             HFVDF(yi,ii+1,ti) = HFVDF(yi,ii,ti) + sum(VDF(inds));
% $$$             HFKNL(yi,ii+1,ti) = HFKNL(yi,ii,ti) + sum(KNL(inds));
% $$$             HFFRZ(yi,ii+1,ti) = HFFRZ(yi,ii,ti) + sum(FRZ(inds));
% $$$             HFETS(yi,ii+1,ti) = HFETS(yi,ii,ti) + sum(ETS(inds));
% $$$         end
% $$$         %Non-tripolar region:
% $$$         HFTEN(1:(tplast-2),ii+1,ti) = HFTEN(1:(tplast-2),ii,ti) + cumsum(sum(TEN(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(TEN(:,(tplast-1):end),1),2);
% $$$         HFADV(1:(tplast-2),ii+1,ti) = HFADV(1:(tplast-2),ii,ti) + cumsum(sum(ADV(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(ADV(:,(tplast-1):end),1),2);
% $$$         HFSUB(1:(tplast-2),ii+1,ti) = HFSUB(1:(tplast-2),ii,ti) + cumsum(sum(SUB(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(SUB(:,(tplast-1):end),1),2);
% $$$         HFPME(1:(tplast-2),ii+1,ti) = HFPME(1:(tplast-2),ii,ti) + cumsum(sum(PME(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(PME(:,(tplast-1):end),1),2);
% $$$         HFRMX(1:(tplast-2),ii+1,ti) = HFRMX(1:(tplast-2),ii,ti) + cumsum(sum(RMX(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(RMX(:,(tplast-1):end),1),2);
% $$$         HFVDS(1:(tplast-2),ii+1,ti) = HFVDS(1:(tplast-2),ii,ti) + cumsum(sum(VDS(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(VDS(:,(tplast-1):end),1),2);
% $$$         HFSWH(1:(tplast-2),ii+1,ti) = HFSWH(1:(tplast-2),ii,ti) + cumsum(sum(SWH(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(SWH(:,(tplast-1):end),1),2);
% $$$         HFVDF(1:(tplast-2),ii+1,ti) = HFVDF(1:(tplast-2),ii,ti) + cumsum(sum(VDF(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(VDF(:,(tplast-1):end),1),2);
% $$$         HFKNL(1:(tplast-2),ii+1,ti) = HFKNL(1:(tplast-2),ii,ti) + cumsum(sum(KNL(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(KNL(:,(tplast-1):end),1),2);
% $$$         HFFRZ(1:(tplast-2),ii+1,ti) = HFFRZ(1:(tplast-2),ii,ti) + cumsum(sum(FRZ(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(FRZ(:,(tplast-1):end),1),2);
% $$$         HFETS(1:(tplast-2),ii+1,ti) = HFETS(1:(tplast-2),ii,ti) + cumsum(sum(ETS(:,1:(tplast-2)),1),2,'reverse')' + sum(sum(ETS(:,(tplast-1):end),1),2);
% $$$ % $$$ 
% $$$ % $$$         for yi=(tplast-2):-1:1
% $$$ % $$$             sprintf('Calculating heat function lat %03d of %03d, time %03d of %03d, temp %03d of %03d',yL-yi+1,yL,ti,tL,ii,TL)
% $$$ % $$$             HFTEN(yi,ii+1,ti) = HFTEN(yi,ii,ti) + nansum(nansum(TEN(:,yi:end),1),2);
% $$$ % $$$             HFADV(yi,ii+1,ti) = HFADV(yi,ii,ti) + nansum(nansum(ADV(:,yi:end),1),2);
% $$$ % $$$             HFSUB(yi,ii+1,ti) = HFSUB(yi,ii,ti) + nansum(nansum(SUB(:,yi:end),1),2);
% $$$ % $$$             HFPME(yi,ii+1,ti) = HFPME(yi,ii,ti) + nansum(nansum(PME(:,yi:end),1),2);
% $$$ % $$$             HFRMX(yi,ii+1,ti) = HFRMX(yi,ii,ti) + nansum(nansum(RMX(:,yi:end),1),2);
% $$$ % $$$             HFVDS(yi,ii+1,ti) = HFVDS(yi,ii,ti) + nansum(nansum(VDS(:,yi:end),1),2);
% $$$ % $$$             HFSWH(yi,ii+1,ti) = HFSWH(yi,ii,ti) + nansum(nansum(SWH(:,yi:end),1),2);
% $$$ % $$$             HFVDF(yi,ii+1,ti) = HFVDF(yi,ii,ti) + nansum(nansum(VDF(:,yi:end),1),2);
% $$$ % $$$             HFKNL(yi,ii+1,ti) = HFKNL(yi,ii,ti) + nansum(nansum(KNL(:,yi:end),1),2);
% $$$ % $$$             HFFRZ(yi,ii+1,ti) = HFFRZ(yi,ii,ti) + nansum(nansum(FRZ(:,yi:end),1),2);
% $$$ % $$$             HFETS(yi,ii+1,ti) = HFETS(yi,ii,ti) + nansum(nansum(ETS(:,yi:end),1),2);
% $$$ % $$$         end
% $$$     end
% $$$ end
% $$$ save([outD model sprintf('_output%03d',output) '_HFunc.mat'],'HFSWH','HFVDS','HFRMX','HFPME','HFFRZ', ...
% $$$      'HFETS','HFSUB','HFVDF','HFKNL','HFADV','HFTEN');

% $$$ end
