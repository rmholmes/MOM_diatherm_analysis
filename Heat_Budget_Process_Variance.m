% This script extracts specified data from MOM025 and MOM01 runs 

% $$$ baseL = '/short/e14/rmh561/mom/archive/';
baseL = '/g/data/e14/rmh561/mom/archive/';
% $$$ baseL = '/short/e14/rmh561/access-om2/archive/';
% $$$ baseL = '/srv/ccrc/data03/z3500785/';
% $$$ types = {'kds50','gfdl50','kds75','kds100','kds135'};

model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/'];

outD = [baseD];

Ts = [10 15 22.5];
TL = length(Ts);

% $$$ output = 95;
for output=101:110

%% file-names and grid properties:
base = [baseD sprintf('output%03d/',output)];
fname = [base 'ocean.nc'];
gname = [base 'ocean_grid.nc'];
         
lon = ncread(gname,'geolon_t');lat = ncread(gname,'geolat_t');
area = ncread(gname,'area_t');[xL,yL] = size(lon);
lonu = ncread(gname,'geolon_c');latu = ncread(gname,'geolat_c');

z = ncread(fname,'st_ocean');zL = length(z);

time = ncread(fname,'time');
tL = length(time);

wvarA = NaN*zeros(xL,yL,tL,TL);
EKEA  = NaN*zeros(xL,yL,tL,TL);
TdxsqA  = NaN*zeros(xL,yL,tL,TL);
TdysqA  = NaN*zeros(xL,yL,tL,TL);
TdzsqA  = NaN*zeros(xL,yL,tL,TL);

for ti=1:tL

    % Load Eulerian variables:
    temp = ncread(fname,'temp',[1 1 1 ti],[xL yL zL 1]);
    
    EKEE  = ncread(fname,'u_sq',[1 1 1 ti],[xL yL zL 1])-ncread(fname,'u',[1 1 1 ti],[xL yL zL 1]).^2 + ...
            ncread(fname,'v_sq',[1 1 1 ti],[xL yL zL 1])-ncread(fname,'v',[1 1 1 ti],[xL yL zL 1]).^2; %u-pts
    TdxsqE  = ncread(fname,'temp_dxsq',[1 1 1 ti],[xL yL zL 1]); %flux-x
    TdysqE  = ncread(fname,'temp_dysq',[1 1 1 ti],[xL yL zL 1]); %flux-y

    wvarE = ncread(fname,'wt_sq',[1 1 1 ti],[xL yL zL 1])-ncread(fname,'wt',[1 1 1 ti],[xL yL zL 1]).^2;%wt
    TdzsqE  = ncread(fname,'temp_dzsq',[1 1 1 ti],[xL yL zL 1]);%wt
    
    wvarE = cat(3,zeros(xL,yL,1),wvarE);
    wvarE = (wvarE(:,:,2:end)+wvarE(:,:,1:(end-1)))/2; %t-pts
    
    TdzsqE = cat(3,zeros(xL,yL,1),TdzsqE);
    TdzsqE = (TdzsqE(:,:,2:end)+TdzsqE(:,:,1:(end-1)))/2; %t-pts

    % Do 1D interpolations:
    for xi = 1:xL-1
        sprintf('Doing time %03d of %03d, x %03d of %03d',ti,tL,xi,xL)
        for yi = 1:yL-1
                
            tvecU = squeeze(temp(xi,yi,:)+temp(xi+1,yi,:)+temp(xi,yi+1,:)+temp(xi+1,yi+1,:))/4;
                tvecX = squeeze(temp(xi,yi,:)+temp(xi+1,yi,:))/2;
                tvecY = squeeze(temp(xi,yi,:)+temp(xi,yi+1,:))/2;
            tvecT = squeeze(temp(xi,yi,:));
            for Ti = 1:TL

                ziU = find(tvecU<Ts(Ti),1,'first');
                ziX = find(tvecX<Ts(Ti),1,'first');
                ziY = find(tvecY<Ts(Ti),1,'first');
                ziT = find(tvecT<Ts(Ti),1,'first');

                if (length(ziU)>0)
                    if (ziU>1)
% $$$                         EKEA(xi,yi,ti,Ti) = interp1(tvecU(ziU-1:ziU),squeeze(EKEE(xi,yi,ziU-1:ziU)),Ts(Ti),'linear');
                        EKEA(xi,yi,ti,Ti) = (EKEE(xi,yi,ziU)-EKEE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Ts(Ti)-tvecU(ziU-1))+EKEE(xi,yi,ziU-1);
                    end
                end

                if (length(ziX)>0)
                    if (ziX>1)
% $$$                         TdxsqA(xi,yi,ti,Ti) = interp1(tvecX(ziX-1:ziX),squeeze(TdxsqE(xi,yi,ziX-1:ziX)),Ts(Ti),'linear');
                        TdxsqA(xi,yi,ti,Ti) = (TdxsqE(xi,yi,ziX)-TdxsqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Ts(Ti)-tvecX(ziX-1))+TdxsqE(xi,yi,ziX-1);
                    end
                end

                if (length(ziY)>0)
                    if (ziY>1)
% $$$                         TdysqA(xi,yi,ti,Ti) = interp1(tvecY(ziY-1:ziY),squeeze(TdysqE(xi,yi,ziY-1:ziY)),Ts(Ti),'linear');
                        TdysqA(xi,yi,ti,Ti) = (TdysqE(xi,yi,ziY)-TdysqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Ts(Ti)-tvecY(ziY-1))+TdysqE(xi,yi,ziY-1);
                    end
                end

                if (length(ziT)>0)
                    if (ziT>1)
% $$$                         wvarA(xi,yi,ti,Ti) = interp1(tvecT(ziT-1:ziT),squeeze(wvarE(xi,yi,ziT-1:ziT)),Ts(Ti),'linear');
                        wvarA(xi,yi,ti,Ti) = (wvarE(xi,yi,ziT)-wvarE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Ts(Ti)-tvecT(ziT-1))+wvarE(xi,yi,ziT-1);
                        TdzsqA(xi,yi,ti,Ti) = (TdzsqE(xi,yi,ziT)-TdzsqE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Ts(Ti)-tvecT(ziT-1))+TdzsqE(xi,yi,ziT-1);
                    end
                end
            end
        end
    end
end

for Ti=1:TL
    wvar = wvarA(:,:,:,Ti);
    EKE = EKEA(:,:,:,Ti);
    Tdxsq = TdxsqA(:,:,:,Ti);
    Tdysq = TdysqA(:,:,:,Ti);
    Tdzsq = TdzsqA(:,:,:,Ti);
    
    name = [outD 'mat_data/' model sprintf('_output%03d',output) '_variances_T' strrep(num2str(Ts(Ti)),'.','p') 'C.mat']
    save(name,'lon','lat','wvar','EKE','Tdxsq','Tdysq','Tdzsq','time');
end


end
