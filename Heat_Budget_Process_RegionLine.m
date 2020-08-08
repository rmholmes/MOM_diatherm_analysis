% This script processes a particular region for num-mix and
% variance variance terms as a function of temperature.

region = [-215 -180 30 45];
regname = 'Kuroshio';

baseL = '/scratch/e14/rmh561/access-om2/archive/';
baseL = '/g/data/e14/rmh561/access-om2/archive/';

model = 'ACCESS-OM2_025deg_jra55_ryf_norediGM';
baseD = [baseL '025deg_jra55_ryf_norediGM/']; %Data Directory.

outD = [baseD 'mat_data/'];

post = 'ocean/';

%%% Do num-mix curves:
outputs = [76 80];

output = outputs(1);
load([outD model sprintf('_output%03d',output) '_BaseVars.mat']);

[tmp x1] = min(abs(xt-region(1)));
[tmp x2] = min(abs(xt-region(2)));
[tmp y1] = min(abs(xt-region(3)));
[tmp y2] = min(abs(xt-region(4)));

NUM = zeros(TL+1,tL,length(outputs));
for oi = 1:length(outputs)
    output = outputs(oi);
    base = [baseD sprintf('output%03d/',output) post];
    wname = [base 'ocean_grid.nc'];

    for ti=1:tL
        sprintf('Num-mix out %03d of %03d, time %03d of %03d',oi,length(outputs),ti,tL)
        for Ti = (TL+1):-1:1
            NUM(Ti,ti,oi) = nansum(nansum(area(x1:x2,y1:y2).*ncread(wname, ...
                         'temp_numdiff_heat_on_nrho',[x1 y1 Ti ti],[x2-x1+1 y2-y1+1 1 1]),1),2);
        end
    end
end

%%% Do variances:
outputs = [81];
output = outputs(1);
base = [baseD sprintf('output%03d/',output) post];
fname = [base 'ocean_month.nc'];
gname = [base 'ocean_grid.nc'];
         
dxt = ncread(gname,'dxt',[x1 y1],[x2-x1+1 y2-y1+1]);
dyt = ncread(gname,'dyt',[x1 y1],[x2-x1+1 y2-y1+1]);
dxu = ncread(gname,'dxu',[x1 y1],[x2-x1+1 y2-y1+1]);
dyu = ncread(gname,'dyu',[x1 y1],[x2-x1+1 y2-y1+1]);
area = ncread(gname,'area_t',[x1 y1],[x2-x1+1 y2-y1+1]);
area_T = nansum(nansum(area));
arear = area/area_T;

time = ncread(fname,'time');
ndays = ncread(fname,'average_DT');
tL = length(time);

% $$$ wvarA = NaN*zeros(tL,TL,length(outputs));
EKE  = NaN*zeros(tL,TL,length(outputs));
Tdxsq  = NaN*zeros(tL,TL,length(outputs));
Tdysq  = NaN*zeros(tL,TL,length(outputs));
Tdzsq  = NaN*zeros(tL,TL,length(outputs));
aiso = NaN*zeros(tL,TL,length(outputs));
rey_bih = NaN*zeros(tL,TL,length(outputs));
udxsq  = NaN*zeros(tL,TL,length(outputs));
udysq  = NaN*zeros(tL,TL,length(outputs));
vdxsq  = NaN*zeros(tL,TL,length(outputs));
vdysq  = NaN*zeros(tL,TL,length(outputs));

for oi = 1:length(outputs)
    output = outputs(oi);
    base = [baseD sprintf('output%03d/',output) post];
    fname = [base 'ocean_month.nc'];
for ti=1:tL
    ['Doing time ' num2str(ti) ' of ' num2str(tL)]
    
    % Load Eulerian variables:
    temp = ncread(fname,'temp',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]);
    if (max(max(max(temp)))>120);temp = temp-273.15;end;
    
    EKEE  = ncread(fname,'u_rms',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2-ncread(fname,'u',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2 + ...
            ncread(fname,'v_rms',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2-ncread(fname,'v',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2; %u-pts
    TdxsqE  = ncread(fname,'temp_dxsq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dxu.^2,[1 1 zL]); %flux-x
    TdysqE  = ncread(fname,'temp_dysq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dyu.^2,[1 1 zL]); %flux-y

% $$$     wvarE = ncread(fname,'wt_rms',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2-ncread(fname,'wt',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).^2;%wt
    TdzsqE  = ncread(fname,'temp_dzsq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]);%wt
    aisoE = ncread(fname,'aiso_bih',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]);
    rey_bihE = ncread(fname,'rey_bih',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]);

    udxsqE = ncread(fname,'u_dxsq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dxt.^2,[1 1 zL]); % flux-y
    udysqE = ncread(fname,'u_dysq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dyt.^2,[1 1 zL]); % flux-x
    vdxsqE = ncread(fname,'v_dxsq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dxt.^2,[1 1 zL]); % flux-y
    vdysqE = ncread(fname,'v_dysq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dyt.^2,[1 1 zL]); % flux-x
    
% $$$     wvarE = cat(3,zeros(xL,yL,1),wvarE);
% $$$     wvarE = (wvarE(:,:,2:end)+wvarE(:,:,1:(end-1)))/2; %t-pts
    
    TdzsqE = cat(3,zeros(xL,yL,1),TdzsqE);
    TdzsqE = (TdzsqE(:,:,2:end)+TdzsqE(:,:,1:(end-1)))/2; %t-pts

    % Do 1D interpolations:
    for xi = 1:(x2-x1)
        sprintf('Doing time %03d of %03d, x %03d of %03d',ti,tL,xi,xL)
        for yi = 1:(y2-y1)
                
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
% $$$                         EKE(ti,Ti,oi) = interp1(tvecU(ziU-1:ziU),squeeze(EKEE(xi,yi,ziU-1:ziU)),Ts(Ti),'linear');
                        EKE(ti,Ti,oi) = nansum(nansum(arear.*((EKEE(xi,yi,ziU)-EKEE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Ts(Ti)-tvecU(ziU-1))+EKEE(xi,yi,ziU-1)),1),2);
                        aiso(ti,Ti,oi) = nansum(nansum(arear.*((aisoE(xi,yi,ziU)-aisoE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Ts(Ti)-tvecU(ziU-1))+aisoE(xi,yi,ziU-1)),1),2);
                        rey_bih(ti,Ti,oi) = nansum(nansum(arear.*((rey_bihE(xi,yi,ziU)-rey_bihE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Ts(Ti)-tvecU(ziU-1))+rey_bihE(xi,yi,ziU-1)),1),2);
                    end
                end

                if (length(ziX)>0)
                    if (ziX>1)
% $$$                         Tdxsq(ti,Ti,oi) = interp1(tvecX(ziX-1:ziX),squeeze(TdxsqE(xi,yi,ziX-1:ziX)),Ts(Ti),'linear');
                        Tdxsq(ti,Ti,oi) = nansum(nansum(arear.*((TdxsqE(xi,yi,ziX)-TdxsqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Ts(Ti)-tvecX(ziX-1))+TdxsqE(xi,yi,ziX-1)),1),2);
                        udysq(ti,Ti,oi) = nansum(nansum(arear.*((udysqE(xi,yi,ziX)-udysqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Ts(Ti)-tvecX(ziX-1))+udysqE(xi,yi,ziX-1)),1),2);
                        vdysq(ti,Ti,oi) = nansum(nansum(arear.*((vdysqE(xi,yi,ziX)-vdysqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Ts(Ti)-tvecX(ziX-1))+vdysqE(xi,yi,ziX-1)),1),2);
                    end
                end

                if (length(ziY)>0)
                    if (ziY>1)
% $$$                         Tdysq(ti,Ti,oi) = interp1(tvecY(ziY-1:ziY),squeeze(TdysqE(xi,yi,ziY-1:ziY)),Ts(Ti),'linear');
                        Tdysq(ti,Ti,oi) = nansum(nansum(arear.*((TdysqE(xi,yi,ziY)-TdysqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Ts(Ti)-tvecY(ziY-1))+TdysqE(xi,yi,ziY-1)),1),2);
                        udxsq(ti,Ti,oi) = nansum(nansum(arear.*((udxsqE(xi,yi,ziY)-udxsqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Ts(Ti)-tvecY(ziY-1))+udxsqE(xi,yi,ziY-1)),1),2);
                        vdxsq(ti,Ti,oi) = nansum(nansum(arear.*((vdxsqE(xi,yi,ziY)-vdxsqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Ts(Ti)-tvecY(ziY-1))+vdxsqE(xi,yi,ziY-1)),1),2);
                    end
                end

                if (length(ziT)>0)
                    if (ziT>1)
% $$$                         wvar(ti,Ti,oi) = interp1(tvecT(ziT-1:ziT),squeeze(wvarE(xi,yi,ziT-1:ziT)),Ts(Ti),'linear');
% $$$                         wvar(ti,Ti,oi) = (wvarE(xi,yi,ziT)-wvarE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Ts(Ti)-tvecT(ziT-1))+wvarE(xi,yi,ziT-1)),1),2);
                        Tdzsq(ti,Ti,oi) = nansum(nansum(arear.*((TdzsqE(xi,yi,ziT)-TdzsqE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Ts(Ti)-tvecT(ziT-1))+TdzsqE(xi,yi,ziT-1)),1),2);
                    end
                end
            end
        end
    end
end

name = [outD model '_' region '_nummixprofiles.mat'];
save(name,'EKE','Tdxsq','Tdysq','Tdzsq','aiso','rey_bih','time','udxsq','vdxsq','udysq','vdysq','NUM','ndays','T','Te','TL');


