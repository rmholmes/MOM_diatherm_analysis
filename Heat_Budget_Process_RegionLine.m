% This script processes a particular region for num-mix and
% variance variance terms as a function of temperature.

region = [-215 -180 30 45];
regname = 'Kuroshio';

baseL = '/g/data/e14/rmh561/access-om2/archive/';

model = 'ACCESS-OM2_025deg_jra55_ryf_norediGM';
baseD = [baseL '025deg_jra55_ryf_norediGM/']; %Data Directory.

outD = [baseD 'mat_data/'];
 
post = 'ocean/';

%%% Do num-mix curves:
outputs = [76:80];

output = outputs(1);
load([outD model sprintf('_output%03d',output) '_BaseVars.mat']);

[tmp x1] = min(abs(xt-region(1)));
[tmp x2] = min(abs(xt-region(2)));
[tmp y1] = min(abs(yt-region(3)));
[tmp y2] = min(abs(yt-region(4)));

NUM = zeros(TL+1,tL,length(outputs));
for oi = 1:length(outputs)
    output = outputs(oi);
    base = [baseD sprintf('output%03d/',output) post];
    wname = [base 'ocean_wmass.nc'];

    for ti=1:tL
        sprintf('Num-mix out %03d of %03d, time %03d of %03d',oi,length(outputs),ti,tL)
        for Ti = (TL+1):-1:1
            NUM(Ti,ti,oi) = nansum(nansum(area(x1:x2,y1:y2).*ncread(wname, ...
                         'temp_numdiff_heat_on_nrho',[x1 y1 Ti ti],[x2-x1+1 y2-y1+1 1 1]),1),2);
        end
    end
end

%%% Do variances:
baseL = '/scratch/e14/rmh561/access-om2/archive/';
baseD = [baseL '025deg_jra55_ryf_norediGM/']; %Data Directory.

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
% $$$ udysq  = NaN*zeros(tL,TL,length(outputs));
% $$$ vdxsq  = NaN*zeros(tL,TL,length(outputs));
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
% $$$     udysqE = ncread(fname,'u_dysq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dyt.^2,[1 1 zL]); % flux-x
% $$$     vdxsqE = ncread(fname,'v_dxsq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dxt.^2,[1 1 zL]); % flux-y
    vdysqE = ncread(fname,'v_dysq',[x1 y1 1 ti],[x2-x1+1 y2-y1+1 zL 1]).*repmat(dyt.^2,[1 1 zL]); % flux-x
    
% $$$     wvarE = cat(3,zeros(xL,yL,1),wvarE);
% $$$     wvarE = (wvarE(:,:,2:end)+wvarE(:,:,1:(end-1)))/2; %t-pts
    
    TdzsqE = cat(3,zeros(x2-x1+1,y2-y1+1,1),TdzsqE);
    TdzsqE = (TdzsqE(:,:,2:end)+TdzsqE(:,:,1:(end-1)))/2; %t-pts

    % Do 1D interpolations:
    EKEs = NaN*zeros(x2-x1,y2-y1,TL);
    aisos = NaN*zeros(x2-x1,y2-y1,TL);
    rey_bihs = NaN*zeros(x2-x1,y2-y1,TL);
    Tdxsqs = NaN*zeros(x2-x1,y2-y1,TL);
% $$$     udysqs = NaN*zeros(x2-x1,y2-y1,TL);
    vdysqs = NaN*zeros(x2-x1,y2-y1,TL);
    Tdysqs = NaN*zeros(x2-x1,y2-y1,TL);
    udxsqs = NaN*zeros(x2-x1,y2-y1,TL);
% $$$     vdxsqs = NaN*zeros(x2-x1,y2-y1,TL);
    Tdzsqs = NaN*zeros(x2-x1,y2-y1,TL);
    
    for xi = 1:(x2-x1)
        sprintf('Doing time %03d of %03d, x %03d of %03d',ti,tL,xi,x2-x1+1)
        for yi = 1:(y2-y1)
                
            tvecU = squeeze(temp(xi,yi,:)+temp(xi+1,yi,:)+temp(xi,yi+1,:)+temp(xi+1,yi+1,:))/4;
                tvecX = squeeze(temp(xi,yi,:)+temp(xi+1,yi,:))/2;
                tvecY = squeeze(temp(xi,yi,:)+temp(xi,yi+1,:))/2;
            tvecT = squeeze(temp(xi,yi,:));
            for Ti = 1:TL+1

                ziU = find(tvecU<Te(Ti),1,'first');
                ziX = find(tvecX<Te(Ti),1,'first');
                ziY = find(tvecY<Te(Ti),1,'first');
                ziT = find(tvecT<Te(Ti),1,'first');

                if (length(ziU)>0)
                    if (ziU>1)
% $$$                         EKEs(xi,yi,Ti) = interp1(tvecU(ziU-1:ziU),squeeze(EKEE(xi,yi,ziU-1:ziU)),Te(Ti),'linear');
                        EKEs(xi,yi,Ti) = (EKEE(xi,yi,ziU)-EKEE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Te(Ti)-tvecU(ziU-1))+EKEE(xi,yi,ziU-1);
                        aisos(xi,yi,Ti) = (aisoE(xi,yi,ziU)-aisoE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Te(Ti)-tvecU(ziU-1))+aisoE(xi,yi,ziU-1);
                        rey_bihs(xi,yi,Ti) = (rey_bihE(xi,yi,ziU)-rey_bihE(xi,yi,ziU-1))/(tvecU(ziU)-tvecU(ziU-1))*(Te(Ti)-tvecU(ziU-1))+rey_bihE(xi,yi,ziU-1);
                    end
                end

                if (length(ziX)>0)
                    if (ziX>1)
% $$$                         Tdxsqs(xi,yi,Ti) = interp1(tvecX(ziX-1:ziX),squeeze(TdxsqE(xi,yi,ziX-1:ziX)),Te(Ti),'linear');
                        Tdxsqs(xi,yi,Ti) = (TdxsqE(xi,yi,ziX)-TdxsqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Te(Ti)-tvecX(ziX-1))+TdxsqE(xi,yi,ziX-1);
% $$$                         udysqs(xi,yi,Ti) = (udysqE(xi,yi,ziX)-udysqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Te(Ti)-tvecX(ziX-1))+udysqE(xi,yi,ziX-1);
                        vdysqs(xi,yi,Ti) = (vdysqE(xi,yi,ziX)-vdysqE(xi,yi,ziX-1))/(tvecX(ziX)-tvecX(ziX-1))*(Te(Ti)-tvecX(ziX-1))+vdysqE(xi,yi,ziX-1);
                    end
                end

                if (length(ziY)>0)
                    if (ziY>1)
% $$$                         Tdysqs(xi,yi,Ti) = interp1(tvecY(ziY-1:ziY),squeeze(TdysqE(xi,yi,ziY-1:ziY)),Te(Ti),'linear');
                        Tdysqs(xi,yi,Ti) = (TdysqE(xi,yi,ziY)-TdysqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Te(Ti)-tvecY(ziY-1))+TdysqE(xi,yi,ziY-1);
                        udxsqs(xi,yi,Ti) = (udxsqE(xi,yi,ziY)-udxsqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Te(Ti)-tvecY(ziY-1))+udxsqE(xi,yi,ziY-1);
% $$$                         vdxsqs(xi,yi,Ti) = (vdxsqE(xi,yi,ziY)-vdxsqE(xi,yi,ziY-1))/(tvecY(ziY)-tvecY(ziY-1))*(Te(Ti)-tvecY(ziY-1))+vdxsqE(xi,yi,ziY-1);
                    end
                end

                if (length(ziT)>0)
                    if (ziT>1)
% $$$                         wvars(xi,yi,Ti) = interp1(tvecT(ziT-1:ziT),squeeze(wvarE(xi,yi,ziT-1:ziT)),Te(Ti),'linear');
% $$$                         wvars(xi,yi,Ti) = (wvarE(xi,yi,ziT)-wvarE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Te(Ti)-tvecT(ziT-1))+wvarE(xi,yi,ziT-1);
                        Tdzsqs(xi,yi,Ti) = (TdzsqE(xi,yi,ziT)-TdzsqE(xi,yi,ziT-1))/(tvecT(ziT)-tvecT(ziT-1))*(Te(Ti)-tvecT(ziT-1))+TdzsqE(xi,yi,ziT-1);
                    end
                end
            end
        end
    end

    EKE(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*EKEs,1),2);
    aiso(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*aisos,1),2);
    rey_bih(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*rey_bihs,1),2);
    Tdxsq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*Tdxsqs,1),2);
% $$$     udysq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*udysqs,1),2);
    vdysq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*vdysqs,1),2);
    Tdysq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*Tdysqs,1),2);
    udxsq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*udxsqs,1),2);
% $$$     vdxsq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*vdxsqs,1),2);
    Tdzsq(ti,:,oi) = nansum(nansum(repmat(arear(1:end-1,1:end-1),[1 1 TL]).*Tdzsqs,1),2);
end
end

name = [outD model '_' region '_nummixprofiles.mat'];
save(name,'EKE','Tdxsq','Tdysq','Tdzsq','aiso','rey_bih','time','udxsq','vdxsq','udysq','vdysq','NUM','ndays','T','Te','TL');


