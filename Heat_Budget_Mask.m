function [mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask(region,gname,fname,wname)
% This function provides masks for different regions for analysing
% the heat budget. 

% Mask information:  -----------------------------------------
%
% Notes on the B-grid:
% - All the heat budget diagnostics are on t-cells, so the mask
%   throughout the Pacific mask uses T-cells.
% - tx_trans_nrho is on the east-face of the T-cells, or the
%   south-face of the u-cells.
% - tx_trans_nrho is on the north-frace of the T-cells, or the
%   west-face of the u-cells. 
%
% The masks below were checked using these grid locations and the
% vertical sum of the transport.

lonv_u = ncread(gname,'xu_ocean');
latv_u = ncread(gname,'yu_ocean');
T = ncread(wname,'neutral');
TL = length(T);
lon = ncread(gname,'geolon_t');
lat = ncread(gname,'geolat_t');
[xL,yL] = size(lon);

mask_t = ones(xL,yL); %Pacific Mask: 1 = Pacific water, 0 = Elsewhere
mask_Ny = 0*mask_t; %North y-trans mask
mask_Nx = 0*mask_t; %North x-trans mask
mask_Sy = 0*mask_t; %South y-trans mask
mask_Sx = 0*mask_t; %South x-trans mask
mask_Wy = 0*mask_t; %West y-trans mask
mask_Wx = 0*mask_t; %West x-trans mask

%%%%%%%% PACIFIC REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(region,'Pacific'))

% ITF segments:
% 114.9E:
[~, t1] = min(abs(lonv_u+245.1));[~, t2] = min(abs(latv_u+23));[~, t3] = min(abs(latv_u+8.25));
mask_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;
% 256.9E:
[~, t1] = min(abs(lonv_u+256.9));[~, t2] = min(abs(latv_u+0.875));[~, t3] = min(abs(latv_u-4.121));
mask_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;
% 254.25E:
[~, t1] = min(abs(lonv_u+254.25));[~, t2] = min(abs(latv_u+6.362));[~, t3] = min(abs(latv_u+4.371));
mask_Wx(t1,t2:t3) = 1;mask_t(1:t1,t2:t3) = 0;

%Checked: YES! Pacific_MASK_CheckITF.png where magenta points are
%txtrans(mask_Wx).

% SOUTH Pacific:
% 46S:
[~, t1] = min(abs(latv_u+46));[~, t2] = min(abs(lonv_u+213.5));[~, t3] = min(abs(lonv_u+73));
mask_Sy((t2+1):t3,t1) = 1;mask_t(:,1:t1) = 0;

% 213.5E:
[~, t1] = min(abs(lonv_u+213.5));[~, t2] = min(abs(latv_u+46));[~, t3] = min(abs(latv_u+37));
mask_Sx(t1,(t2+1):t3) = 1;mask_t(1:t1,t2:t3) = 0;

%Checked: YES! Pacific_MASK_CheckTassieNZ.png where magenta crosses
%are txtrans(mask_Sx) and magenta diamonds are tytrans(mask_Sy).

% NORTH Pacific:
[~, t1] = min(abs(latv_u-66));[~, t2] = min(abs(lonv_u+172));[~, t3] = min(abs(lonv_u+166));
mask_Ny(t2:t3,t1) = 1;
mask_t(:,(t1+1):end) = 0;

%Checked: YES! Pacific_MASK_CheckBering.png 

% Pacific Mask:
[~, t1] = min(abs(lonv_u+68));mask_t(t1:end,:) = 0;
[~, t1] = min(abs(lonv_u+98));[~, t2] = min(abs(latv_u-18));mask_t(t1:end,t2:end) = 0;
[~, t1] = min(abs(lonv_u+89));[~, t2] = min(abs(latv_u-15));
mask_t(t1:end,t2:end) = 0;[~, t1] = min(abs(lonv_u+84));
[~, t2] = min(abs(latv_u-10));mask_t(t1:end,t2:end) = 0;
[~, t1] = min(abs(lonv_u+82.5));[~, t2] = min(abs(latv_u-9.5));mask_t(t1:end,t2:end) = 0;
[~, t2] = min(abs(lonv_u+80.5));[~, t3] = min(abs(latv_u-9));
mask_t(t1:t2,t3:end) = 0;[~, t1] = min(abs(lonv_u+78.25));
[~, t2] = min(abs(latv_u-9));mask_t(t1:end,t2:end) = 0;
[~, t1] = min(abs(lonv_u+77.5));[~, t2] = min(abs(latv_u-7.25));mask_t(t1:end,t2:end) = 0;
[~, t1] = min(abs(lonv_u+217));[~, t2] = min(abs(latv_u+22.29));mask_t(1:t1,1:t2) = 0;
[~, t1] = min(abs(lonv_u+261));mask_t(1:t1,:) = 0;
[~, t1] = min(abs(lonv_u+256.75));[~, t2] = min(abs(latv_u+1));mask_t(1:t1,1:t2) = 0;
[~, t1] = min(abs(lonv_u+252));[~, t2] = min(abs(latv_u+6.5));mask_t(1:t1,1:t2) = 0;
[~, t1] = min(abs(lonv_u+248));[~, t2] = min(abs(latv_u+7.25));mask_t(1:t1,1:t2) = 0;
[~, t1] = min(abs(lonv_u+260));[~, t2] = min(abs(latv_u-8));mask_t(1:t1,1:t2) = 0;
[~, t1] = min(abs(lonv_u+258.5));[~, t2] = min(abs(latv_u-6.5));mask_t(1:t1,1:t2) = 0;
mask_t(81,525) = 0;mask_t(799,535)=0;

SST = ncread(fname,'temp',[1 1 1 1],[xL yL 1 1]);
mask_t(isnan(SST)) = 0;

% $$$ txtrans = sum(ncread(wname,'tx_trans_nrho',[1 1 1 1],[xL yL TL 1]),3);
% $$$ tytrans = sum(ncread(wname,'ty_trans_nrho',[1 1 1 1],[xL yL TL 1]),3);
% $$$ 
% $$$ % ITF check:
% $$$ % $$$ [tmp ln1] = min(abs(lonv_u+265));
% $$$ % $$$ [tmp ln2] = min(abs(lonv_u+242));
% $$$ % $$$ [tmp lt1] = min(abs(latv_u+25));
% $$$ % $$$ [tmp lt2] = min(abs(latv_u-15));
% $$$ [tmp ln1] = min(abs(lonv_u+220));
% $$$ [tmp ln2] = min(abs(lonv_u+190));
% $$$ % $$$ [tmp ln1] = min(abs(lonv_u+190));
% $$$ % $$$ [tmp ln2] = min(abs(lonv_u+70));
% $$$ [tmp lt1] = min(abs(latv_u+50));
% $$$ [tmp lt2] = min(abs(latv_u+36));
% $$$ % $$$ [tmp ln1] = min(abs(lonv_u+180));
% $$$ % $$$ [tmp ln2] = min(abs(lonv_u+160));
% $$$ % $$$ [tmp lt1] = min(abs(latv_u-60));
% $$$ % $$$ [tmp lt2] = min(abs(latv_u-70));
% $$$ SST = SST(ln1:ln2,lt1:lt2);
% $$$ lon = lon(ln1:ln2,lt1:lt2);
% $$$ lat = lat(ln1:ln2,lt1:lt2);
% $$$ mask_t = mask_t(ln1:ln2,lt1:lt2);
% $$$ txtrans = txtrans(ln1:ln2,lt1:lt2);
% $$$ tytrans = tytrans(ln1:ln2,lt1:lt2);
% $$$ 
% $$$ % $$$ mask_Wx = mask_Wx(ln1:ln2,lt1:lt2);
% $$$ % $$$ txtrans_Wx = txtrans; txtrans_Wx(mask_Wx~=1) = 0;
% $$$ % $$$ mask_Wy = mask_Wy(ln1:ln2,lt1:lt2);
% $$$ % $$$ tytrans_Wy = tytrans; tytrans_Wy(mask_Wy~=1) = 0;
% $$$ mask_Sx = mask_Sx(ln1:ln2,lt1:lt2);
% $$$ txtrans_Sx = txtrans; txtrans_Sx(mask_Sx~=1) = 0;
% $$$ mask_Sy = mask_Sy(ln1:ln2,lt1:lt2);
% $$$ tytrans_Sy = tytrans; tytrans_Sy(mask_Sy~=1) = 0;
% $$$ % $$$ mask_Nx = mask_Nx(ln1:ln2,lt1:lt2);
% $$$ % $$$ txtrans_Nx = txtrans; txtrans_Nx(mask_Nx~=1) = 0;
% $$$ % $$$ mask_Ny = mask_Ny(ln1:ln2,lt1:lt2);
% $$$ % $$$ tytrans_Ny = tytrans; tytrans_Ny(mask_Ny~=1) = 0;
% $$$ 
% $$$ clf;
% $$$ plot(lon(isnan(SST)),lat(isnan(SST)),'or');
% $$$ hold on;
% $$$ plot(lon(~isnan(SST)),lat(~isnan(SST)),'ok');
% $$$ plot(lon(mask_t==1),lat(mask_t==1),'oc');
% $$$ plot(lon(txtrans~=0)+0.125,lat(txtrans~=0),'xg');
% $$$ hold on;
% $$$ plot(lon(tytrans~=0),lat(tytrans~=0)+0.05,'dy');
% $$$ plot(lon(txtrans_Sx~=0)+0.125,lat(txtrans_Sx~=0),'xm');
% $$$ plot(lon(tytrans_Sy~=0),lat(tytrans_Sy~=0)+0.05,'dm');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

