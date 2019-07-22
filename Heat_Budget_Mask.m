function [mask_t,mask_u,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask(region,gname,wname,outD,model)
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

% Mask already exists in saved file? 
outname = [outD model '_RegionMask_' region '.mat'];
if (exist(outname))
    ['Loading mask that already exists from file ' outname]
    load(outname);
else    

lonv_u = ncread(gname,'xu_ocean');
latv_u = ncread(gname,'yu_ocean');
if (exist(wname))
    T = ncread(wname,'neutral');
    TL = length(T);
end
lon = ncread(gname,'geolon_t');
lat = ncread(gname,'geolat_t');
mask_t_full = ~isnan(ncread(gname,'kmt'));
mask_u_full = ~isnan(ncread(gname,'kmu')); %kmu mask is the same as
                                           %kmt mask
[xL,yL] = size(lon);

mask_t = ones(xL,yL); %Mask: 1 = water in region, 0 = outside
                      %region/not water
mask_u = ones(xL,yL); % mask on u-points
mask_Ny = 0*mask_t; %North y-trans mask
mask_Nx = 0*mask_t; %North x-trans mask
mask_Sy = 0*mask_t; %South y-trans mask
mask_Sx = 0*mask_t; %South x-trans mask
mask_Wy = 0*mask_t; %West y-trans mask
mask_Wx = 0*mask_t; %West x-trans mask

made_mask = 0;

%%%%%%%% PACIFIC REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(region,'Pacific'))
    ['Generating new mask to file ' outname]
    made_mask = 1;

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

mask_t(mask_t_full == 0) = 0;

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

%%%%%%%% INDO-PACIFIC REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(region,'IndoPacific'))
    made_mask = 1;
    % Get Pacific mask:
    [mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('Pacific',gname,wname,outD,model);

    % Add eastern Indian:
    [~, t1] = min(abs(latv_u+46));[~, t3] = min(abs(lonv_u+73));
    mask_Sy(1:t3,t1) = 1;
    mask_Wx = 0*mask_Wx;
    mask_Wy = 0*mask_Wy;
    
    [~, t2] = min(abs(latv_u-40));[~, t3] = min(abs(lonv_u+200));
    mask_t(1:t3,(t1+1):t2) = 1;

    % Add western Indian:
    [~, t3] = min(abs(lonv_u-20));[~, t2] = min(abs(latv_u-30));
    mask_t(t3:end,(t1+1):t2) = 1;
    [~, t3] = min(abs(lonv_u-47));[~, t2] = min(abs(latv_u-32));
    mask_t(t3:end,(t1+1):t2) = 1;

    mask_t(mask_t_full == 0) = 0;
%%%%%%%% ATLANTIC REGION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(region,'Atlantic'))
    made_mask = 1;
    % Get Indo-Pacific mask:
    [mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('IndoPacific',gname,wname,outD,model);
    
    % Invert:
    tmp = mask_t;
    tmp(mask_t == 1) = 0;
    tmp(mask_t == 0) = 1;
    mask_t = tmp;
    [~, t1] = min(abs(latv_u+46));[~, t2] = min(abs(latv_u-66));
    
    mask_t(:,1:t1) = 0;
    mask_t(:,(t2+1):end) = 0;

    mask_t(mask_t_full == 0) = 0;
    

%%%%%%% Zonal-average two-basin masks: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(region,'IndoPacific2BAS'))
    ['Generating new mask to file ' outname]
    made_mask = 1;

    % Northern boundary at 66N:
    [~, t1] = min(abs(latv_u-66));
    mask_t(:,(t1+1):end) = 0;
    mask_u(:,(t1+1):end) = 0;
    
    % Southern boundary at 34S:
    [~, t1] = min(abs(latv_u+34));
    mask_t(:,1:t1) = 0;
    mask_u(:,1:(t1-1)) = 0;
    
% $$$     % Fix Bering straight bit:
% $$$     [~, t1] = min(abs(latv_u-66.5));
% $$$     [~, t2] = min(abs(latv_u-65));
% $$$     [~, t3] = min(abs(lonv_u+180));
% $$$     [~, t4] = min(abs(lonv_u+178));
% $$$     mask_t(t3:t4,t2:t1) = 1;
% $$$     mask_u((t3-1):t4,(t2-1):t1) = 1;

    % Zonal definitions:
    [~, t1] = min(abs(lonv_u+68));mask_t(t1:end,:) = 0;mask_u(t1:end,:) = 0;
    [~, t1] = min(abs(lonv_u+98));[~, t2] = min(abs(latv_u-18));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    [~, t1] = min(abs(lonv_u+89));[~, t2] = min(abs(latv_u-15));
    mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    [~, t1] = min(abs(lonv_u+84));
    [~, t2] = min(abs(latv_u-10));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    [~, t1] = min(abs(lonv_u+82.5));[~, t2] = min(abs(latv_u-9.5));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    [~, t2] = min(abs(lonv_u+80.5));[~, t3] = min(abs(latv_u-9));
    mask_t(t1:t2,t3:end) = 0;mask_u(t1:(t2-1),t3:end) = 0;
    [~, t1] = min(abs(lonv_u+78.25));
    [~, t2] = min(abs(latv_u-9));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    [~, t1] = min(abs(lonv_u+77.5));[~, t2] = min(abs(latv_u-7.25));mask_t(t1:end,t2:end) = 0;mask_u(t1:end,t2:end) = 0;
    mask_t(799,535)=0;
    mask_u(798:800,534:536)=0;
    
    % Add western Indian:
    [~, t1] = min(abs(latv_u+34));
    [~, t3] = min(abs(lonv_u-20));[~, t2] = min(abs(latv_u-30));
    mask_t(t3:end,(t1+1):t2) = 1;
    mask_u((t3-1):end,t1:t2) = 1;
    [~, t3] = min(abs(lonv_u-47));[~, t2] = min(abs(latv_u-32));
    mask_t(t3:end,(t1+1):t2) = 1;
    mask_u((t3-1):end,t1:t2) = 1;

    mask_t(mask_t_full == 0) = 0;
    mask_u(mask_u_full == 0) = 0;
    
elseif (strcmp(region,'Atlantic2BAS'))
    ['Generating new mask to file ' outname]
    made_mask = 1;

    % Get IndoPacific2BAS:
    [mask_t,mask_u,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('IndoPacific2BAS',gname,wname,outD,model);
    
    % Invert:
    tmp = mask_t;
    tmp(mask_t == 1) = 0;
    tmp(mask_t == 0) = 1;
    mask_t = tmp;
    tmp = mask_u;
    tmp(mask_u == 1) = 0;
    tmp(mask_u == 0) = 1;
    mask_u = tmp;
    
    % Only North of 34S:
    [~, t1] = min(abs(latv_u+34));
    mask_t(:,1:t1) = 0;
    mask_u(:,1:(t1-1)) = 0;

    % For u-mask add Bering Strait u-points:
    [~, t1] = min(abs(latv_u-66));
% $$$     [~, t2] = min(abs(lonv_u+178));
% $$$     mask_u(t2:end,t1) = 1;
    mask_u(:,t1) = 1;
    
    mask_t(mask_t_full == 0) = 0;
    mask_u(mask_u_full == 0) = 0;

elseif (strcmp(region,'AtlanticNZ'))
    ['Generating new mask to file ' outname]
    made_mask = 1;
    % Get Atlantic mask:
    [mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('Atlantic',gname,wname,outD,model);
    
    [~, t1] = min(abs(latv_u+34));
    
    mask_t(:,1:t1) = 0;
elseif (strcmp(region,'IndoPacificNZ'))
    ['Generating new mask to file ' outname]
    made_mask = 1;
    [mask_t,mask_Ny,mask_Nx,mask_Sx,mask_Sy,mask_Wx,mask_Wy] = ...
    Heat_Budget_Mask('IndoPacific',gname,wname,outD,model);
    
    [~, t1] = min(abs(latv_u+34));
    
    mask_t(:,1:t1) = 0;
elseif (strcmp(region,'SO_IndoPacific'))
    ['Generating new mask to file ' outname]
    made_mask = 1;

    % Northern boundary at 34S:
    [~, t1] = min(abs(latv_u+34));
    mask_t(:,(t1+1):end) = 0;
    mask_u(:,(t1+1):end) = 0;

    % Zonal definitions:
    [~, t1] = min(abs(lonv_u-23));
    [~, t2] = min(abs(lonv_u+70.5));
    mask_t((t2+1):(t1-1),:) = 0;mask_u(t2:(t1-1),:) = 0;

    mask_t(mask_t_full == 0) = 0;
    mask_u(mask_u_full == 0) = 0;
elseif (strcmp(region,'SO_Atlantic'))
    ['Generating new mask to file ' outname]
    made_mask = 1;

    % Northern boundary at 34S:
    [~, t1] = min(abs(latv_u+34));
    mask_t(:,(t1+1):end) = 0;
    mask_u(:,(t1+1):end) = 0;

    % Zonal definitions:
    [~, t1] = min(abs(lonv_u-23));mask_t(t1:end,:) = 0;mask_u(t1:end,:) = 0;
    [~, t1] = min(abs(lonv_u+70.5));mask_t(1:t1,:) = 0;mask_u(1:(t1-1),:) = 0;

    mask_t(mask_t_full == 0) = 0;
    mask_u(mask_u_full == 0) = 0;
end

if (made_mask)
    save(outname,'mask_t','mask_u','mask_Ny','mask_Nx','mask_Sx','mask_Sy', ...
         'mask_Wx','mask_Wy','-v7.3');
end
end

end

