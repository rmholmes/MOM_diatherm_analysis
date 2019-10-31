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
if (strcmp(region,'IndoPacific2BAS'))
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
    mask_u(:,t1) = 1;
    
    mask_t(mask_t_full == 0) = 0;
    mask_u(mask_u_full == 0) = 0;

end

if (made_mask)
    save(outname,'mask_t','mask_u','mask_Ny','mask_Nx','mask_Sx','mask_Sy', ...
         'mask_Wx','mask_Wy','-v7.3');
end
end

end

