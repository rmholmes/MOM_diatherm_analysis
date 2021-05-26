% This script makes plots of the spatial structure of the
% diathermal fluxes in the MOM simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[81]}, ...
       };
cols = {'b','r','k','m','g'};

rr = 1;
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
        ndays = ndays(1:12);
    end
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
    nyrs = tL/12;
    if (round(nyrs)~=nyrs)
        anavg = 1;
        nyrs = tL;
        months = {[1:1]};
    else
        anavg = 0;
        months = {[1:12]};
    end
    yrs = 1:nyrs;
    ycur = 1;


%% Variances:
Tl = 15;

% FAKE IT, with ndays:
ndays = [31 28 31 30 31 30 31 31 30 31 30 31];

% $$$ VARS = {'udxsq','vdxsq','udysq','vdysq','rey_bih','aiso','EKE','wvar','Tdxsq','Tdysq','Tdzsq'};
VARS = {'udxsq','vdxsq','udysq','vdysq','Tdxsq','Tdysq','Tdzsq'};
TYPE = 'variances';
name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
for ii=1:length(VARS)
    VAR = VARS{ii}
    eval(['load(name,''' VAR ''');']);
    eval([VAR '(isnan(' VAR ')) = 0.0;']);
    if (length(outputs)==1)
        eval([VAR ' = reshape(' VAR ',[length(' VAR '(:,1,1)) length(' VAR '(1,:,1)) 12 nyrs]);']);
    else
        eval([VAR 'a = ' VAR ';']);
        for i=2:length(outputs)
            name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
            eval(['load(name,''' VAR ''');']);
            eval([VAR '(isnan(' VAR ')) = 0.0;']);
            eval([VAR 'a = ' VAR 'a + ' VAR ';']);
        end
        eval([VAR ' = ' VAR 'a/length(outputs);']);
    end
    eval([VAR ' = mean(monmean(' VAR ',3,ndays),4);']);
    eval([VAR '(' VAR '==0) = NaN;']);
end
Tdh = sqrt(0.5*(Tdxsq+Tdysq));
Tdv = sqrt(Tdzsq);
udh = sqrt(0.25*(udxsq+vdxsq+udysq+vdysq));
udhd = sqrt(0.5*(udxsq+vdysq));
udhc = sqrt(0.5*(udysq+vdxsq));

VARS = {'udhd','udhc','Tdh','Tdv'};
names = {'(c) $\sqrt{\frac{1}{2}\left(|\Delta_x u|^2+|\Delta_y v|^2\right)}$', ...
         '(d) $\sqrt{\frac{1}{2}\left(|\Delta_y u|^2+|\Delta_x v|^2\right)}$', ...
         '(c) $\sqrt{\frac{1}{2}\left(|\Delta_x \Theta|^2+|\Delta_y \Theta|^2\right)}$', ...
         '(e) $\sqrt{|\Delta_z \Theta|^2}$'};
units = {'$ms^{-1}$','$ms^{-1}$','$^\circ C$','$^\circ C$'};
clims = {[0 0.06],[0 0.06],[0 1],[0 6]};
    
%%% Plot spatial pattern:

    try
        obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
        LAND = obj.SST(:,:,1);
    catch
        LAND = zeros(size(FlM(:,:,1)));
    end

    [xL,yL] = size(lon);
    xvec = 1:4:xL;
    yvec = 1:4:yL;

    %Colormaps:
% $$$     clims = {[0 0.02],[0 0.3e-7],[0 4e-9],[0 20]};
    nlv = 50;

    cmapbase = parula(nlv-3);
    cmapbase(end,:) = [0.97 0.97 0.8];
    cmapbase(end-1,:) = (cmapbase(end-1,:)+cmapbase(end,:))/2;
    cmapbase = flipud(cmapbase);
    for i=1:length(VARS)
        sp = (clims{i}(2)-clims{i}(1))/(nlv-3);
        cpts{i} = [-1e10 clims{i}(1):sp:clims{i}(2) 1e10];
        cmap{i} = cmapbase;
        LANDmask{i} = LAND;
        LANDmask{i}(isnan(LAND)) = clims{i}(1)-sp/2;
        LANDmask{i}(~isnan(LAND)) = NaN;
        cmap{i}(2:(end+1),:) = cmapbase;
        cmap{i}(1,:) = [0 0 0];
        climns{i} = [clims{i}(1)-sp clims{i}(2)];
    end
    
%Mean of all months:
figure;
set(gcf,'Position',[1921           1        1920        1005]);
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

    poss = [0.1300    0.58      0.3548    0.3692; ...
            0.5700    0.58      0.3548    0.3692; ...
            0.1300    0.1100    0.3548    0.3692; ...
            0.5700    0.1100    0.3548    0.3692];
            
for i=1:length(VARS)
    subplot(2,2,i);
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    eval(['Z = ' VARS{i} '(xvec,yvec);']);
    Z(Z<clims{i}(1)) = clims{i}(1);
    contourf(X,Y,Z,cpts{i},'linestyle','none');
    hold on;    
    contourf(X,Y,LANDmask{i}(xvec,yvec),climns{i},'linestyle','none');
    caxis(climns{i});
    cb = colorbar;
    ylabel(cb,units{i});
    ylim(cb,clims{i});
    if (i>=3)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i==1 | i == 3)
        ylabel('Latitude ($^\circ$N)');
    end
    set(gca,'xtick',[-270:30:60]);
    set(gca,'ytick',[-75:15:75]);
    ylim([-60 60]);
    colormap(gca,cmap{i});
    title('$15^\circ$C isotherm');
    text(-275,54,[names{i}],'BackgroundColor','w','FontSize',13);% ' on ' num2str(Tl) '$^\circ$C isotherm']);
    set(gca,'Position',poss(i,:));
end

