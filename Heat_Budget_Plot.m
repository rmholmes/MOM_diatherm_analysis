% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
% ACCESS-OM2 1-degree:
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may',[0]}, ... % KDS50 base
         {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may',[36]}, ... % KDS50 base
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_kb1em5',[0]}, ... %KDS50 with kb = 1e-5 background diffusivity
       };

% $$$ rr = 1;
for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    ndays = diff(time_snap);
    ndays = ndays(1:12);
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
% $$$ region = 'Pacific';
    
    ycur = 1;
    %% Global Calculations:
    for i=1:length(outputs)
        
% $$$     % Annual or Monthly offline Binning:
% $$$     load([base model sprintf('_output%03d_',outputs(i)) 'GlobalHBud_MonAnBin.mat']);
% $$$     GWB = GWBann;

        load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
        nyrs = tL/12;szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
        
        % Fluxes:
        P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
        F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
        M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
        if (isfield(GWB,'VDFkppiw')) % Vertical mixing components
            VDFkppiw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw,szTe);
            VDFkppish(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppish,szTe);
            VDFkppicon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppicon,szTe);
            VDFkppbl(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppbl,szTe);
            VDFkppdd(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppdd,szTe);
            VDFwave(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFwave,szTe);
            KPPnloc(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
            VDFsum(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw+GWB.VDFkppish+GWB.VDFkppicon+ ...
                GWB.VDFkppbl+GWB.VDFkppdd+GWB.VDFwave+GWB.KNL,szTe);
        end
        if (isfield(GWB,'RED')) % Redi Diffusion
            R(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.RED+GWB.K33,szTe); % Redi diffusion (W)
        else
            R(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end
        if (isfield(GWB,'NGM')) % GM parameterization
            GM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NGM,szTe); % GM (W)
        else
            GM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'MDS')) % Mix-downslope
            MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);; % GM (W)
            M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
        else
            MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
            NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
        else
            NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV-GWB.SUB,szTe)-GM(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
        SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
        JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux
        SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);

        % Pacific Interior fluxes:
        if (strcmp(region,'Pacific'))
            JI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.JBS+GWB.JSP+GWB.JITF,szTe); %Combined volume flux out
            QI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.QBS+GWB.QSP+GWB.QITF,szTe); %Combined heat flux out
        else
            QI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
            JI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end

        % Snapshot fields:
        dVdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dVdt,szTe); % V Change (m3s-1)
        dHdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dHdt,szTe); % H Change (W)

        % Water-mass transformation:
        G(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)) - JS(:,:,ycur:(ycur+nyrs-1)) + JI(:,:,ycur:(ycur+nyrs-1)); %Water-mass transformation (m3s-1)

        % Surface Volume flux base flux (not P!)
        JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % Interior heat source P:
        PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));

        % Interior heat source Q:
        QII(:,:,ycur:(ycur+nyrs-1)) = QI(:,:,ycur:(ycur+nyrs-1)) - JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % Across-isotherm advective heat flux:
        CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % External HC Tendency:
        EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;

        % Internal HC Tendency:
        N(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1)) - EHC(:,:,ycur:(ycur+nyrs-1));
% $$$ N(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe);

% $$$ % Alternative method 1 for the N calculation:
% $$$ N(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*cumsum(dVdt*dT,1,'reverse');
% However, this does not work well, potentially because the
% integral should conceptually then be defined on Tcenters as
% opposed to Tedges. In any case, this method gives a non-zero
% total heat flux due to implicit mixing and is much noisier. 

        % Implicit mixing:
        I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1));

        % Non-advective flux into volume:
        B(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+I(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));

        % Monthly binned Internal HC Tendency:
        if (isfield('GWB','TENMON'))
            Nmon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TENMON,szTe);
        end

        % Alternative method 2 (which is what I was using for
        % the GRL submit - but the above method is simpler to explain). 
% $$$ Ialt(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1))-D(:,:,ycur:(ycur+nyrs-1))-CIA(:,:,ycur:(ycur+nyrs-1)) + QI(:,:,ycur:(ycur+nyrs-1));
% $$$ Balt(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+Ialt(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));
% $$$ Nalt(:,:,ycur:(ycur+nyrs-1)) = Balt(:,:,ycur:(ycur+nyrs-1)) + PI(:,:,ycur:(ycur+nyrs-1)) - QII(:,:,ycur:(ycur+nyrs-1));
% Gives very close results. The difference between N and Nalt, and
% I and Ialt (i.e. (Ialt-I)./I, (Nalt-N)./N) is less than 1e-5 in
% all cases (except where I==0, which gives Inf).

        % WMT from B:
        WMTM(:,:,ycur:(ycur+nyrs-1)) = -diff(M(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTF(:,:,ycur:(ycur+nyrs-1)) = -diff(F(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTI(:,:,ycur:(ycur+nyrs-1)) = -diff(I(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTR(:,:,ycur:(ycur+nyrs-1)) = -diff(R(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMT(:,:,ycur:(ycur+nyrs-1)) = WMTM(:,:,ycur:(ycur+nyrs-1))+WMTF(:,:,ycur:(ycur+nyrs-1))+WMTI(:,:,ycur:(ycur+nyrs-1))+WMTR(:,:,ycur:(ycur+nyrs-1));

        % WMT HB from B:
        HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
        HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
        HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
        HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
        HWMT(:,:,ycur:(ycur+nyrs-1)) = HWMTM(:,:,ycur:(ycur+nyrs-1))+HWMTF(:,:,ycur:(ycur+nyrs-1))+HWMTI(:,:,ycur:(ycur+nyrs-1))+HWMTR(:,:,ycur:(ycur+nyrs-1));

        % Alternative method 3 I from Volume budget (as for spatial structure calc FlI):
% $$$ WMTI(:,:,ycur:(ycur+nyrs-1)) = avg(dVdt(:,:,ycur:(ycur+nyrs-1)),1) - avg(JS(:,:,ycur:(ycur+nyrs-1)),1)-WMTM(:,:,ycur:(ycur+nyrs-1))-WMTF(:,:,ycur:(ycur+nyrs-1))-WMTR(:,:,ycur:(ycur+nyrs-1));
% $$$ I(:,:,ycur:(ycur+nyrs-1)) = zeros(size(I(:,:,ycur:(ycur+nyrs-1))));
% $$$ I(1:(end-1),:,ycur:(ycur+nyrs-1)) = rho0*Cp*cumsum(WMTI(:,:,ycur:(ycur+nyrs-1))*dT,1,'reverse');
% This also gives pretty consistent results, although has similar
% noisy problems and a non-zero total heat flux at the cooler
% temperatures as alternative method 1 above. This method is
% slightly better at the warmest temperatures.
        ycur = ycur+nyrs;
    end
    months = [1:length(P(1,:,1))];
    yrs = [1:length(P(1,1,:))];
%%%%Heat Flux:
% Production fields:
fields = { ...
          {N(:,months,yrs), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
          {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
          {M(:,months,yrs)+R(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{R}+\mathcal{I}$','r',2,'-'}, ...
          {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
          {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
          {NUM(:,months,yrs), 'Numerical Mixing Direct $\mathcal{I}$','c',2,'-'}, ...
% $$$           {SW(:,months,yrs), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$           {dHdt(:,months,yrs), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',2,'--'}, ...
% $$$           {PI(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}_I$',[0.49 0.18 0.56],2,'--'}, ...
% $$$           {F(:,months,yrs), 'Surface Heat Fluxes $\mathcal{F}$','k',2,'-'}, ...
% $$$           {P(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}$',[0.49 0.18 0.56],2,'-'}, ...
% $$$           {GM(:,months,yrs), 'GM $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$           {SUB(:,months,yrs), 'SUB $\mathcal{S}$',[0 0.5 0],2,':'}, ...
% $$$           {HWMTI(:,months,yrs), 'Advective Implicit Mixing','b',2,'--'}, ...
% $$$           {HWMTM(:,months,yrs), 'Advective Vertical Mixing','r',2,'--'}, ...
% $$$           {HWMTF(:,months,yrs), 'Advective Surface Forcing','k',2,'--'}, ...
% $$$           {M(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           {Nmon(:,months,yrs), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$           {CIA(:,months,yrs), 'Across-Isotherm Advection $\mathcal{G}\Theta\rho_0C_p$',[0.49 0.18 0.56],2,'--'}, ...
% $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
          };
% $$$ fields = { ...
% $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {VDFkppiw(:,months,yrs), 'Vertical Diffusion KPPIW','b',2,'-'}, ...
% $$$           {VDFkppish(:,months,yrs), 'Vertical Diffusion KPPISH',[0 0.5 0],2,'-'}, ...
% $$$           {VDFkppbl(:,months,yrs), 'Vertical Diffusion KPPBL','m',2,'-'}, ...
% $$$           {VDFwave(:,months,yrs), 'Vertical Diffusion WAVE','k',2,'-'}, ...
% $$$ % $$$           {VDFkppicon(:,months,yrs), 'Vertical Diffusion KPPICON','y',2,'-'}, ...
% $$$ % $$$           {KPPnloc(:,months,yrs), 'KPP Non-local',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {VDFkppdd(:,months,yrs), 'Vertical Diffusion KPPDD','m',2,'-'}, ...
% $$$           {VDFkppdd(:,months,yrs)+VDFkppicon(:,months,yrs)+KPPnloc(:, ...
% $$$                                                   months,yrs),'DD + KPP non-local + Int. Convection',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$           {VDFsum(:,months,yrs), 'Vertical Mixing SUM','r',2,'--'}, ...
% $$$           };

Fscale = 1/1e15;

yrtyps = {'-','--','-.',':','-'}; % line-types for different years

%Fluxes only:
figure;
set(gcf,'Position',[207          97        1609         815]);
leg = {};
legh = [];
for i=1:length(fields)
    hold on;
    if (length(fields{i}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    for j=1:length(yrs) % Plot years separately
        h = plot(Te,monmean(fields{i}{1}(:,:,yrs(j)),2,ndays(months))*Fscale,yrtyps{j}, 'color',fields{i}{3} ...
             ,'linewidth',3);
        if (j == 1)
            legh(i) = h;
        end
    end
% $$$     % Average years together:
% $$$     legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
    leg{i} = fields{i}{2};
% $$$     leg{i} = strrep(RUNS{rr}{1},'_',' ');
end
ylim([-1.5 1.5]);
xlim([-3 31]);
box on; 
grid on;
ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
xlabel('Temperature $\Theta$ ($^\circ$C)');
lg = legend(legh,leg);
set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

%%% Temperature vs. time:
months = [1:12];
fields = { ...
% $$$           {N(:,months,yrs), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
          {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$           {M(:,months,yrs)+R(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{R}+\mathcal{I}$','r',2,'-'}, ...
          };

% Fluxes:
scale = 1/1e15;label = '(PW)';x = Te;
caxs = [-1.5 0];
sp = 0.05;
caxs = [-5 5];
sp = 0.5;

% $$$ % Transformations:
% $$$ scale = 1/1e6;label = '(Sv)';
% $$$ caxs = [-250 250];
% $$$ sp = 25;
% $$$ caxs = [-70 70];
% $$$ sp = 3.5;

cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

figure;
for ii=1:length(fields)
    subplot(1,length(fields),ii);
    % V = mean(fields{ii}{1},3)'*scale;
    V = reshape(fields{ii}{1}*scale,[TL+1 12*length(yrs)])';
    if (length(fields{ii}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
% $$$     [X,Y] = ndgrid(months,x);
    [X,Y] = ndgrid(1:(12*length(yrs)),x);
    contourf(X,Y,V,cint);%,'linestyle','none');
    cb = colorbar('Location','NorthOutside','FontSize',25);    
    set(gca,'ytick',-5:5:35);
    set(gca,'xtick',[1:(12*length(yrs))]);
    ylim([-3 31]);
    grid on;
    caxis(caxs);
    xlabel('Month');
    ylabel('Temperature ($^\circ$C)');
    xlabel(cb,strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
    xlabel(cb,[model ' ' fields{ii}{2} ' ' ...
               label],'FontSize',20);
    set(gca,'FontSize',15);
end
cmap = redblue((length(cint)-3)*2);
cmap = cmap(1:(length(cint)-3),:);
colormap(cmap);
colormap(redblue);


%%% Spatial Structure:

VAR = 'FlM';
TYPE = 'VertInt';
Tl = 22.5;
yr = 1;
name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
eval(['load(name,''' VAR ''');']);
eval([VAR '(isnan(' VAR ')) = 0.0;']);
eval([VAR 'a = ' VAR '(:,:,((yr-1)*12+1):yr*12);']);

% $$$     for i=2:length(outputs)
% $$$         name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$         eval(['load(name,''' VAR ''');']);
% $$$         eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$         eval([VAR 'a = ' VAR 'a + ' VAR ';']);
% $$$     end
% $$$     eval([VAR ' = ' VAR 'a/length(outputs);']);
% $$$     eval([VAR '(' VAR '==0) = NaN;']);
% $$$     eval(['FlM = ' VAR ';']);
%%% Plot spatial pattern:

    try
        obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
        LAND = obj.SST(:,:,1);
    catch
        LAND = zeros(size(FlM(:,:,1)));
    end

    %If MOM01, fix NaN's in grid:
    if (strfind(model,'01'))
        lon = repmat(lon(:,500),[1 yL]);
        latv = nanmean(lat,1);
        lat = repmat(latv,[xL 1]);
    end

    [xL,yL] = size(lon);
    xvec = 1:1:xL;
    yvec = 1:1:yL;
    txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

    months = {[1:12], ...
              [3], ...
              [7], ...
              [11]};
    labels = {'(a) Annual', ...
              '(b) March', ...
              '(c) July', ...
              '(d) November'};
% $$$     months = {[1:12]}; labels = {'Annual'};
    %Colormap:
    clim = [-150 0];
    sp = 1;

    cpts = [-1e10 clim(1):sp:clim(2) 1e10];
    npts = length(cpts)
    cmap = flipud(lbmap(npts-3,'RedBlue'));
    cmap = redblue(npts-3);
    doWMT = 0;
    if (doWMT)
    for i=1:(npts-3)
        if (cmap(i,:) == 1.0)
            cmap(i,:) = [0.94 0.94 0.94];
        end
    end
    else
        cmap = parula(npts-3);
        cmap = parula(npts-3);
        cmap(end,:) = [0.97 0.97 0.8];
        cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
    end

    tmp = LAND;
    tmp(isnan(LAND)) = clim(1)-sp/2;
    tmp(~isnan(LAND)) = NaN;
    LAND = tmp;
    cmap(2:(end+1),:) = cmap;
    cmap(1,:) = [0 0 0];

    climn = [clim(1)-sp clim(2)];
    
figure;
set(gcf,'Position',[3          59        1916         914]);
set(gcf,'defaulttextfontsize',20);
set(gcf,'defaultaxesfontsize',20);

poss = [0.1300    0.4553    0.7693    0.4697; ...
        0.1300    0.1389    0.2343    0.2680; ...
        0.3951    0.1389    0.2343    0.2680; ...
        0.6681    0.1389    0.2343    0.2680];
for i=1:length(months)
        if (i == 1)
            subplot(5,3,[1 9]);
        else
            subplot(5,3,[10 13]+(i-2));
        end
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    if (length(months{i})>1)
        tmp = FlM;
        tmp(isnan(tmp)) = 0.0;
        Z = monmean(tmp(:,:,months{i}),3,ndays(months{i}));
        Z(Z == 0) = NaN;
    else
        Z = FlM(:,:,months{i});
    end
    Z = Z(xvec,yvec);
    if (doWMT)
        Z = Z*86400;
    end
    
    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z,cpts,'linestyle','none');
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    if (i==1)
        cb = colorbar;
        if (~doWMT)
            ylabel(cb,'Wm$^{-2}$');
        else
            ylabel(cb,'m/day');
        end            
        ylim(cb,clim);
    end
    hold on;
    if (i>1)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i<=2)
        ylabel('Latitude ($^\circ$N)');
    end
    
    set(gca,'Position',[poss(i,:)]);
    ylim([-45 45]);
    set(gca,'ytick',[-45:15:45]);
    set(gca,'FontSize',20);
    colormap(cmap);
end

%%% Plot spatial pattern of net heat flux and SST:

% Load Variable and calculate mean:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
yr = 1;
shflux = shflux(:,:,((yr-1)*12+1):yr*12);
SST = SST(:,:,((yr-1)*12+1):yr*12);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux = shfluxa/length(outputs);
% $$$ SST = SSTa/length(outputs);
if (max(max(SST))>120);SST = SST-273.15;end;

%If MOM01, fix NaN's in grid:
if (strfind(model,'01'))
    lon = repmat(lon(:,500),[1 yL]);
    latv = nanmean(lat,1);
    lat = repmat(latv,[xL 1]);
end

%Sum of all positives:
shfluxA = mean(shflux,3);
shfluxP = nansum(shfluxA(shfluxA>0).*area(shfluxA>0))/1e15
shfluxM = nansum(shfluxA(shfluxA<0).*area(shfluxA<0))/1e15

%Sum of mean flux above mean isotherm:
SSTA = mean(SST,3);
isot = 21.5;
shfluxP = nansum(shfluxA(SSTA>=isot).*area(SSTA>=isot))/1e15
shfluxM = nansum(shfluxA(SSTA<isot).*area(SSTA<isot))/1e15

[xL,yL] = size(lon);
xvec = 1:1:xL;
yvec = 1:1:yL;
txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

months = {[1:12], ...
                            [3], ...
                            [7], ...
                            [11]};

labels = {'(a) Annual', ...
          '(b) March', ...
          '(c) July', ...
          '(d) November'};

clim = [-200 200];
sp = 20;

cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)
cmap = redblue(npts-3);

figure;
set(gcf,'Position',[3          59        1916         914]);
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',20);
poss = [0.1300    0.415  0.7693    0.56; ...
        0.1300    0.1    0.2343    0.2680; ...
        0.3951    0.1    0.2343    0.2680; ...
        0.6681    0.1    0.2343    0.2680];
for i=1:length(months)
    if (i == 1)
        subplot(5,3,[1 9]);
    else
        subplot(5,3,[10 13]+(i-2));
    end
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
    Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
    Z = Z(xvec,yvec);
    Z2 = Z2(xvec,yvec);
    Z(Z==0) = NaN;
    
    % Map projection:
% $$$     map = 1;
% $$$     ax = worldmap('World');
% $$$     setm(ax, 'Origin', [0 260 0])
% $$$     contourfm(Y,X,Z,[-1e10 -500:20:500 1e10],'linestyle','none');
    contourf(X,Y,Z,cpts,'linestyle','none');
    hold on;
% $$$     [c,h] = contourm(Y,X,Z2,[-3:2:21 25:2:35],'-k');
    [c,h] = contour(X,Y,Z2,[-3:2:35],'-k');
    clabel(c,h);        
% $$$     hand = clabelm(c,h);        
% $$$     set(hand,'FontSize',10,'BackgroundColor','none');
    if (i == 1)
        [c,h] = contour(X,Y,Z2,[21.5 21.5],'-k','linewidth',2);
        clabel(c,h)
% $$$         [c,h] = contourm(Y,X,Z2,[23 23],'-k','linewidth',2);
% $$$         hand = clabelm(c,h);
% $$$         set(hand,'FontSize',10,'BackgroundColor','none');
    end
    caxis(clim);
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
    end
    ylim([-75 75]);
    if (i>1)
        xlabel('Longitude ($^\circ$E)');
    end
    if (i<=2)
        ylabel('Latitude ($^\circ$N)');
    end
    if (i>1)
% $$$         text(-277,53,labels{i},'BackgroundColor','w');
        text(-277,64,labels{i},'BackgroundColor','w','margin',0.01);
        set(gca,'xtick',[-240:60:60]);
    else
% $$$         text(-279,55,labels{i},'BackgroundColor','w');
        text(-279,69,labels{i},'BackgroundColor','w','margin',0.01);
        set(gca,'xtick',[-270:30:60]);
    end        
    set(gca,'ytick',[-60:30:60]);
    set(gca,'Position',[poss(i,:)]);
    set(gca,'color','k');

% $$$     setm(gca, 'MlabelParallel', 'south');
% $$$     setm(gca, 'MlabelParallel', 'south');   
% $$$     land = shaperead('landareas.shp', 'UseGeoCoords', true);
% $$$     geoshow(land, 'FaceColor', [0 0 0])
end 
colormap(cmap);

end
% $$$ 
% $$$ %%%%%%%%%%%% Supp FIGURES
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ %% Plot SST difference and shflux difference plots:
% $$$ 
% $$$ % $$$ % Load MOM01 data:
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/';
% $$$ % $$$ model = 'MOM01';
% $$$ % $$$ outputs = [111 222];
% $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ ndays = diff(time_snap);
% $$$ % $$$ 
% $$$ % $$$ [xL,yL] = size(lon)
% $$$ % $$$ %If MOM01, fix NaN's in grid:
% $$$ % $$$ if (strfind(model,'01'))
% $$$ % $$$     lonv = lon(:,500);
% $$$ % $$$     latv = nanmean(lat,1);
% $$$ % $$$     [tmp ind] = min(abs(latv - 60))
% $$$ % $$$     [lon,lat] = ndgrid(lonv,latv(1:ind));
% $$$ % $$$ end
% $$$ % $$$ MOM01lon = lon;
% $$$ % $$$ MOM01lat = lat;
% $$$ % $$$ 
% $$$ % $$$ load([base model sprintf('_output%03d_SurfaceVars.mat', ...
% $$$ % $$$                          outputs(1))]);
% $$$ % $$$ MOM01shflux = shflux(:,1:ind,:);
% $$$ % $$$ MOM01SST = SST(:,1:ind,:);
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = 94;
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ outputs = [14]
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em6/mat_data/';
% $$$ model = 'MOM025_kb1em6';
% $$$ outputs = 30;
% $$$ 
% $$$ % Load ACCESS-OM2:
% $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ outputs = 3;
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ SSTa = SSTa(:,:,13:24);
% $$$ ndays = ndays(1:12);
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux = shfluxa/length(outputs);
% $$$ SST = SSTa/length(outputs);
% $$$ 
% $$$ if (max(max(max(SST)))>100)
% $$$     SST = SST-273.15;
% $$$ end
% $$$ 
% $$$ % WOA13 SST:
% $$$ WOAname = '/srv/ccrc/data03/z3500785/WOA13/woa13_decav_t00_04v2.nc';
% $$$ WOASST = ncread(WOAname,'t_an',[1 1 1 1],[1440 720 1 1]);
% $$$ [WOAlon,WOAlat] = ndgrid(ncread(WOAname,'lon'),ncread(WOAname,'lat'));
% $$$ 
% $$$ %Shift longitudes:
% $$$ [tmp ind] = min(abs(WOAlon(:,1)-80));
% $$$ WOASST = cat(1,WOASST(ind+1:end,:),WOASST(1:ind,:));
% $$$ WOAlon = cat(1,WOAlon(ind+1:end,:)-360,WOAlon(1:ind,:));
% $$$ WOAlat = cat(1,WOAlat(ind+1:end,:),WOAlat(1:ind,:));
% $$$ 
% $$$ % Calculate bias from WOA:
% $$$ SSTbias = monmean(SST,3,ndays)-interp2(WOAlon',WOAlat',WOASST',lon,lat,'linear');
% $$$ % $$$ for i=1:12
% $$$ % $$$     SST(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01SST(:,:,i)',lon,lat,'linear')-SST(:,:,i);
% $$$ % $$$     shflux(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01shflux(:,:,i)',lon,lat,'linear')-shflux(:,:,i);
% $$$ % $$$ end
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:2:xL;
% $$$ yvec = 1:2:yL;
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ poss = [0.15  0.56 0.3 0.34; ...
% $$$         0.48  0.56 0.3 0.34; ...
% $$$         0.15  0.1 0.3 0.34; ...
% $$$         0.48  0.1 0.3 0.34];
% $$$ 
% $$$ contourf(lon(xvec,yvec),lat(xvec,yvec),SSTbias(xvec,yvec),[-1e10 -5:0.25:5 1e10],'linestyle', ...
% $$$          'none');
% $$$ hold on;
% $$$ % $$$ contour(lon,lat,SSTbias,[-3 -2 -1 1 2 3],'-k');
% $$$ contour(lon,lat,SSTbias,[-1:0.5:1],'-k');
% $$$ set(gca,'color','k');
% $$$ title('$\kappa_B=10^{-5}$ - WOA13 SST Year 2 ($^\circ$C)');
% $$$ caxis([-1 1]);
% $$$ colorbar;
% $$$ colormap(redblue(24));
% $$$ set(gca,'FontSize',20);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ % $$$ ylabel('Latitude ($^\circ$N)');
% $$$ set(gca,'Position',poss(4,:));
% $$$ 
% $$$ 
% $$$ %
% $$$ for i=1:length(months)
% $$$     if (i == 1)
% $$$         subplot(5,3,[1 9]);
% $$$     else
% $$$         subplot(5,3,[10 13]+(i-2));
% $$$     end
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
% $$$     Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$     Z = Z(xvec,yvec);
% $$$     Z2 = Z2(xvec,yvec);
% $$$ % $$$     contourf(X,Y,Z2,[-1e10 -5:0.25:5 1e10],'linestyle','none');
% $$$     contourf(X,Y,Z,[-1e10 -500:10:500 1e10],'linestyle','none');
% $$$     hold on;
% $$$ % $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$ % $$$     if (i==1)
% $$$ % $$$         [c,h] = contour(X,Y,Z,[-100:10:-10],'--k');
% $$$ % $$$         [c,h] = contour(X,Y,Z,[10:10:100],'-k');
% $$$ % $$$     end
% $$$ % $$$     else
% $$$ % $$$         [c,h] = contour(X,Y,Z2,[-3:4:35],'-k');
% $$$ % $$$     end
% $$$ % $$$     clabel(c,h);
% $$$     caxis([-100 100]);
% $$$     if (i==1)
% $$$         cb = colorbar;
% $$$ % $$$         ylabel(cb,'$^\circ$C');
% $$$         ylabel(cb,'Wm$^{-2}$');
% $$$     end
% $$$     ylim([-75 60]);
% $$$     if (i>1)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (i<=2)
% $$$         ylabel('Latitude ($^\circ$N)');
% $$$     end
% $$$     if (i>1)
% $$$         text(-276,53,labels{i},'BackgroundColor','w');
% $$$     else
% $$$         text(-278,55,labels{i},'BackgroundColor','w');
% $$$     end        
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'color','k');
% $$$ end 
% $$$ colormap(redblue);
% $$$ 
