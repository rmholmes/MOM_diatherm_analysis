% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
         {'ACCESS-OM2_025deg_jra55_iaf',[17:56]}, ...
       };

rr = 1;
    rr
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model leg legh;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    region = 'Global';
    nyrs = tL/12;
    if (nyrs == round(nyrs))
        szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
        yrs = 1:nyrs;
    else
        nyrs = 1;
        szTe = [TL+1 tL];szT = [TL tL];
    end    
    ycur = 1;
    
    dnum = [];

    %% Load Global Budget:
    for i=1:length(outputs)
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        dnum = [dnum; time];
% $$$     % Annual or Monthly offline Binning:
% $$$     load([base model sprintf('_output%03d_',outputs(i)) 'GlobalHBud_MonAnBin.mat']);
% $$$     GWB = GWBann;

        load([base model sprintf('_output%03d_',outputs(i)) region '_HBud.mat']);
        
        % Fluxes:
        P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
        F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
        Ffz(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.FRZ,szTe); % Surface heat flux (W)
        Fsw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Surface heat flux (W)
        Fsh(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDS,szTe); % Surface heat flux (W)
        Fet(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ETS,szTe); % Surface heat flux (W)
        MNL(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
        M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
        if (isfield(GWB,'VDFkppiw')) % Vertical mixing components
            VDFkppiw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw,szTe);
            VDFkppish(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppish,szTe);
            VDFkppicon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppicon,szTe);
            VDFkppbl(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppbl,szTe);
            VDFkppdd(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppdd,szTe);
            VDFwave(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFwave,szTe);
            VDFnloc(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
            VDFsum(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw+GWB.VDFkppish+GWB.VDFkppicon+ ...
                GWB.VDFkppbl+GWB.VDFkppdd+GWB.VDFwave+GWB.KNL,szTe);
            % Note: May be missing enhanced mixing near rivers
            % (river_diffuse_temp) in ACCESS-OM2
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
            MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);
            M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
        else
            MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'SIG')) % Sigma-diff
            SG(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SIG,szTe);
            M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.SIG,szTe);
        else
            SG(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
            NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
        else
            NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'SUB'))
            SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);
        else
            SUB(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
        end
        D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV,szTe)-GM(:,:,ycur:(ycur+nyrs-1))-SUB(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
        TEN(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe); % Tendency
        ADV(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ADV,szTe); % Advection
        SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
        JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux

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
        JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;

        % Interior heat source P:
        PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));

        % Interior heat source Q:
        QII(:,:,ycur:(ycur+nyrs-1)) = QI(:,:,ycur:(ycur+nyrs-1)) - JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;

        % Across-isotherm advective heat flux:
        CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;

        % External HC Tendency:
        EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;

        % Internal HC Tendency:
        N(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1)) - EHC(:,:,ycur:(ycur+nyrs-1));
% $$$ N(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe);

        % Implicit mixing:
        I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) ...
                                  - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1)) ...
                                  - SUB(:,:,ycur:(ycur+nyrs-1)) - GM(:,:,ycur:(ycur+nyrs-1));

        % Non-advective flux into volume:
        B(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+I(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));

        % Monthly binned Internal HC Tendency:
        if (isfield('GWB','TENMON'))
            Nmon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TENMON,szTe);
        end

        % WMT from B:
        WMTM(:,:,ycur:(ycur+nyrs-1)) = -diff(M(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTF(:,:,ycur:(ycur+nyrs-1)) = -diff(F(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTI(:,:,ycur:(ycur+nyrs-1)) = -diff(I(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMTR(:,:,ycur:(ycur+nyrs-1)) = -diff(R(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        WMT(:,:,ycur:(ycur+nyrs-1)) = WMTM(:,:,ycur:(ycur+nyrs-1))+WMTF(:,:,ycur:(ycur+nyrs-1))+WMTI(:,:,ycur:(ycur+nyrs-1))+WMTR(:,:,ycur:(ycur+nyrs-1));

        % WMT HB from B:
        HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
        HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
        HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
        HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
        HWMT(:,:,ycur:(ycur+nyrs-1)) = HWMTM(:,:,ycur:(ycur+nyrs-1))+HWMTF(:,:,ycur:(ycur+nyrs-1))+HWMTI(:,:,ycur:(ycur+nyrs-1))+HWMTR(:,:,ycur:(ycur+nyrs-1));

        ycur = ycur+nyrs;
    end
    months = [1:length(P(1,:,1))];
    yrs = [1:length(P(1,1,:))];
    
    % Correct time vector zero year:
    dvec = datevec(dnum);
    dvec(:,1) = dvec(:,1) + 1900;
    dnum = datenum(dvec);

    % Load Global V and H:
    Vs = zeros(TL,12,length(yrs));
    Hs = Vs;
    for i=1:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) ...
              '_VHza.mat']);
        Vs(:,:,i) = squeeze(nansum(V,1));
        Hs(:,:,i) = squeeze(nansum(H,1));
    end
    V = cat(1,cumsum(Vs,1,'reverse'),zeros(1,12,length(yrs)));
    H = cat(1,cumsum(Hs,1,'reverse'),zeros(1,12,length(yrs)));
    HE = rho0*Cp*V.*repmat(Te,[1 12 length(yrs)]);
    HI = H - HE;    

%%% Temperature vs. time:
fields = { ...
% $$$           {N, 'Internal HC Tendency $\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
% $$$           {dHdt, 'Total HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$           {EHC, 'External HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$           {F+PI, 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$           {M, 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I, 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {R, 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
% $$$           {M+I+R, 'Total Mixing $\mathcal{M}+\mathcal{I}+\mathcal{R}$',[0 0.5 0],2,'--'}, ...
          {HI, 'Internal Heat Content $\mathcal{H}_I$','m',2,'-'}, ...
          {H, 'Heat Content $\mathcal{H}$','m',2,'-'}, ...
          {HE, 'External Heat Content $\mathcal{H}_E$','m',2,'-'}, ...
          };


% Fluxes:
scale = 1/1e15;label = '(PW)';x = Te;
caxs = [-0.5 0.5];
sp = 0.01;

% Time-integrated fluxes:
Tint = 0;

% Heat Content:
scale = 1/1e23;label = '($10^{23}$J)';x = Te;
caxs = [-1 1];
sp = 0.01;

% Subtract climatology:
Sclim = 1;
climean = [1976 1986];

% Annual average:
AA = 1;

cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

for ii=1:length(fields)
    subplot(1,length(fields),ii);
    tvec = dnum;
    if (length(fields{ii}{1}(:,1)) == length(Te))
        x = Te;
    else
        x = T;
    end
    % Annual average:
    if (AA)
        V = squeeze(monmean(fields{ii}{1},2,ndays))'*scale;
        tvec = unique(dvec(:,1));
        % Subtract climateology:
        if (Sclim)
            V = V-repmat(mean(V(find(tvec>=climean(1) & tvec<=climean(2)),:),1),[length(tvec) 1]);
        end
        % Time-integrate:
        if (Tint == 1)
            V = cumsum(V*86400*365,1);
        end
        [X,Y] = ndgrid(tvec,x);
    end
    contourf(X,Y,V,cint,'linestyle','none');
    cb = colorbar('Location','NorthOutside','FontSize',15);    
    set(gca,'ytick',-5:5:35);
    ylabel(cb,label);
    if (~AA)
        datetick(gca,'x');
    end
    ylim([-3 31]);
    grid on;
    caxis(caxs);
    if (ii == 1)
        ylabel('Temperature ($^\circ$C)');
    end
    xlabel('Date');
    title(fields{ii}{2});
    set(gca,'FontSize',15);
end
colormap(redblue);


%%% ARGO data:

load('ArgoData.mat');
yrsargo = 2004:2014;
Teargo = Tedge;
Hargo = zeros(length(Teargo),length(yrsargo));
Vargo = Hargo;
dvec = datevec(days);
for yi=1:length(yrsargo)
    yr = yrsargo(yi);
    Hargo(:,yi) = mean(Hsnap(:,dvec(:,1)==yr),2);
    Vargo(:,yi) = mean(Vsnap(:,dvec(:,1)==yr),2);
end
HIargo = Hargo-rho0*Cp*Vargo.*repmat(Teargo',[1 length(yrsargo)]);
HEargo = Hargo-HIargo;

[X,Y] = ndgrid(yrsargo,Teargo);

Hargo = Hargo-repmat(mean(Hargo,2),[1 length(yrsargo)]);
HEargo = HEargo-repmat(mean(HEargo,2),[1 length(yrsargo)]);
HIargo = HIargo-repmat(mean(HIargo,2),[1 length(yrsargo)]);

% Fluxes:
scale = 1/1e23;label = '(10$^{23}$J)';
caxs = [-0.5 0.5];
sp = 0.025;

cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

fields = { ...
          {[], 'Internal Heat Content $\mathcal{H}_I$','m',2,'-'}, ...
          {[], 'Heat Content $\mathcal{H}$','m',2,'-'}, ...
          {[], 'External Heat Content $\mathcal{H}_E$','m',2,'-'}, ...
          };

figure;
subplot(1,3,1);
contourf(X,Y,HIargo'*scale,cint,'linestyle','none');
ylabel('Temperature ($^\circ$C)');
cb = colorbar('Location','NorthOutside','FontSize',15);    
set(gca,'ytick',-5:5:35);
ylabel(cb,label);
xlabel('Date');
ylim([-3 31]);
grid on;
caxis(caxs);
xlim([2004 2014]);
title(fields{1}{2});
set(gca,'FontSize',15);
subplot(1,3,2);
contourf(X,Y,Hargo'*scale,cint,'linestyle','none');
cb = colorbar('Location','NorthOutside','FontSize',15);    
set(gca,'ytick',-5:5:35);
ylabel(cb,label);
xlabel('Date');
ylim([-3 31]);
grid on;
caxis(caxs);
xlim([2004 2014]);
title(fields{2}{2});
set(gca,'FontSize',15);
subplot(1,3,3);
contourf(X,Y,HEargo'*scale,cint,'linestyle','none');
cb = colorbar('Location','NorthOutside','FontSize',15);    
set(gca,'ytick',-5:5:35);
ylabel(cb,label);
xlabel('Date');
ylim([-3 31]);
grid on;
caxis(caxs);
xlim([2004 2014]);
title(fields{3}{2});
set(gca,'FontSize',15);
colormap(redblue);
