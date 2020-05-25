% This script makes temperature-time global plots of the heat budget
% in MOM5.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
         {'ACCESS-OM2_1deg_jra55_rdf',[51:55],[1972]}, ...
         {'ACCESS-OM2_1deg_jra55_rdf_pert',[51:55],[1972]}, ...
       };

rr = 1;
    rr
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

% $$$     clearvars -except base RUNS rr outputs model leg legh;
    
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
    DT = [];

    %% Load Global Budget:
    for i=1:length(outputs)
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        dnum = [dnum; time];
        DT = [DT; ndays];
% $$$     % Annual or Monthly offline Binning:
% $$$     load([base model sprintf('_output%03d_',outputs(i)) 'GlobalHBud_MonAnBin.mat']);
% $$$     GWB = GWBann;

% $$$         load([base model sprintf('_output%03d_',outputs(i)) region '_HBud.mat']);
% $$$         
% $$$         % Fluxes:
% $$$         P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
% $$$         F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
% $$$         Ffz(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.FRZ,szTe); % Surface heat flux (W)
% $$$         Fsw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Surface heat flux (W)
% $$$         Fsh(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDS,szTe); % Surface heat flux (W)
% $$$         Fet(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ETS,szTe); % Surface heat flux (W)
% $$$         MNL(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
% $$$         M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
% $$$         if (isfield(GWB,'VDFkppiw')) % Vertical mixing components
% $$$             VDFkppiw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw,szTe);
% $$$             VDFkppish(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppish,szTe);
% $$$             VDFkppicon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppicon,szTe);
% $$$             VDFkppbl(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppbl,szTe);
% $$$             VDFkppdd(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppdd,szTe);
% $$$             VDFwave(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFwave,szTe);
% $$$             VDFnloc(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
% $$$             VDFsum(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw+GWB.VDFkppish+GWB.VDFkppicon+ ...
% $$$                 GWB.VDFkppbl+GWB.VDFkppdd+GWB.VDFwave+GWB.KNL,szTe);
% $$$             % Note: May be missing enhanced mixing near rivers
% $$$             % (river_diffuse_temp) in ACCESS-OM2
% $$$         end
% $$$         if (isfield(GWB,'RED')) % Redi Diffusion
% $$$             R(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.RED+GWB.K33,szTe); % Redi diffusion (W)
% $$$         else
% $$$             R(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end
% $$$         if (isfield(GWB,'NGM')) % GM parameterization
% $$$             GM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NGM,szTe); % GM (W)
% $$$         else
% $$$             GM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'MDS')) % Mix-downslope
% $$$             MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);
% $$$             M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
% $$$         else
% $$$             MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'SIG')) % Sigma-diff
% $$$             SG(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SIG,szTe);
% $$$             M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.SIG,szTe);
% $$$         else
% $$$             SG(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
% $$$             NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
% $$$         else
% $$$             NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'SUB'))
% $$$             SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);
% $$$         else
% $$$             SUB(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end
% $$$         D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV,szTe)-GM(:,:,ycur:(ycur+nyrs-1))-SUB(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
% $$$         TEN(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe); % Tendency
% $$$         ADV(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ADV,szTe); % Advection
% $$$         SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
% $$$         JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux
% $$$ 
% $$$         % Pacific Interior fluxes:
% $$$         if (strcmp(region,'Pacific'))
% $$$             JI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.JBS+GWB.JSP+GWB.JITF,szTe); %Combined volume flux out
% $$$             QI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.QBS+GWB.QSP+GWB.QITF,szTe); %Combined heat flux out
% $$$         else
% $$$             QI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$             JI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end
% $$$ 
% $$$         % Snapshot fields:
% $$$         dVdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dVdt,szTe); % V Change (m3s-1)
% $$$         dHdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dHdt,szTe); % H Change (W)
% $$$ 
% $$$         % Water-mass transformation:
% $$$         G(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)) - JS(:,:,ycur:(ycur+nyrs-1)) + JI(:,:,ycur:(ycur+nyrs-1)); %Water-mass transformation (m3s-1)
% $$$ 
% $$$         % Surface Volume flux base flux (not P!)
% $$$         JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;
% $$$ 
% $$$         % Interior heat source P:
% $$$         PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Interior heat source Q:
% $$$         QII(:,:,ycur:(ycur+nyrs-1)) = QI(:,:,ycur:(ycur+nyrs-1)) - JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;
% $$$ 
% $$$         % Across-isotherm advective heat flux:
% $$$         CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;
% $$$ 
% $$$         % External HC Tendency:
% $$$         EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 tL nyrs])*rho0*Cp;
% $$$ 
% $$$         % Internal HC Tendency:
% $$$         N(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1)) - EHC(:,:,ycur:(ycur+nyrs-1));
% $$$ % $$$ N(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe);
% $$$ 
% $$$         % Implicit mixing:
% $$$         I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) ...
% $$$                                   - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1)) ...
% $$$                                   - SUB(:,:,ycur:(ycur+nyrs-1)) - GM(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Non-advective flux into volume:
% $$$         B(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+I(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Monthly binned Internal HC Tendency:
% $$$         if (isfield('GWB','TENMON'))
% $$$             Nmon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TENMON,szTe);
% $$$         end
% $$$ 
% $$$         % WMT from B:
% $$$         WMTM(:,:,ycur:(ycur+nyrs-1)) = -diff(M(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTF(:,:,ycur:(ycur+nyrs-1)) = -diff(F(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTI(:,:,ycur:(ycur+nyrs-1)) = -diff(I(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTR(:,:,ycur:(ycur+nyrs-1)) = -diff(R(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMT(:,:,ycur:(ycur+nyrs-1)) = WMTM(:,:,ycur:(ycur+nyrs-1))+WMTF(:,:,ycur:(ycur+nyrs-1))+WMTI(:,:,ycur:(ycur+nyrs-1))+WMTR(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % WMT HB from B:
% $$$         HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
% $$$         HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
% $$$         HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
% $$$         HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 tL nyrs]);
% $$$         HWMT(:,:,ycur:(ycur+nyrs-1)) = HWMTM(:,:,ycur:(ycur+nyrs-1))+HWMTF(:,:,ycur:(ycur+nyrs-1))+HWMTI(:,:,ycur:(ycur+nyrs-1))+HWMTR(:,:,ycur:(ycur+nyrs-1));

        ycur = ycur+nyrs;
    end
% $$$     months = [1:length(P(1,:,1))];
% $$$     yrs = [1:length(P(1,1,:))];
    
    % Correct time vector zero year:
    dvec = datevec(dnum);
    if (length(RUNS{rr}{3})== 1)
        dvec(:,1) = dvec(:,1)-dvec(1,1) + RUNS{rr}{3}(1);
    else
        dvec(:,1) = dvec(:,1) + 1900;
    end
    dnum = datenum(dvec);
    yrs = 1:10;

    % Load Global V and H:
    load([base model sprintf('_output%03d',outputs(i)) '_VHza.mat']);
    tLV = length(V(1,1,:));
    Vs = zeros(TL,tLV,length(outputs));
    Hs = Vs;
    for i=1:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) ...
              '_VHza.mat']);
        Vs(:,:,i) = squeeze(nansum(V,1));
        Hs(:,:,i) = squeeze(nansum(H,1));
    end
    if (tLV ~= 12)
        nyrs = tLV*length(outputs);
        V = reshape(Vs,[TL nyrs]);
        H = reshape(Hs,[TL nyrs]);
        yrs = [1:nyrs];
    else
        'ERROR: V & H MONTHLY!'
    end
    V = cat(1,cumsum(V,1,'reverse'),zeros(1,nyrs));
    H = cat(1,cumsum(H,1,'reverse'),zeros(1,nyrs));
    HE = rho0*Cp*V.*repmat(Te,[1 nyrs]);
    HI = H - HE;    
    
    % Remap to ocean percentile:
    Vtot = V(1,:);
    p_ofT = V./repmat(Vtot,[TL+1 1]);

    p = linspace(0,1,1000);
    pl = length(p);
    T_ofp = zeros(pl,length(yrs));
    H_ofp = T_ofp;
    HI_ofp = T_ofp;
    HE_ofp = T_ofp;
    for yi = 1:length(yrs)
        T_ofp(:,yi) = interp1(p_ofT(:,yi)+(1:TL+1)'/1e10,Te,p,'linear');
        T_ofp(1,yi) = Te(end);
        H_ofp(:,yi) = interp1(Te,H(:,yi),T_ofp(:,yi),'linear');
        HE_ofp(:,yi) = interp1(Te,HE(:,yi),T_ofp(:,yi),'linear');
        HI_ofp(:,yi) = interp1(Te,HI(:,yi),T_ofp(:,yi),'linear');
    end

    % Remap to ocean volume:
    Vtot = max(V(1,:));
% $$$     p_ofT = V./repmat(Vtot,[TL+1 1]);

    Vi = linspace(0,Vtot,1000);
    Vl = length(Vi);
    T_ofV = zeros(Vl,nyrs);%length(yrs));
    H_ofV = T_ofV;
    HI_ofV = T_ofV;
    HE_ofV = T_ofV;
    for yi = 1:nyrs
        T_ofV(:,yi) = interp1(flipud(V(:,yi))+(1:(TL+1))'/(TL+1)*Vtot/1e12,flipud(Te),Vi,'linear');
        T_ofV(1,yi) = Te(end);
        H_ofV(:,yi) = interp1(Te,H(:,yi),T_ofV(:,yi),'linear');
        HE_ofV(:,yi) = interp1(Te,HE(:,yi),T_ofV(:,yi),'linear');
        HI_ofV(:,yi) = interp1(Te,HI(:,yi),T_ofV(:,yi),'linear');
    end

% $$$     % Remap ocean percentiles to depth:
% $$$     Vz = zeros(zL,1);%length(V(1,1,:)));
% $$$     for i=1:length(outputs)
% $$$         load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
% $$$         dnum = [dnum; time];
% $$$ 
% $$$         file = load([base model sprintf('_output%03d_',outputs(i)) 'VHofz.mat']);
% $$$         Vz = Vz+monmean(file.V,2,ndays);
% $$$     end
% $$$     Vz = Vz/length(outputs);
% $$$     p_Vz = cumsum(Vz)/sum(Vz);
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
% $$$           {HI_ofp, '$\overline{\Theta}(p,t)-\Theta(p,t) = \mathcal{H}_I/(\rho_0 C_p V_T p)$','m',2,'-'}, ...
% $$$           {H_ofp, '$\overline{\Theta}(p,t) = \mathcal{H}/(\rho_0 C_p V_T p)$','m',2,'-'}, ...
% $$$           {HE_ofp, '$\Theta(p,t) = \mathcal{H}_E/(\rho_0 C_p V_T p)$','m',2,'-'}, ...
          {T_ofV, '$\Theta(\mathcal{V},t)$','m',2,'-'}, ...
          };


% Scales and labels:

% Fluxes:
scale = 1/1e15;label = '(PW)';x = Te;
caxs = [-0.5 0.5];
sp = 0.01;

% Temperature:
scale = 1;label = '$^\circ$C';x = Te;
caxs = [-0.15 0.15];
sp = 0.005;

% Time-integrate fluxes:
Tint = 0;

% remapping for percentiles:
premap = 0;

% Subtract climatology:
Sclim = 1;
climean = [1972 1981];

% Annual average:
AA = 0;

% $$$ % remap to depth:
% $$$ remapz = 0;

cint = [-1e10 caxs(1):sp:caxs(2) 1e10];

for ii=1:length(fields)
    subplot(1,length(fields),ii);
    tvec = dnum;
    if (length(fields{ii}{1}(:,1)) == length(Te))
        x = Te;
        Torp = 1;
    elseif (length(fields{ii}{1}(:,1)) == TL)
        x = T;
        Torp = 1;
    elseif (length(fields{ii}{1}(:,1)) == Vl)
        if (premap)
            x = mean(T_ofp,2);
            Torp = 1;
        else
            x = Vi/Vtot;
            Torp = 0;
        end
    end
    % Annual average:
    if (AA)
        var = squeeze(monmean(fields{ii}{1},2,ndays))'*scale;
    else
        var = fields{ii}{1}*scale;
    end
    yrvec = unique(dvec(:,1));
    % Subtract climatology:
    if (Sclim)
        var = var-repmat(mean(var(:,find(yrvec>=climean(1) & ...
                                       yrvec<=climean(2)),:),2),[1 length(yrvec)]);
    end
    % Time-integrate:
    if (Tint == 1)
        var = cumsum(var*86400*DT,1);
    end
    [X,Y] = ndgrid(yrvec,x);
    contourf(X,Y,var',cint,'linestyle','none');
    
    cb = colorbar('Location','NorthOutside','FontSize',15);    
    ylabel(cb,label);
    if (Torp)
        ylim([-3 31]);
        set(gca,'ytick',-5:5:35);
    else
        ylim([0 1]);
        set(gca,'ydir','reverse');
    end        
    grid on;
    caxis(caxs);
    if (ii == 1)
        if (Torp)
            ylabel('Temperature ($^\circ$C)');
        else
            ylabel('Volume fraction $p$');
        end
    end
    xlabel('Date');
    title(fields{ii}{2});
    set(gca,'FontSize',15);
end
colormap(redblue);


% $$$ %%% ARGO data:
% $$$ 
% $$$ load('ArgoData.mat');
% $$$ yrsargo = 2004:2014;
% $$$ Teargo = Tedge;
% $$$ Hargo = zeros(length(Teargo),length(yrsargo));
% $$$ Vargo = Hargo;
% $$$ dvec = datevec(days);
% $$$ for yi=1:length(yrsargo)
% $$$     yr = yrsargo(yi);
% $$$     Hargo(:,yi) = mean(Hsnap(:,dvec(:,1)==yr),2);
% $$$     Vargo(:,yi) = mean(Vsnap(:,dvec(:,1)==yr),2);
% $$$ end
% $$$ HIargo = Hargo-rho0*Cp*Vargo.*repmat(Teargo',[1 length(yrsargo)]);
% $$$ HEargo = Hargo-HIargo;
% $$$ 
% $$$ [X,Y] = ndgrid(yrsargo,Teargo);
% $$$ 
% $$$ Hargo = Hargo-repmat(mean(Hargo,2),[1 length(yrsargo)]);
% $$$ HEargo = HEargo-repmat(mean(HEargo,2),[1 length(yrsargo)]);
% $$$ HIargo = HIargo-repmat(mean(HIargo,2),[1 length(yrsargo)]);
% $$$ 
% $$$ % Fluxes:
% $$$ scale = 1/1e23;label = '(10$^{23}$J)';
% $$$ caxs = [-0.5 0.5];
% $$$ sp = 0.025;
% $$$ 
% $$$ cint = [-1e10 caxs(1):sp:caxs(2) 1e10];
% $$$ 
% $$$ fields = { ...
% $$$           {[], 'Internal Heat Content $\mathcal{H}_I$','m',2,'-'}, ...
% $$$           {[], 'Heat Content $\mathcal{H}$','m',2,'-'}, ...
% $$$           {[], 'External Heat Content $\mathcal{H}_E$','m',2,'-'}, ...
% $$$           };
% $$$ 
% $$$ figure;
% $$$ subplot(1,3,1);
% $$$ contourf(X,Y,HIargo'*scale,cint,'linestyle','none');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ set(gca,'ytick',-5:5:35);
% $$$ ylabel(cb,label);
% $$$ xlabel('Date');
% $$$ ylim([-3 31]);
% $$$ grid on;
% $$$ caxis(caxs);
% $$$ xlim([2004 2014]);
% $$$ title(fields{1}{2});
% $$$ set(gca,'FontSize',15);
% $$$ subplot(1,3,2);
% $$$ contourf(X,Y,Hargo'*scale,cint,'linestyle','none');
% $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ set(gca,'ytick',-5:5:35);
% $$$ ylabel(cb,label);
% $$$ xlabel('Date');
% $$$ ylim([-3 31]);
% $$$ grid on;
% $$$ caxis(caxs);
% $$$ xlim([2004 2014]);
% $$$ title(fields{2}{2});
% $$$ set(gca,'FontSize',15);
% $$$ subplot(1,3,3);
% $$$ contourf(X,Y,HEargo'*scale,cint,'linestyle','none');
% $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ set(gca,'ytick',-5:5:35);
% $$$ ylabel(cb,label);
% $$$ xlabel('Date');
% $$$ ylim([-3 31]);
% $$$ grid on;
% $$$ caxis(caxs);
% $$$ xlim([2004 2014]);
% $$$ title(fields{3}{2});
% $$$ set(gca,'FontSize',15);
% $$$ colormap(redblue);

%%%% Depth space plot:

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
    
    Vs = zeros(zL,tL,1);
    Hs = zeros(zL,tL,1);

    %% Load Global Budget:
    for i=1:length(outputs)
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        dnum = [dnum; time];

        load([base model sprintf('_output%03d_',outputs(i)) 'VHofz.mat']);
        Vs(:,:,i) = V;
        Hs(:,:,i) = H;
        ycur = ycur+nyrs;
    end
    months = [1:length(Vs(1,:,1))];
    yrs = [1:length(Vs(1,1,:))];
    
    % Correct time vector zero year:
    dvec = datevec(dnum);
    dvec(:,1) = dvec(:,1) + 1900;
    dnum = datenum(dvec);
    
    Ts = Hs/rho0/Cp./Vs;
    Vs = cumsum(Vs,1);
    Hs = cumsum(Hs,1);

% Subtract climatology:
Sclim = 1;
climean = [1976 1986];

% Annual average:
AA = 1;

tvec = dnum;
% Annual average:
if (AA)
    T =  squeeze(monmean(Ts,2,ndays));
    H =  squeeze(monmean(Hs,2,ndays));
    V =  squeeze(monmean(Vs,2,ndays));
    tvec = unique(dvec(:,1));
else
    V = Vs;
    H = Hs;
end
% Subtract climateology:
if (Sclim)
    T = T-repmat(mean(T(:,find(tvec>=climean(1) & tvec<=climean(2))),2),[1 ...
                        length(tvec)]);
    V = V-repmat(mean(V(:,find(tvec>=climean(1) & tvec<=climean(2))),2),[1 ...
                        length(tvec)]);
    H = H-repmat(mean(H(:,find(tvec>=climean(1) & tvec<=climean(2))),2),[1 ...
                        length(tvec)]);
end
[X,Y] = ndgrid(tvec,-z);
% $$$ subplot(1,2,1);
% $$$ pcolPlot(X,Y,H');
% $$$ caxs = [-0.8 0.8]*1e23;
% $$$ sp = 0.005e23
% $$$ cint = [-1e50 caxs(1):sp:caxs(2) 1e50];
% $$$ contourf(X,Y,H',cint,'linestyle','none');
% $$$ title('Heat Content Anomaly above depth level (J)');
caxs = [-0.15 0.15];
sp = 0.005
cint = [-1e50 caxs(1):sp:caxs(2) 1e50];
contourf(X,Y,T',cint,'linestyle','none');
cb = colorbar('Location','NorthOutside','FontSize',15);    
ylabel('Depth (m)');
xlabel('Date');
title('Temperature Anomaly $(^\circ$C)');%Heat Content Anomaly above depth level (J)');
caxis(caxs);%[-0.5 0.5]*1e23);
set(gca,'FontSize',15);
colormap(redblue);
