% This script makes temperature-time global plots of the heat budget
% in MOM5.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = struct( ...
       'model',{'ACCESS-OM2_1deg_jra55_rdf','ACCESS-OM2_1deg_jra55_rdf_pert'},...
       'outputs',[51:55], ...
       'zeroyear',[1972]);

rr = 2;
for rr = 1:length(RUNS)
    rr
    outputs = RUNS(rr).outputs;
    model = RUNS(rr).model;

% $$$     clearvars -except base RUNS rr outputs model leg legh;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    RUNS(rr).ndays = ndays;
    region = 'Global';
    if (ndays(1) < 300) % Monthly data
        nyrs = tL/12;
        nmnt = 12;
        szTe = [TL+1 nmnt nyrs];szT  = [TL nmnt nyrs];
        yrs = 1:nyrs;
    else
        nyrs = tL;
        nmnt = 1;
        szTe = [TL+1 nyrs];szT = [TL nyrs];
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

        load([base model sprintf('_output%03d_',outputs(i)) region '_HBud.mat']);
        
        % Fluxes:
        RUNS(rr).P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
        RUNS(rr).F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
        RUNS(rr).Ffz(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.FRZ,szTe); % Surface heat flux (W)
        RUNS(rr).Fsw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Surface heat flux (W)
        RUNS(rr).Fsh(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDS,szTe); % Surface heat flux (W)
        RUNS(rr).Fet(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ETS,szTe); % Surface heat flux (W)
        RUNS(rr).MNL(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
        RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
        if (isfield(GWB,'VDFkppiw')) % Vertical mixing components
            RUNS(rr).VDFkppiw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw,szTe);
            RUNS(rr).VDFkppish(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppish,szTe);
            RUNS(rr).VDFkppicon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppicon,szTe);
            RUNS(rr).VDFkppbl(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppbl,szTe);
            RUNS(rr).VDFkppdd(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppdd,szTe);
            RUNS(rr).VDFwave(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFwave,szTe);
            RUNS(rr).VDFnloc(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
            RUNS(rr).VDFsum(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw+GWB.VDFkppish+GWB.VDFkppicon+ ...
                GWB.VDFkppbl+GWB.VDFkppdd+GWB.VDFwave+GWB.KNL,szTe);
            % Note: May be missing enhanced mixing near rivers
            % (river_diffuse_temp) in ACCESS-OM2
        end
        if (isfield(GWB,'RED')) % Redi Diffusion
            RUNS(rr).R(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.RED+GWB.K33,szTe); % Redi diffusion (W)
        else
            RUNS(rr).R(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end
        if (isfield(GWB,'NGM')) % GM parameterization
            RUNS(rr).GM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NGM,szTe); % GM (W)
        else
            RUNS(rr).GM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'MDS')) % Mix-downslope
            RUNS(rr).MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);
            RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
        else
            RUNS(rr).MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'SIG')) % Sigma-diff
            RUNS(rr).SG(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SIG,szTe);
            RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.SIG,szTe);
        else
            RUNS(rr).SG(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
            RUNS(rr).NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
        else
            RUNS(rr).NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end    
        if (isfield(GWB,'SUB'))
            RUNS(rr).SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);
        else
            RUNS(rr).SUB(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end
        RUNS(rr).D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV,szTe)-RUNS(rr).GM(:,:,ycur:(ycur+nyrs-1))-RUNS(rr).SUB(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
        RUNS(rr).TEN(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe); % Tendency
        RUNS(rr).ADV(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.ADV,szTe); % Advection
        RUNS(rr).SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
        RUNS(rr).JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux

        % Pacific Interior fluxes:
        if (strcmp(region,'Pacific'))
            RUNS(rr).JI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.JBS+GWB.JSP+GWB.JITF,szTe); %Combined volume flux out
            RUNS(rr).QI(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.QBS+GWB.QSP+GWB.QITF,szTe); %Combined heat flux out
        else
            RUNS(rr).QI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
            RUNS(rr).JI(:,:,ycur:(ycur+nyrs-1)) = zeros(size(RUNS(rr).P(:,:,ycur:(ycur+nyrs-1))));
        end

        % Snapshot fields:
        RUNS(rr).dVdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dVdt,szTe); % V Change (m3s-1)
        RUNS(rr).dHdt(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.dHdt,szTe); % H Change (W)

        % Water-mass transformation:
        RUNS(rr).G(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).dVdt(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).JS(:,:,ycur:(ycur+nyrs-1)) + RUNS(rr).JI(:,:,ycur:(ycur+nyrs-1)); %Water-mass transformation (m3s-1)

        % Surface Volume flux base flux (not P!)
        RUNS(rr).JSH(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 nmnt nyrs])*rho0*Cp;

        % Interior heat source P:
        RUNS(rr).PI(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).P(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).JSH(:,:,ycur:(ycur+nyrs-1));

        % Interior heat source Q:
        RUNS(rr).QII(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).QI(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 nmnt nyrs])*rho0*Cp;

        % Across-isotherm advective heat flux:
        RUNS(rr).CIA(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 nmnt nyrs])*rho0*Cp;

        % External HC Tendency:
        RUNS(rr).EHC(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 nmnt nyrs])*rho0*Cp;

        % Internal HC Tendency:
        RUNS(rr).N(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).dHdt(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).EHC(:,:,ycur:(ycur+nyrs-1));
% $$$ N(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe);

        % Implicit mixing:
        RUNS(rr).I(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).N(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).F(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).P(:,:,ycur:(ycur+nyrs-1)) ...
                                  - RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).R(:,:,ycur:(ycur+nyrs-1)) + RUNS(rr).JSH(:,:,ycur:(ycur+nyrs-1)) ...
                                  - RUNS(rr).SUB(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).GM(:,:,ycur:(ycur+nyrs-1));

        % Non-advective flux into volume:
        RUNS(rr).B(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).F(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).M(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).I(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).R(:,:,ycur:(ycur+nyrs-1));

        % Monthly binned Internal HC Tendency:
        if (isfield('GWB','TENMON'))
            RUNS(rr).Nmon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TENMON,szTe);
        end

        % WMT from B:
        RUNS(rr).WMTM(:,:,ycur:(ycur+nyrs-1)) = -diff(RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        RUNS(rr).WMTF(:,:,ycur:(ycur+nyrs-1)) = -diff(RUNS(rr).F(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        RUNS(rr).WMTI(:,:,ycur:(ycur+nyrs-1)) = -diff(RUNS(rr).I(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        RUNS(rr).WMTR(:,:,ycur:(ycur+nyrs-1)) = -diff(RUNS(rr).R(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
        RUNS(rr).WMT(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).WMTM(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).WMTF(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).WMTI(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).WMTR(:,:,ycur:(ycur+nyrs-1));

        % WMT HB from B:
        RUNS(rr).HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*RUNS(rr).WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 nmnt nyrs]);
        RUNS(rr).HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*RUNS(rr).WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 nmnt nyrs]);
        RUNS(rr).HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*RUNS(rr).WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 nmnt nyrs]);
        RUNS(rr).HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*RUNS(rr).WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 nmnt nyrs]);
        RUNS(rr).HWMT(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).HWMTM(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).HWMTF(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).HWMTI(:,:,ycur:(ycur+nyrs-1))+RUNS(rr).HWMTR(:,:,ycur:(ycur+nyrs-1));

        ycur = ycur+nyrs;
    end
% $$$     months = [1:length(P(1,:,1))];
% $$$     yrs = [1:length(P(1,1,:))];
    
    % Correct time vector zero year:
    dvec = datevec(dnum);
    if (isfield(RUNS(1),'zeroyear'))
        dvec(:,1) = dvec(:,1)-dvec(1,1) + RUNS(rr).zeroyear;
    else
        dvec(:,1) = dvec(:,1) + 1900;
    end
    RUNS(rr).dvec = dvec;
    RUNS(rr).dnum = datenum(dvec);
    RUNS(rr).yrs = 1:10;

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
    RUNS(rr).V = cat(1,cumsum(V,1,'reverse'),zeros(1,nyrs));
    RUNS(rr).H = cat(1,cumsum(H,1,'reverse'),zeros(1,nyrs));
    RUNS(rr).HE = rho0*Cp*RUNS(rr).V.*repmat(Te,[1 nyrs]);
    RUNS(rr).HI = RUNS(rr).H - RUNS(rr).HE;    
    
% $$$     % Remap to ocean percentile:
% $$$     Vtot = V(1,:);
% $$$     p_ofT = V./repmat(Vtot,[TL+1 1]);
% $$$ 
% $$$     p = linspace(0,1,1000);
% $$$     pl = length(p);
% $$$     T_ofp = zeros(pl,length(yrs));
% $$$     H_ofp = T_ofp;
% $$$     HI_ofp = T_ofp;
% $$$     HE_ofp = T_ofp;
% $$$     for yi = 1:length(yrs)
% $$$         T_ofp(:,yi) = interp1(p_ofT(:,yi)+(1:TL+1)'/1e10,Te,p,'linear');
% $$$         T_ofp(1,yi) = Te(end);
% $$$         H_ofp(:,yi) = interp1(Te,H(:,yi),T_ofp(:,yi),'linear');
% $$$         HE_ofp(:,yi) = interp1(Te,HE(:,yi),T_ofp(:,yi),'linear');
% $$$         HI_ofp(:,yi) = interp1(Te,HI(:,yi),T_ofp(:,yi),'linear');
% $$$     end

    % Remap to ocean volume:
    Vtot = max(RUNS(rr).V(1,:));
% $$$     p_ofT = V./repmat(Vtot,[TL+1 1]);

    Vi = linspace(0,Vtot,1000);
    Vl = length(Vi);
    T_ofV = zeros(Vl,nyrs);%length(yrs));
    H_ofV = T_ofV;
    HI_ofV = T_ofV;
    HE_ofV = T_ofV;
    for yi = 1:nyrs
        T_ofV(:,yi) = interp1(flipud(RUNS(rr).V(:,yi))+(1:(TL+1))'/(TL+1)*Vtot/1e12,flipud(Te),Vi,'linear');
        T_ofV(1,yi) = Te(end);
        H_ofV(:,yi) = interp1(Te,RUNS(rr).H(:,yi),T_ofV(:,yi),'linear');
        HE_ofV(:,yi) = interp1(Te,RUNS(rr).HE(:,yi),T_ofV(:,yi),'linear');
        HI_ofV(:,yi) = interp1(Te,RUNS(rr).HI(:,yi),T_ofV(:,yi),'linear');
    end
    RUNS(rr).Vi = Vi;
    RUNS(rr).Vl = Vl;
    RUNS(rr).T_ofV = T_ofV;
    RUNS(rr).H_ofV = H_ofV;
    RUNS(rr).HI_ofV = HI_ofV;
    RUNS(rr).HE_ofV = HE_ofV;
    RUNS(rr).Vtot = Vtot;
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

% Take annual mean to make things easier:
for rr = 1:length(RUNS)
    names = fieldnames(RUNS(rr));
    for vi = 1:length(names)
        eval(['sz = size(RUNS(rr).' names{vi} ');']);
        ind = find(sz==12);
        if (length(ind)==1)
            eval(['RUNS(rr).' names{vi} '=squeeze(monmean(RUNS(rr).' names{vi} ...
                  ',ind,RUNS(rr).ndays(1:12)));']);
        end
    end
end
          
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
% $$$           {HI, '$\mathcal{H}_I(\Theta)$','m',2,'-'}, ...
% $$$           {H, '$\mathcal{H}(\Theta)$','m',2,'-'}, ...
% $$$           {HE, '$\mathcal{H}_E(\Theta)$','m',2,'-'}, ...
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

% Heat content:
scale = 1/1e23;label = '$10^{23}J$';x = Te;
caxs = [-2 2];
sp = 0.01;

% Time-integrate fluxes:
Tint = 0;

% remapping for percentiles:
premap = 1;

% Subtract climatology:
Sclim = 0;
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
         {'ACCESS-OM2_1deg_jra55_rdf',[51:55],[1972]}, ...
         {'ACCESS-OM2_1deg_jra55_rdf_pert',[51:55],[1972]}, ...
       };

rr = 2;
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
    
    Vs = zeros(zL,nyrs,1);
    Hs = zeros(zL,nyrs,1);

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
    if (length(RUNS{rr}{3})== 1)
        dvec(:,1) = dvec(:,1)-dvec(1,1) + RUNS{rr}{3}(1);
    else
        dvec(:,1) = dvec(:,1) + 1900;
    end
    dnum = datenum(dvec);
    yrs = 1:10;
    
    Ts = Hs/rho0/Cp./Vs;
    Vs = cumsum(Vs,1);
    Hs = cumsum(Hs,1);

% Subtract climatology:
Sclim = 0;
climean = [1972 1981];

% Annual average:
AA = 0;

% Annual average:
if (AA)
    T =  squeeze(monmean(Ts,2,ndays));
    H =  squeeze(monmean(Hs,2,ndays));
    V =  squeeze(monmean(Vs,2,ndays));
    tvec = unique(dvec(:,1));
else
    V = reshape(Vs,[length(Vs(:,1,1))  nyrs*length(Vs(1,1,:))]);
    H = reshape(Hs,[length(Vs(:,1,1)) nyrs*length(Hs(1,1,:))]);
    T = reshape(Ts,[length(Vs(:,1,1)) nyrs*length(Hs(1,1,:))]);
end

% $$$ Tcont = T;
T = T - Tcont;

% Subtract climateology:
yrvec = unique(dvec(:,1));
% Subtract climatology:
if (Sclim)
    T = T-repmat(mean(T(:,find(yrvec>=climean(1) & ...
                                     yrvec<=climean(2)),:),2),[1 length(yrvec)]);
    V = V-repmat(mean(V(:,find(yrvec>=climean(1) & ...
                                     yrvec<=climean(2)),:),2),[1 length(yrvec)]);
    H = H-repmat(mean(H(:,find(yrvec>=climean(1) & ...
                                     yrvec<=climean(2)),:),2),[1 length(yrvec)]);
end
[X,Y] = ndgrid(yrvec,-z);
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
xlim([1972 2018]);
