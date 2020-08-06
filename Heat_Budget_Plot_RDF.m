% This script makes temperature-time global plots of the heat budget
% in MOM5.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = struct( ...
       'model',{'ACCESS-OM2_1deg_jra55_rdf','ACCESS-OM2_1deg_jra55_rdf_pert'},...
       'outputs',[51:55], ...
       'zeroyear',[1972]);

% Determine absolute maximum volume to use for grid for both RUNS:
Vtot = 0;
for rr = 1:length(RUNS)
    outputs = RUNS(rr).outputs;
    model = RUNS(rr).model;    
    for i=1:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) ...
              '_VHza.mat']);
        Vtot = max(Vtot,max(nansum(nansum(V,1),2)));
    end
end

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
    region = 'Global';
    if (ndays(1) < 300) % Monthly data
        nyrs = tL/12;
        nmnt = 12;
        szTe = [TL+1 nmnt nyrs];szT  = [TL nmnt nyrs];
        szZ = [zL nmnt nyrs];
        yrs = 1:nyrs;
    else
        nyrs = tL;
        nmnt = 1;
        szTe = [TL+1 nyrs];szT = [TL nyrs];
        szZ = [zL nyrs];
    end    
    ycur = 1;
    
    dnum = [];
    dnum_snap = [];
    DT = [];

    %% Load Global Budget:
    for i=1:length(outputs)
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        dnum = [dnum; time];
        dnum_snap = [dnum_snap; time_snap];
        DT = [DT; ndays];

        load([base model sprintf('_output%03d_',outputs(i)) region '_HBud.mat']);

        % Replace by monthly binning data:
        load([base model sprintf('_output%03d_',outputs(i)) region '_HBud_MonAnBin.mat']);
        names = fieldnames(GWBmon);
        for vi=1:length(names)
            eval(['GWB.' names{vi} ' = GWBmon.' names{vi} ';']);
        end
        
        % Fluxes:
        RUNS(rr).P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
        RUNS(rr).F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
        RUNS(rr).SWP(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe);
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

        % Implicit mixing (includes numerical mixing from SUB and GM!!):
        RUNS(rr).I(:,:,ycur:(ycur+nyrs-1)) = RUNS(rr).N(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).F(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).P(:,:,ycur:(ycur+nyrs-1)) ...
                                  - RUNS(rr).M(:,:,ycur:(ycur+nyrs-1)) - RUNS(rr).R(:,:,ycur:(ycur+nyrs-1)) + RUNS(rr).JSH(:,:,ycur:(ycur+nyrs-1));

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
    
    % Total surface forcing, mixing and tendency:
    RUNS(rr).SUF = RUNS(rr).F+RUNS(rr).PI;
    RUNS(rr).MIX = RUNS(rr).I+RUNS(rr).M+RUNS(rr).R;
    
    % Timing:
    dvec = datevec(dnum);
    dvec_snap = datevec(dnum_snap);
    if (isfield(RUNS(1),'zeroyear'))
        dvec(:,1) = dvec(:,1)-dvec(1,1) + RUNS(rr).zeroyear;
        dvec_snap(:,1) = dvec_snap(:,1)-dvec_snap(1,1) + RUNS(rr).zeroyear;
    else
        dvec(:,1) = dvec(:,1) + 1900;
        dvec_snap(:,1) = dvec_snap(:,1) + 1900;
    end
    RUNS(rr).dvec = dvec;
    RUNS(rr).dnum = datenum(dvec);
    RUNS(rr).dvec_snap = dvec_snap;
    RUNS(rr).dnum_snap = datenum(dvec_snap);
    RUNS(rr).DT = reshape(DT,[nmnt nyrs*length(outputs)]);

    % Load Global V and H from averages:
    ycur = 1;
    for i=1:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) ...
              '_VHza.mat']);
        RUNS(rr).V(:,:,ycur:(ycur+nyrs-1)) = reshape(squeeze(nansum(V,1)),szT);
        RUNS(rr).H(:,:,ycur:(ycur+nyrs-1)) = reshape(squeeze(nansum(H,1)),szT);
        ycur = ycur+nyrs;
    end
    sz = size(RUNS(rr).V);
    sz(1) = 1;
    RUNS(rr).V = cat(1,cumsum(RUNS(rr).V,1,'reverse'),zeros(sz));
    RUNS(rr).H = cat(1,cumsum(RUNS(rr).H,1,'reverse'),zeros(sz));
    RUNS(rr).HE = rho0*Cp*RUNS(rr).V.*repmat(Te,sz);%[1 nyrs]);
    RUNS(rr).HI = RUNS(rr).H - RUNS(rr).HE;    

    % Load Global V and H from snapshots:
    for i=1:length(outputs)
        load([base model sprintf('_output%03d',outputs(i)) '_VHzaSNAP.mat']);
        if (i == 1)
            RUNS(rr).Vsnap = squeeze(nansum(V,1));
            RUNS(rr).Hsnap = squeeze(nansum(H,1));
        else
            RUNS(rr).Vsnap = cat(2,RUNS(rr).Vsnap,squeeze(nansum(V(:,:,2:end),1)));
            RUNS(rr).Hsnap = cat(2,RUNS(rr).Hsnap,squeeze(nansum(H(:,:,2:end),1)));
        end
    end
    sz = size(RUNS(rr).Vsnap);
    sz(1) = 1;
    RUNS(rr).Vsnap = cat(1,cumsum(RUNS(rr).Vsnap,1,'reverse'),zeros(sz));
    RUNS(rr).Hsnap = cat(1,cumsum(RUNS(rr).Hsnap,1,'reverse'),zeros(sz));
    RUNS(rr).HEsnap = rho0*Cp*RUNS(rr).Vsnap.*repmat(Te,sz);%[1 nyrs]);
    RUNS(rr).HIsnap = RUNS(rr).Hsnap - RUNS(rr).HEsnap;    

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
% $$$     Vtot = max(max(RUNS(rr).V(1,:,:)));
    Vi = linspace(0,Vtot,2000);
    Vl = length(Vi);
    sz = size(RUNS(rr).V);
    sz_snap = size(RUNS(rr).Vsnap);
    sz(1) = Vl;sz_snap(1) = Vl;
    T_ofV = zeros(sz);
    Tsnap_ofV = zeros(sz_snap);
    vars = {'H','HE','HI','SUF','MIX','N','dHdt','EHC','SWP','I','M','R'};
    varsSNAP = {'Hsnap','HEsnap','HIsnap'};
    for vi = 1:length(vars)
        eval(['RUNS(rr).' vars{vi} '_ofV = T_ofV;']);
    end
    for vi = 1:length(varsSNAP)
        eval(['RUNS(rr).' varsSNAP{vi} '_ofV = Tsnap_ofV;']);
    end
    for mi = 1:sz(2)
    for yi = 1:sz(3)
        T_ofV(:,mi,yi) = interp1(flipud(RUNS(rr).V(:,mi,yi))+(0:TL)'/TL*Vtot/1e12,flipud(Te),Vi,'pchip',NaN);
        T_ofV(1,mi,yi) = Te(end);
        for vi=1:length(vars)
            eval(['RUNS(rr).' vars{vi} '_ofV(:,mi,yi) = interp1(Te,RUNS(rr).' ...
                  vars{vi} '(:,mi,yi),T_ofV(:,mi,yi),''pchip'',NaN);']);
            eval(['RUNS(rr).' vars{vi} '_ofV(end,mi,yi) = RUNS(rr).' ...
                  vars{vi} '(1,mi,yi);']);
        end
    end
    end
    for mi = 1:sz_snap(2)
        Tsnap_ofV(:,mi) = interp1(flipud(RUNS(rr).Vsnap(:,mi))+(0:TL)'/TL*Vtot/1e12,flipud(Te),Vi,'pchip',NaN);
        Tsnap_ofV(1,mi) = Te(end);
        for vi=1:length(varsSNAP)
            eval(['RUNS(rr).' varsSNAP{vi} '_ofV(:,mi) = interp1(Te,RUNS(rr).' ...
                  varsSNAP{vi} '(:,mi),Tsnap_ofV(:,mi),''pchip'',NaN);']);
            eval(['RUNS(rr).' varsSNAP{vi} '_ofV(end,mi) = RUNS(rr).' ...
                  varsSNAP{vi} '(1,mi);']);
        end
    end
    RUNS(rr).T_ofV = T_ofV;
    RUNS(rr).Tsnap_ofV = Tsnap_ofV;
    RUNS(rr).Vtot = Vtot;
    
    % Load depth space data:
    ycur = 1;
    for i=1:length(outputs)
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        load([base model sprintf('_output%03d_',outputs(i)) 'VHofz.mat']);
        RUNS(rr).V_ofz(:,:,ycur:(ycur+nyrs-1)) = reshape(Vz,szZ);
        RUNS(rr).H_ofz(:,:,ycur:(ycur+nyrs-1)) = reshape(Hz,szZ);
        ycur = ycur+nyrs;
    end
    RUNS(rr).T_ofz = RUNS(rr).H_ofz/rho0/Cp./RUNS(rr).V_ofz;
    RUNS(rr).V_ofz = cumsum(RUNS(rr).V_ofz,1);
    RUNS(rr).H_ofz = cumsum(RUNS(rr).H_ofz,1);
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
        if (ind == 2)
            eval(['var = RUNS(rr).' names{vi} ';']);
            var = squeeze(sum(var.*repmat(permute( ...
                RUNS(rr).DT./repmat(sum(RUNS(rr).DT,1),[12 1]) ...
                ,[3 1 2]),[sz(1) 1 1]),2));
            eval(['RUNS(rr).' names{vi} '=var;']);
        end
    end
    dnum = reshape(RUNS(rr).dnum,[12 length(RUNS(rr).dnum)/12]);
    RUNS(rr).dnum = squeeze(sum(dnum.*RUNS(rr).DT./repmat(sum(RUNS(rr).DT,1),[12 1]),1));
    RUNS(rr).dvec = datevec(RUNS(rr).dnum);
    RUNS(rr).dnum_snap = cat(1,RUNS(rr).dnum_snap(1)-RUNS(rr).DT(1,1),RUNS(rr).dnum_snap(12:12:end));
    RUNS(rr).dvec_snap = datevec(RUNS(rr).dnum_snap);
    RUNS(rr).DT = squeeze(sum(RUNS(rr).DT,1));
    RUNS(rr).Vsnap = RUNS(rr).Vsnap(:,1:12:end);
    RUNS(rr).Hsnap = RUNS(rr).Hsnap(:,1:12:end);
    RUNS(rr).HEsnap = RUNS(rr).HEsnap(:,1:12:end);
    RUNS(rr).HIsnap = RUNS(rr).HIsnap(:,1:12:end);
    RUNS(rr).Hsnap_ofV = RUNS(rr).Hsnap_ofV(:,1:12:end);
    RUNS(rr).HEsnap_ofV = RUNS(rr).HEsnap_ofV(:,1:12:end);
    RUNS(rr).HIsnap_ofV = RUNS(rr).HIsnap_ofV(:,1:12:end);
    RUNS(rr).Tsnap_ofV = RUNS(rr).Tsnap_ofV(:,1:12:end);
end

% Define some timings and drop the last years:
tinds = find(RUNS(2).dvec(:,1) <= 2017);
tinds_snap = [tinds; tinds(end)+1];
tvec = RUNS(2).dnum(tinds)';
tvecSNAP = RUNS(2).dnum_snap(tinds_snap);
[X,Y] = ndgrid(tvec,Vi/Vi(end));

% $$$ %%% Plot Theta(V,t) and Theta(z,t) time series:
% $$$ 
% $$$ % Linear fit control:
% $$$ ThetaVlin = linfit(tvec,RUNS(1).T_ofV(:,tinds)');
% $$$ ThetaZlin = linfit(tvec,RUNS(1).T_ofz(:,tinds)');
% $$$ 
% $$$ % $$$ % Plot/check a couple points:
% $$$ % $$$ pts = [50];
% $$$ % $$$ figure;
% $$$ % $$$ for i=1:length(pts)
% $$$ % $$$     plot(tvec,ThetaV(:,pts(i)),'-k');
% $$$ % $$$     hold on;
% $$$ % $$$     plot(tvec,ThetaVlin(:,pts(i)),'-r');
% $$$ % $$$ end
% $$$ 
% $$$ % Calculate anomalies:
% $$$ ThetaV = RUNS(2).T_ofV(:,tinds)'-ThetaVlin;
% $$$ ThetaZ = RUNS(2).T_ofz(:,tinds)'-ThetaZlin;
% $$$ [Xz,Yz] = ndgrid(tvec,z);
% $$$ [XT,YT] = ndgrid(tvec,mean(RUNS(1).T_ofV(:,tinds),2));
% $$$ 
% $$$ caxs = [-0.15 0.15];
% $$$ cpts = [-1 -1:0.0025:1 1];
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ subplot(1,3,1);
% $$$ contourf(XT,YT,ThetaV,cpts,'linestyle','none');
% $$$ caxis([-0.3 0.3]);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\Theta(\mathcal{V},t)$');
% $$$ %ylim(cb,[-0.15 0.15001]);
% $$$ xlabel('Year');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,2);
% $$$ contourf(X,Y,ThetaV,cpts,'linestyle','none');
% $$$ caxis(caxs);%[-0.2 0.2]);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\Theta(\mathcal{V},t)$');
% $$$ ylim(cb,[-0.15 0.15001]);
% $$$ set(gca,'ydir','reverse');
% $$$ xlabel('Year');
% $$$ ylabel('Volume Fraction');
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,3);
% $$$ contourf(Xz,Yz,ThetaZ,cpts,'linestyle','none');
% $$$ caxis(caxs);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','NorthOutside');
% $$$ ylim(cb,[-0.15 0.15001]);
% $$$ ylabel(cb,'$\Theta(z,t)$');
% $$$ set(gca,'ydir','reverse');
% $$$ xlabel('Year');
% $$$ ylabel('Depth (m)');%Volume (%)')
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ 
% $$$ %%% Surface forcing, mixing and heat content time series's:
% $$$ 
% $$$ % Check budget:
% $$$ figure;
% $$$ plot(Y(1,:),mean(RUNS(1).SUF_ofV,2),'-k');
% $$$ hold on;
% $$$ plot(Y(1,:),mean(RUNS(1).MIX_ofV,2),'-b');
% $$$ plot(Y(1,:),mean(RUNS(1).N_ofV,2),'-m');
% $$$ plot(Y(1,:),mean(RUNS(2).SUF_ofV,2),'--k');
% $$$ plot(Y(1,:),mean(RUNS(2).MIX_ofV,2),'--b');
% $$$ plot(Y(1,:),mean(RUNS(2).N_ofV,2),'--m');
% $$$ title('Average budgets in Vol space');
% $$$ legend('F','M','dH/dt');%RDF_pert');

% Load vars simple:
Hc = RUNS(1).Hsnap_ofV(:,tinds_snap)';
Hp = RUNS(2).Hsnap_ofV(:,tinds_snap)';
Hcavg = RUNS(1).H_ofV(:,tinds)';
Hpavg = RUNS(2).H_ofV(:,tinds)';

Fc = RUNS(1).SUF_ofV(:,tinds)';
Mc = RUNS(1).MIX_ofV(:,tinds)';
Nc = RUNS(1).N_ofV(:,tinds)';
Sc = RUNS(1).SWP_ofV(:,tinds)';
Fp = RUNS(2).SUF_ofV(:,tinds)';
Mp = RUNS(2).MIX_ofV(:,tinds)';
Np = RUNS(2).N_ofV(:,tinds)';
Sp = RUNS(2).SWP_ofV(:,tinds)';

MIc = RUNS(1).I_ofV(:,tinds)';
MVc = RUNS(1).M_ofV(:,tinds)';
MRc = RUNS(1).R_ofV(:,tinds)';
MIp = RUNS(2).I_ofV(:,tinds)';
MVp = RUNS(2).M_ofV(:,tinds)';
MRp = RUNS(2).R_ofV(:,tinds)';

DTc = RUNS(1).DT(tinds)'*86400;
DTp = RUNS(2).DT(tinds)'*86400;
Hc_fromN = cumsumt(Hc(1,:),Nc,DTc);
Hc_fromF = cumsumt(Hc(1,:),Fc,DTc);
Hc_fromM = cumsumt(Hc(1,:),Mc,DTc);
Hc_fromS = cumsumt(Hc(1,:),Sc,DTc);
Hp_fromN = cumsumt(Hp(1,:),Np,DTp);
Hp_fromF = cumsumt(Hp(1,:),Fp,DTp);
Hp_fromM = cumsumt(Hp(1,:),Mp,DTp);
Hp_fromS = cumsumt(Hp(1,:),Sp,DTp);

% Fits to control:
% linear regression:
[Hl,b] = linfit(tvecSNAP,Hc);
xvec = [ones(length(tvec),1) tvec];
Hlavg = xvec*b; % linear fit to Hcavg with same coefficients.

% end-to-end:
tsecs = sum(DTc);
m = (Hc(end,:)-Hc(1,:))/tsecs;
Hee = repmat(Hc(1,:),[length(tvecSNAP) 1])+ ...
     repmat(m,[length(tvecSNAP) 1]).* ...
     repmat(linspace(0,tsecs,length(tvecSNAP))',[1 Vl]);

% $$$ % Check budget closures and total heat before calculating perturbations:
% $$$ figure;
% $$$ inds = [500 250];
% $$$ for ii=1:length(inds)
% $$$     subplot(length(inds),1,ii);
% $$$     ind = inds(ii);
% $$$     plot(tvecSNAP,Hc(:,ind),'-k','linewidth',2);
% $$$     hold on;
% $$$     plot(tvecSNAP,Hc_fromN(:,ind),'--m','linewidth',2);
% $$$     plot(tvec,Hcavg(:,ind),'--k','linewidth',1);
% $$$     plot(tvecSNAP,Hl(:,ind),'-b','linewidth',2);
% $$$ % $$$     plot(tvecSNAP,Hl_fromN(:,ind)+Hl(1,ind),'--b','linewidth',2);
% $$$     plot(tvecSNAP,Hee(:,ind),'-r','linewidth',2);
% $$$     plot(tvecSNAP,Hp(:,ind),'-k','linewidth',2);
% $$$     plot(tvecSNAP,Hp_fromN(:,ind),'--m','linewidth',2);
% $$$     plot(tvec,Hpavg(:,ind),'--k','linewidth',1);
% $$$     legend('$H(\mathcal{V},t_{snap})$','Integrated dH/dt','$H(\mathcal{V},t)$','Linear Fit','End-to-end fit');
% $$$     title(['Heat Content above Vol Frac = ' num2str(ind/Vl)]);
% $$$     xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$     datetick(gca,'x','YYYY','keeplimits');
% $$$ end
% $$$ 
% $$$ 
% $$$ % $$$ %% Option 1 - no fitting to control:
% $$$ % $$$ % $$$ Hp = Hp - Hc;
% $$$ % $$$ % $$$ Hpavg = Hpavg - Hcavg;
% $$$ % $$$ % $$$ Fp = Fp - Fc;
% $$$ % $$$ % $$$ Np = Np - Nc;
% $$$ % $$$ % $$$ Mp = Mp - Mc;
% $$$ % $$$ % $$$ Sp = Sp - Sc;
% $$$ % $$$ Hp_fromN = cumsumt(Hp(1,:),Np,DTp);
% $$$ % $$$ Hp_fromF = cumsumt(Hp(1,:),Fp,DTp);
% $$$ % $$$ Hp_fromM = cumsumt(Hp(1,:),Mp,DTp);
% $$$ % $$$ Hp_fromS = cumsumt(Hp(1,:),Sp,DTp);

%% Option 2 - linear fit to time-integrated values (This is the
%% best one to use!):
Hd = Hp - Hl;

[Hl_fromN,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),Nc,DTc));
[Hl_fromF,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),Fc,DTc));
[Hl_fromM,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),Mc,DTc));
[Hl_fromS,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),Sc,DTc));
[Hl_fromMI,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),MIc,DTc));
[Hl_fromMV,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),MVc,DTc));
[Hl_fromMR,b] = linfit(tvecSNAP,cumsumt(zeros(size(Hc(1,:))),MRc,DTc));
Hl_fromN = Hl_fromN-repmat(Hl_fromN(1,:),[length(tvecSNAP) 1]);
Hl_fromF = Hl_fromF-repmat(Hl_fromF(1,:),[length(tvecSNAP) 1]);
Hl_fromM = Hl_fromM-repmat(Hl_fromM(1,:),[length(tvecSNAP) 1]);
Hl_fromS = Hl_fromS-repmat(Hl_fromS(1,:),[length(tvecSNAP) 1]);
Hl_fromMI = Hl_fromMI-repmat(Hl_fromMI(1,:),[length(tvecSNAP) 1]);
Hl_fromMV = Hl_fromMV-repmat(Hl_fromMV(1,:),[length(tvecSNAP) 1]);
Hl_fromMR = Hl_fromMR-repmat(Hl_fromMR(1,:),[length(tvecSNAP) 1]);

Hd_fromN = cumsumt(Hd(1,:),Np,DTp)-Hl_fromN;
Hd_fromF = cumsumt(Hd(1,:),Fp,DTp)-Hl_fromF;
Hd_fromM = cumsumt(Hd(1,:),Mp,DTp)-Hl_fromM;
Hd_fromS = cumsumt(Hd(1,:),Sp,DTp)-Hl_fromS;
Hd_fromMI = cumsumt(Hd(1,:),MIp,DTp)-Hl_fromMI;
Hd_fromMV = cumsumt(Hd(1,:),MVp,DTp)-Hl_fromMV;
Hd_fromMR = cumsumt(Hd(1,:),MRp,DTp)-Hl_fromMR;

% $$$ % $$$ %% Option 3: End-to-end difference instead of linear fit:
% $$$ % $$$ Hd = Hp - Hee;
% $$$ % $$$ Nd = Np - repmat(sum(Nc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Fd = Fp - repmat(sum(Fc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Md = Mp - repmat(sum(Mc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Sd = Sp - repmat(sum(Sc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Hd_fromN = cumsumt(Hd(1,:),Nd,DTp);
% $$$ % $$$ Hd_fromF = cumsumt(Hd(1,:),Fd,DTp);
% $$$ % $$$ Hd_fromM = cumsumt(Hd(1,:),Md,DTp);
% $$$ % $$$ Hd_fromS = cumsumt(Hd(1,:),Sd,DTp);
% $$$ % $$$ 
% $$$ % $$$ MId = MIp - repmat(sum(MIc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Hd_fromMI = cumsumt(Hd(1,:),MId,DTp);
% $$$ % $$$ MVd = MVp - repmat(sum(MVc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Hd_fromMV = cumsumt(Hd(1,:),MVd,DTp);
% $$$ % $$$ MRd = MRp - repmat(sum(MRc.*repmat(DTc/sum(DTc),[1 Vl]),1),[length(tvec) 1]);
% $$$ % $$$ Hd_fromMR = cumsumt(Hd(1,:),MRd,DTp);
% $$$ 
% $$$ % Check net time series:
% $$$ figure;
% $$$ inds = [500 250];
% $$$ for ii=1:length(inds)
% $$$     subplot(length(inds),1,ii);
% $$$     ind = inds(ii);
% $$$     plot(tvecSNAP,Hd(:,ind),'--k','linewidth',2);
% $$$     hold on;
% $$$     plot(tvecSNAP,Hd_fromN(:,ind),'--m','linewidth',2);
% $$$     plot(tvecSNAP,Hd_fromF(:,ind),'--b','linewidth',2);
% $$$     plot(tvecSNAP,Hd_fromM(:,ind),'--r','linewidth',2);
% $$$     plot(tvecSNAP,Hd_fromS(:,ind),'--g','linewidth',2);
% $$$     plot(tvecSNAP,Hd_fromF(:,ind)+Hd_fromM(:,ind)-Hd_fromM(1,ind),':c','linewidth',1);
% $$$     legend('$H(\mathcal{V},t_{snap})$','Integrated dH/dt','Integrated Surface Forcing','Integrated Mixing','Shortwave Redistribution');%,'Left-over linfit');
% $$$     title(['Heat Content above Vol Frac = ' num2str(ind/Vl)]);
% $$$     xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$     datetick(gca,'x','YYYY','keeplimits');
% $$$ end
% $$$ 
% $$$ 
% $$$ clim = [-2 2]*1e23;
% $$$ sp = 0.025e23;
% $$$ cpts = [clim(1):sp:clim(2)];
% $$$ 
% $$$ % Color plots:
% $$$ VF = 1;
% $$$ if (VF)
% $$$     [Xs,Ys] = ndgrid(tvecSNAP,Vi/Vi(end));
% $$$     lab = 'Volume Fraction';
% $$$     rev = 1;
% $$$ else
% $$$     [Xs,Ys] = ndgrid(tvecSNAP,mean(RUNS(1).T_ofV(:,tinds),2));
% $$$     lab = 'Temperature ($^\circ$C)';
% $$$     rev = 0;
% $$$ end
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ subplot(1,3,1);
% $$$ contourf(Xs,Ys,Hd,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$H(\mathcal{V},t_{snap}) [J]$');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,2);
% $$$ contourf(Xs,Ys,Hd_fromF,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\int F(\mathcal{V},t)dt$ [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,3);
% $$$ contourf(Xs,Ys,Hd_fromM,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\int M(\mathcal{V},t)dt$ [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ subplot(1,3,1);
% $$$ contourf(Xs,Ys,Hd_fromN,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\int \frac{\partial H}{\partial t}(\mathcal{V},t)dt$ [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,2);
% $$$ contourf(Xs,Ys,Hd-Hd_fromN,cpts,'linestyle','none');
% $$$ caxis(clim/10);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$H-\int \frac{\partial H}{\partial t}(\mathcal{V},t)dt$ [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,3);
% $$$ contourf(Xs,Ys,Hd_fromS,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'$\int SW(\mathcal{V},t)dt$ [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ 
% $$$ clim = [-0.5 0.5]*1e23;
% $$$ sp = 0.01e23;
% $$$ cpts = [-1e50 clim(1):sp:clim(2) 1e50];
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ subplot(1,3,1);
% $$$ contourf(Xs,Ys,Hd_fromMI,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'Numerical Mixing [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,2);
% $$$ contourf(Xs,Ys,Hd_fromMV,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'Vertical Mixing [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');
% $$$ subplot(1,3,3);
% $$$ contourf(Xs,Ys,Hd_fromMR,cpts,'linestyle','none');
% $$$ caxis(clim);
% $$$ colormap(redblue);
% $$$ cb = colorbar('Location','northoutside');
% $$$ ylabel(cb,'Redi Mixing [J]');
% $$$ xlabel('Year');
% $$$ ylabel(lab);
% $$$ if (rev); set(gca,'ydir','reverse'); end;
% $$$ xlim([tvecSNAP(1) tvecSNAP(end)]);
% $$$ datetick(gca,'x','YYYY','keeplimits');

% Linear trend line plots:
[ts,b] = linfit(tvecSNAP*86400,Hd/1e15);
Hdtr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromN/1e15);
Hd_fromNtr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromF/1e15);
Hd_fromFtr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromM/1e15);
Hd_fromMtr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromS/1e15);
Hd_fromStr = b(2,:)';

[ts,b] = linfit(tvecSNAP*86400,Hd_fromMI/1e15);
Hd_fromMItr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromMV/1e15);
Hd_fromMVtr = b(2,:)';
[ts,b] = linfit(tvecSNAP*86400,Hd_fromMR/1e15);
Hd_fromMRtr = b(2,:)';

figure;
x = Vi/Vi(end);
x = mean(RUNS(2).T_ofV(:,tinds),2);
plot(Hdtr,x,'-k','linewidth',2);
hold on;
plot(Hd_fromNtr,x,'--k','linewidth',2);
plot(Hd_fromFtr,x,'--b','linewidth',2);
plot(Hd_fromMtr,x,'--r','linewidth',2);
% $$$ plot(Hd_fromFtr+Hd_fromMtr,x,':k','linewidth',2);
plot(Hd_fromStr,x,'--g','linewidth',2);
% $$$ plot(Hd_fromFtr-Hd_fromStr,x,':b','linewidth',2);
plot(Hd_fromMItr,x,':b','linewidth',2);
plot(Hd_fromMVtr,x,':r','linewidth',2);
plot(Hd_fromMRtr,x,':','color',[0 0.5 0],'linewidth',2);
% $$$ set(gca,'ydir','reverse');
legend('$H(\mathcal{V},t_{snap})$','Integrated dH/dt', ...
       'Integrated Surface Forcing','Integrated Mixing', ...
% $$$        'Integrated Surface Forcing + Mixing', ...
       'Shortwave Redistribution','Numerical Mixing', ...
       'Vertical Mixing','Redi Mixing');%,'Left-over linfit');
xlabel('Heat Flux (PW) 2017-1972');
ylabel('Temperature ($^\circ$C)');
xlim([-0.08 0.2]);
ylim([-2 34]);
grid on;

% $$$ %%% Save variables for Taimoor:
% $$$ VTmaps.time = tvec;
% $$$ VTmaps.time_snap = tvecSNAP;
% $$$ VTmaps.vol = Vi;
% $$$ VTmaps.temp = mean(RUNS(1).T_ofV(:,tinds),2);
% $$$ VTmaps.z = z;
% $$$ VTmaps.thetaV = ThetaV;
% $$$ VTmaps.thetaZ = ThetaZ;
% $$$ VTmaps.heat_snap = Hd;
% $$$ VTmaps.heat_snap_fromSFR = Hd_fromF;
% $$$ VTmaps.heat_snap_fromMIX = Hd_fromM;
% $$$ VTmaps.heat_snap_fromDHDT = Hd_fromN;
% $$$ VTmaps.heat_snap_fromSWRD = Hd_fromS;
% $$$ VTmaps.heat_snap_fromNUM_MIX = Hd_fromMI;
% $$$ VTmaps.heat_snap_fromVER_MIX = Hd_fromMV;
% $$$ VTmaps.heat_snap_fromRED_MIX = Hd_fromMR;
% $$$ 
% $$$ VTlintr.vol = Vi;
% $$$ VTlintr.temp = mean(RUNS(1).T_ofV(:,tinds),2);
% $$$ VTlintr.heat_snap = Hdtr;
% $$$ VTlintr.heat_snap_fromSFR = Hd_fromFtr;
% $$$ VTlintr.heat_snap_fromMIX = Hd_fromMtr;
% $$$ VTlintr.heat_snap_fromDHDT = Hd_fromNtr;
% $$$ VTlintr.heat_snap_fromSWRD = Hd_fromStr;
% $$$ VTlintr.heat_snap_fromNUM_MIX = Hd_fromMItr;
% $$$ VTlintr.heat_snap_fromVER_MIX = Hd_fromMVtr;
% $$$ VTlintr.heat_snap_fromRED_MIX = Hd_fromMRtr;
% $$$ 
% $$$ save('ACCESS-OM2_1deg_rdf_linear-fit-perturbation_data.mat','VTmaps','VTlintr');


% $$$ %%% Temperature vs. time:
% $$$ fields = { ...
% $$$ % $$$           {N, 'Internal HC Tendency $\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
% $$$ % $$$           {dHdt, 'Total HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$ % $$$           {EHC, 'External HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$ % $$$           {F+PI, 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$ % $$$           {M, 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {I, 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$           {R, 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {M+I+R, 'Total Mixing $\mathcal{M}+\mathcal{I}+\mathcal{R}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {HI, '$\mathcal{H}_I(\Theta)$','m',2,'-'}, ...
% $$$ % $$$           {H, '$\mathcal{H}(\Theta)$','m',2,'-'}, ...
% $$$ % $$$           {HE, '$\mathcal{H}_E(\Theta)$','m',2,'-'}, ...
% $$$           {T_ofV, '$\Theta(\mathcal{V},t)$','m',2,'-'}, ...
% $$$           };
% $$$ 
% $$$ 
% $$$ % Scales and labels:
% $$$ 
% $$$ % Fluxes:
% $$$ scale = 1/1e15;label = '(PW)';x = Te;
% $$$ caxs = [-0.5 0.5];
% $$$ sp = 0.01;
% $$$ 
% $$$ % Temperature:
% $$$ scale = 1;label = '$^\circ$C';x = Te;
% $$$ caxs = [-0.15 0.15];
% $$$ sp = 0.005;
% $$$ 
% $$$ % Heat content:
% $$$ scale = 1/1e23;label = '$10^{23}J$';x = Te;
% $$$ caxs = [-2 2];
% $$$ sp = 0.01;
% $$$ 
% $$$ % Time-integrate fluxes:
% $$$ Tint = 0;
% $$$ 
% $$$ % remapping for percentiles:
% $$$ premap = 1;
% $$$ 
% $$$ % Subtract climatology:
% $$$ Sclim = 0;
% $$$ climean = [1972 1981];
% $$$ 
% $$$ % Annual average:
% $$$ AA = 0;
% $$$ 
% $$$ % $$$ % remap to depth:
% $$$ % $$$ remapz = 0;
% $$$ 
% $$$ cint = [-1e10 caxs(1):sp:caxs(2) 1e10];
% $$$ 
% $$$ for ii=1:length(fields)
% $$$     subplot(1,length(fields),ii);
% $$$     tvec = dnum;
% $$$     if (length(fields{ii}{1}(:,1)) == length(Te))
% $$$         x = Te;
% $$$         Torp = 1;
% $$$     elseif (length(fields{ii}{1}(:,1)) == TL)
% $$$         x = T;
% $$$         Torp = 1;
% $$$     elseif (length(fields{ii}{1}(:,1)) == Vl)
% $$$         if (premap)
% $$$             x = mean(T_ofp,2);
% $$$             Torp = 1;
% $$$         else
% $$$             x = Vi/Vtot;
% $$$             Torp = 0;
% $$$         end
% $$$     end
% $$$     % Annual average:
% $$$     if (AA)
% $$$         var = squeeze(monmean(fields{ii}{1},2,ndays))'*scale;
% $$$     else
% $$$         var = fields{ii}{1}*scale;
% $$$     end
% $$$     yrvec = unique(dvec(:,1));
% $$$     % Subtract climatology:
% $$$     if (Sclim)
% $$$         var = var-repmat(mean(var(:,find(yrvec>=climean(1) & ...
% $$$                                        yrvec<=climean(2)),:),2),[1 length(yrvec)]);
% $$$     end
% $$$     % Time-integrate:
% $$$     if (Tint == 1)
% $$$         var = cumsum(var*86400*DT,1);
% $$$     end
% $$$     [X,Y] = ndgrid(yrvec,x);
% $$$     contourf(X,Y,var',cint,'linestyle','none');
% $$$     
% $$$     cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$     ylabel(cb,label);
% $$$     if (Torp)
% $$$         ylim([-3 31]);
% $$$         set(gca,'ytick',-5:5:35);
% $$$     else
% $$$         ylim([0 1]);
% $$$         set(gca,'ydir','reverse');
% $$$     end        
% $$$     grid on;
% $$$     caxis(caxs);
% $$$     if (ii == 1)
% $$$         if (Torp)
% $$$             ylabel('Temperature ($^\circ$C)');
% $$$         else
% $$$             ylabel('Volume fraction $p$');
% $$$         end
% $$$     end
% $$$     xlabel('Date');
% $$$     title(fields{ii}{2});
% $$$     set(gca,'FontSize',15);
% $$$ end
% $$$ colormap(redblue);
% $$$ 
% $$$ 
% $$$ % $$$ %%% ARGO data:
% $$$ % $$$ 
% $$$ % $$$ load('ArgoData.mat');
% $$$ % $$$ yrsargo = 2004:2014;
% $$$ % $$$ Teargo = Tedge;
% $$$ % $$$ Hargo = zeros(length(Teargo),length(yrsargo));
% $$$ % $$$ Vargo = Hargo;
% $$$ % $$$ dvec = datevec(days);
% $$$ % $$$ for yi=1:length(yrsargo)
% $$$ % $$$     yr = yrsargo(yi);
% $$$ % $$$     Hargo(:,yi) = mean(Hsnap(:,dvec(:,1)==yr),2);
% $$$ % $$$     Vargo(:,yi) = mean(Vsnap(:,dvec(:,1)==yr),2);
% $$$ % $$$ end
% $$$ % $$$ HIargo = Hargo-rho0*Cp*Vargo.*repmat(Teargo',[1 length(yrsargo)]);
% $$$ % $$$ HEargo = Hargo-HIargo;
% $$$ % $$$ 
% $$$ % $$$ [X,Y] = ndgrid(yrsargo,Teargo);
% $$$ % $$$ 
% $$$ % $$$ Hargo = Hargo-repmat(mean(Hargo,2),[1 length(yrsargo)]);
% $$$ % $$$ HEargo = HEargo-repmat(mean(HEargo,2),[1 length(yrsargo)]);
% $$$ % $$$ HIargo = HIargo-repmat(mean(HIargo,2),[1 length(yrsargo)]);
% $$$ % $$$ 
% $$$ % $$$ % Fluxes:
% $$$ % $$$ scale = 1/1e23;label = '(10$^{23}$J)';
% $$$ % $$$ caxs = [-0.5 0.5];
% $$$ % $$$ sp = 0.025;
% $$$ % $$$ 
% $$$ % $$$ cint = [-1e10 caxs(1):sp:caxs(2) 1e10];
% $$$ % $$$ 
% $$$ % $$$ fields = { ...
% $$$ % $$$           {[], 'Internal Heat Content $\mathcal{H}_I$','m',2,'-'}, ...
% $$$ % $$$           {[], 'Heat Content $\mathcal{H}$','m',2,'-'}, ...
% $$$ % $$$           {[], 'External Heat Content $\mathcal{H}_E$','m',2,'-'}, ...
% $$$ % $$$           };
% $$$ % $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ subplot(1,3,1);
% $$$ % $$$ contourf(X,Y,HIargo'*scale,cint,'linestyle','none');
% $$$ % $$$ ylabel('Temperature ($^\circ$C)');
% $$$ % $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ % $$$ set(gca,'ytick',-5:5:35);
% $$$ % $$$ ylabel(cb,label);
% $$$ % $$$ xlabel('Date');
% $$$ % $$$ ylim([-3 31]);
% $$$ % $$$ grid on;
% $$$ % $$$ caxis(caxs);
% $$$ % $$$ xlim([2004 2014]);
% $$$ % $$$ title(fields{1}{2});
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ % $$$ subplot(1,3,2);
% $$$ % $$$ contourf(X,Y,Hargo'*scale,cint,'linestyle','none');
% $$$ % $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ % $$$ set(gca,'ytick',-5:5:35);
% $$$ % $$$ ylabel(cb,label);
% $$$ % $$$ xlabel('Date');
% $$$ % $$$ ylim([-3 31]);
% $$$ % $$$ grid on;
% $$$ % $$$ caxis(caxs);
% $$$ % $$$ xlim([2004 2014]);
% $$$ % $$$ title(fields{2}{2});
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ % $$$ subplot(1,3,3);
% $$$ % $$$ contourf(X,Y,HEargo'*scale,cint,'linestyle','none');
% $$$ % $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ % $$$ set(gca,'ytick',-5:5:35);
% $$$ % $$$ ylabel(cb,label);
% $$$ % $$$ xlabel('Date');
% $$$ % $$$ ylim([-3 31]);
% $$$ % $$$ grid on;
% $$$ % $$$ caxis(caxs);
% $$$ % $$$ xlim([2004 2014]);
% $$$ % $$$ title(fields{3}{2});
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ % $$$ colormap(redblue);
% $$$ 
% $$$ %%%% Depth space plot:
% $$$ 
% $$$ % This script makes plots of the heat budget in the MOM
% $$$ % simulations.
% $$$ 
% $$$ close all;
% $$$ clear all;
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ 
% $$$ RUNS = { ...
% $$$          {'ACCESS-OM2_1deg_jra55_rdf',[51:55],[1972]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_rdf_pert',[51:55],[1972]}, ...
% $$$        };
% $$$ 
% $$$ rr = 2;
% $$$     rr
% $$$     outputs = RUNS{rr}{2};
% $$$     model = RUNS{rr}{1};
% $$$ 
% $$$ % $$$     clearvars -except base RUNS rr outputs model leg legh;
% $$$     
% $$$     load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$     if (~exist('ndays'))
% $$$         ndays = diff(time_snap);
% $$$     end
% $$$     region = 'Global';
% $$$     nyrs = tL/12;
% $$$     if (nyrs == round(nyrs))
% $$$         szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
% $$$         yrs = 1:nyrs;
% $$$     else
% $$$         nyrs = 1;
% $$$         szTe = [TL+1 tL];szT = [TL tL];
% $$$     end    
% $$$     ycur = 1;
% $$$     
% $$$     dnum = [];
% $$$     
% $$$     Vs = zeros(zL,nyrs,1);
% $$$     Hs = zeros(zL,nyrs,1);
% $$$ 
% $$$     %% Load Global Budget:
% $$$     for i=1:length(outputs)
% $$$         load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
% $$$         dnum = [dnum; time];
% $$$ 
% $$$         load([base model sprintf('_output%03d_',outputs(i)) 'VHofz.mat']);
% $$$         Vs(:,:,i) = V;
% $$$         Hs(:,:,i) = H;
% $$$         ycur = ycur+nyrs;
% $$$     end
% $$$     months = [1:length(Vs(1,:,1))];
% $$$     yrs = [1:length(Vs(1,1,:))];
% $$$     
% $$$     % Correct time vector zero year:
% $$$     dvec = datevec(dnum);
% $$$     if (length(RUNS{rr}{3})== 1)
% $$$         dvec(:,1) = dvec(:,1)-dvec(1,1) + RUNS{rr}{3}(1);
% $$$     else
% $$$         dvec(:,1) = dvec(:,1) + 1900;
% $$$     end
% $$$     dnum = datenum(dvec);
% $$$     yrs = 1:10;
% $$$     
% $$$     Ts = Hs/rho0/Cp./Vs;
% $$$     Vs = cumsum(Vs,1);
% $$$     Hs = cumsum(Hs,1);
% $$$ 
% $$$ % Subtract climatology:
% $$$ Sclim = 0;
% $$$ climean = [1972 1981];
% $$$ 
% $$$ % Annual average:
% $$$ AA = 0;
% $$$ 
% $$$ % Annual average:
% $$$ if (AA)
% $$$     T =  squeeze(monmean(Ts,2,ndays));
% $$$     H =  squeeze(monmean(Hs,2,ndays));
% $$$     V =  squeeze(monmean(Vs,2,ndays));
% $$$     tvec = unique(dvec(:,1));
% $$$ else
% $$$     V = reshape(Vs,[length(Vs(:,1,1))  nyrs*length(Vs(1,1,:))]);
% $$$     H = reshape(Hs,[length(Vs(:,1,1)) nyrs*length(Hs(1,1,:))]);
% $$$     T = reshape(Ts,[length(Vs(:,1,1)) nyrs*length(Hs(1,1,:))]);
% $$$ end
% $$$ 
% $$$ % $$$ Tcont = T;
% $$$ T = T - Tcont;
% $$$ 
% $$$ % Subtract climateology:
% $$$ yrvec = unique(dvec(:,1));
% $$$ % Subtract climatology:
% $$$ if (Sclim)
% $$$     T = T-repmat(mean(T(:,find(yrvec>=climean(1) & ...
% $$$                                      yrvec<=climean(2)),:),2),[1 length(yrvec)]);
% $$$     V = V-repmat(mean(V(:,find(yrvec>=climean(1) & ...
% $$$                                      yrvec<=climean(2)),:),2),[1 length(yrvec)]);
% $$$     H = H-repmat(mean(H(:,find(yrvec>=climean(1) & ...
% $$$                                      yrvec<=climean(2)),:),2),[1 length(yrvec)]);
% $$$ end
% $$$ [X,Y] = ndgrid(yrvec,-z);
% $$$ % $$$ subplot(1,2,1);
% $$$ % $$$ pcolPlot(X,Y,H');
% $$$ % $$$ caxs = [-0.8 0.8]*1e23;
% $$$ % $$$ sp = 0.005e23
% $$$ % $$$ cint = [-1e50 caxs(1):sp:caxs(2) 1e50];
% $$$ % $$$ contourf(X,Y,H',cint,'linestyle','none');
% $$$ % $$$ title('Heat Content Anomaly above depth level (J)');
% $$$ caxs = [-0.15 0.15];
% $$$ sp = 0.005
% $$$ cint = [-1e50 caxs(1):sp:caxs(2) 1e50];
% $$$ contourf(X,Y,T',cint,'linestyle','none');
% $$$ cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$ ylabel('Depth (m)');
% $$$ xlabel('Date');
% $$$ title('Temperature Anomaly $(^\circ$C)');%Heat Content Anomaly above depth level (J)');
% $$$ caxis(caxs);%[-0.5 0.5]*1e23);
% $$$ set(gca,'FontSize',15);
% $$$ colormap(redblue);
% $$$ xlim([1972 2018]);
