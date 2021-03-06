% This script makes plots of the heat budget in the MOM
% simulations.

% $$$ close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';
% $$$ base = 'C://Users//holme//data//mom//';

RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[4567]}, ...
% $$$ % MOM025-SIS:
%    {'MOM025_kb3seg',[101120]}, ...
% $$$     {'MOM025_kb3seg',[120]}, ...
% $$$     {'MOM025',[15:19]}, ...
% $$$     {'MOM025_kb1em5',[95:99]}, ...
% $$$     {'MOM025_kb1em6',[30]}, ...
% $$$ % ACCESS-OM2 Gadi runs:
% $$$ % $$$ % 1-degree
% $$$          {'ACCESS-OM2_1deg_jra55_ryf',[31]}, ...
% $$$ % $$$          {'ACCESS-OM2_1deg_jra55_ryf_sgl',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135]}, ...
% $$$ % 1/4-degree
% $$$          {'ACCESS-OM2_025deg_jra55_ryf',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[82]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM_dt2',[81]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_norediGM_smoothkppbl',[81]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5',[7781]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kbvar',[7781]}, ...
% $$$          {'ACCESS-OM2_025deg_jra55_ryf_kds75',[7680]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf8485_gmredi',[73]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[80]}, ...
% $$$ % $$$          {'ACCESS-OM2_025deg_jra55_ryf',[300]}, ...
% $$$          {'ACCESS-OM2_01deg_jra55_ryf',[636643]}, ...
% $$$ % $$$          {'ACCESS-OM2_1deg_jra55_rdf',[51:55]}, ...
% $$$ % $$$          {'ACCESS-OM2_1deg_jra55_rdf_pert',[51:55]}, ...
% $$$          {'ACCESS-OM2_01deg_jra55_ryf_k_smag_iso3',[640643]}, ...
         {'ACCESS-OM2_01deg_jra55_ryf_k_smag_iso7p9',[648:655]}, ...
       };
Inetstr = [];

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

% $$$ rr = 1;
for rr = 1:length(RUNS);
    rr
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model leg legh Inetstr;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    region = 'Global';
% $$$     region = 'IndoPacific';
    if (mod(tL,12) == 0) % monthly output
        nyrs = tL/12
        szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
        yrs = 1:nyrs;
    else
        nyrs = tL;
        ndays = ndays./ndays;
        szTe = [TL+1 1 tL];szT = [TL 1 tL];
    end    
    ycur = 1;

    %% Global Calculations:
    for i=1:length(outputs)
        
% $$$     % Annual or Monthly offline Binning:
% $$$     load([base model sprintf('_output%03d_',outputs(i)) 'GlobalHBud_MonAnBin.mat']);
% $$$     GWB = GWBann;
        try
            load([base model sprintf('_output%03d_',outputs(i)) region '_HBud.mat']);
        catch        
            load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
        end
        
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
        if (isfield(GWB,'NUMDISS')) % Burchard numerical mixing
            NUMDISS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUMDISS,szT); % NUMDISS (W)
        else
            NUMDISS(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
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
        JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 szTe(2) szTe(3)])*rho0*Cp;

        % Interior heat source P:
        PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));

        % Interior heat source Q:
        QII(:,:,ycur:(ycur+nyrs-1)) = QI(:,:,ycur:(ycur+nyrs-1)) - JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 szTe(2) szTe(3)])*rho0*Cp;

        % Across-isotherm advective heat flux:
        CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 szTe(2) szTe(3)])*rho0*Cp;

        % External HC Tendency:
        EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 szTe(2) szTe(3)])*rho0*Cp;

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
        I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) ...
            - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1));% ...
% $$$             - SUB(:,:,ycur:(ycur+nyrs-1)) - GM(:,:,ycur:(ycur+nyrs-1));

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
        HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 szTe(2) szTe(3)]);
        HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 szTe(2) szTe(3)]);
        HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 szTe(2) szTe(3)]);
        HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 szTe(2) szTe(3)]);
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
    
    % Print some overall numbers on numerical mixing:
    Inet = sum(I*dT,1);
    Mnet = sum(M*dT,1)+sum(R*dT,1);
% $$$     [tmp ind] = min(abs(Te-5));
% $$$     Inetstr = [Inetstr model sprintf(' Avg. I(5C) = %5.2f PW',mean(monmean(I(ind,:,:),2,ndays(1:length(Inet(1,:,1)))),3)/1e15) ' \n '];
% $$$     [tmp ind] = min(abs(Te-15));
% $$$     Inetstr = [Inetstr model sprintf(' Avg. I(15C) = %5.2f PW',mean(monmean(I(ind,:,:),2,ndays(1:length(Inet(1,:,1)))),3)/1e15) ' \n '];
% $$$     [tmp ind] = min(abs(Te-22.5));
% $$$     Inetstr = [Inetstr model sprintf(' Avg. I(22.5C) = %5.2f PW',mean(monmean(I(ind,:,:),2,ndays(1:length(Inet(1,:,1)))),3)/1e15) ' \n '];
    Inetstr = [Inetstr sprintf(' Inet = %5.1f PWdegC, Mnet = %5.1f PWdegC',mean(monmean(Inet,2,ndays(1:length(Inet(1,:,1)))),3)/1e15,mean(monmean(Mnet,2,ndays(1:length(Mnet(1,:,1)))),3)/1e15) ' ' model ' \n '];

% $$$ end
% $$$ fprintf(Inetstr)


% $$$ 
% $$$     yrs = [1 5];
% $$$
% $$$
    %%%%Heat Flux: ---------------------------------------------------------------------------------------------
% $$$ % Production fields:
    fields = { ...
        {N(:,months,yrs), 'Tendency $\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
        {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
        {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
        {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
        {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
        {NUMDISS(:,months,yrs), 'Numerical Mixing $\mathcal{I}$ Burchard','b',2,'--'}, ...
% $$$           {NUM(:,months,yrs), '3D Numerical Mixing Global Sum','g',2,':'}, ...
% $$$           {GM(:,months,yrs)+SUB(:,months,yrs), 'Submesoscale','c',2,'--'}, ...
% $$$           {TEN(:,months,yrs), 'Eulerian Tendency','m',1,'--'}, ...
% $$$           {dHdt(:,months,yrs), 'dH/dt','g',1,'--'}, ...
% $$$           {ADV(:,months,yrs), 'Eulerian Advection','b',1,'--'}
             };

% $$$     % Eulerian checks:
% $$$     fields = { ...
% $$$         {N(:,months,yrs), '$\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
% $$$         {F(:,months,yrs)+M(:,months,yrs)+R(:,months,yrs), 'Forcing + Explicit Mixing $ECR=\mathcal{M}+\mathcal{F}+\mathcal{R}$','k',2,'-'}, ...
% $$$         {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$         {PI(:,months,yrs), 'Surface Volume Internal $\mathcal{P}_I$','r',2,'-'}, ...
% $$$         {P(:,months,yrs), 'Surface Volume $\mathcal{P}$','r',2,'--'}, ...
% $$$         {TEN(:,months,yrs), 'Eulerian TEN','m',2,'--'}, ...
% $$$         {dHdt(:,months,yrs), 'dH/dt','m',2,':'}, ...
% $$$         {ADV(:,months,yrs), 'Eulerian ADV','b',2,'--'}, ...
% $$$         {TEN(:,months,yrs)-ADV(:,months,yrs)+JSH(:,months,yrs), 'TEN-ADV-JSH','m',1,'--'}, ...
% $$$              };

    % Save global budget:
% $$$     GloDiaHB.temp = Te;
% $$$     GloDiaHB.surface_forcing = mean(monmean(F(:,months,yrs)+PI(:,months,yrs),2,ndays(months)),3)/1e15;
% $$$     GloDiaHB.tendency = mean(monmean(N(:,months,yrs),2,ndays(months)),3)/1e15;
% $$$     GloDiaHB.vertical_mixing = mean(monmean(M(:,months,yrs),2,ndays(months)),3)/1e15;
% $$$     GloDiaHB.numerical_mixing = mean(monmean(R(:,months,yrs),2,ndays(months)),3)/1e15;
% $$$     GloDiaHB.redi_mixing = mean(monmean(I(:,months,yrs),2,ndays(months)),3)/1e15;

    fields = { ...
        {F(:,months,yrs)+PI(:,months,yrs)-N(:,months,yrs), 'Forcing - Tendency (PW)','k',2,'-'}, ...
        {R(:,months,yrs)+M(:,months,yrs)+MD(:,months,yrs)+SG(:,months,yrs), 'Explicit Mixing (PW)','r',2,':'}, ...
% $$$         {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,':'}, ...
        {I(:,months,yrs), 'Numerical Mixing (PW)','b',2,'--'}, ...
% $$$           {NUM(:,months,yrs), '3D Numerical Mixing Global Sum','g',2,':'}, ...
% $$$           {GM(:,months,yrs)+SUB(:,months,yrs), 'Submesoscale','c',2,'--'}, ...
% $$$           {TEN(:,months,yrs), 'Eulerian Tendency','m',1,'--'}, ...
% $$$           {dHdt(:,months,yrs), 'dH/dt','g',1,'--'}, ...
% $$$           {ADV(:,months,yrs), 'Eulerian Advection','b',1,'--'}
             };

% $$$ 
    Fscale = 1/1e15;
% $$$ 
% $$$ % $$$ yrtyps = {'-','--','-.',':','-d','-s','--d','--s',':d',':s','-.d','-.s','-o'}; % line-types for different years
% $$$ typs = {'-','-','--',':','-.','-','--',':','-','-','--',':','-.'}; % line-types for different years
% $$$ cols = {'m','k','k','k','k','r','r','r',[0.302 0.7451 0.9333],'b','b','b','b'};
% $$$ % $$$ wids = {2,2,2,2,2,2,2,2,2,2,2,2,2,2};
% $$$ typs = {'-','-','--',':','-.','-','--',':','-','-','--',':','-.'}; % line-types for different years
% $$$ wids = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    typs = {'-','--',':','-.','-','-',':','-.','--',':'};
    wids = {2,2,2,2,2,2,1,1,1,1,1,1,1,1};
    cols = {'b','r','k','m','g'};
    %Fluxes only:
% $$$ figure;
% $$$ set(gcf,'Position',[207          97        1609         815]);
% $$$ leg = {};
% $$$ legh = [];
    for i=1:length(fields)
        hold on;
        if (length(fields{i}{1}(:,1)) == length(Te))
            x = Te;
        else
            x = T;
        end
        
        % Plot years from a single run separately:
% $$$     for j=1:length(yrs) 
% $$$         h = plot(Te,monmean(fields{i}{1}(:,:,yrs(j)),2,ndays(months))*Fscale,yrtyps{j}, 'color',fields{i}{3} ...
% $$$              ,'linewidth',3);
% $$$         if (j == 1)
% $$$             legh(i) = h;
% $$$         end
% $$$     end
% $$$     leg{i} = fields{i}{2};

% $$$     % Average years together for a single run:
% $$$     legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$     leg{i} = fields{i}{2};
        
        % Average years together for multiple runs:
        tmp = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),typs{rr}, 'color',fields{i}{3} ...
                   ,'linewidth',wids{rr});
% $$$         tmp = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',cols{rr} ...
% $$$                    ,'linewidth',wids{rr});
% $$$         tmp = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),typs{rr}, 'color',cols{rr}, 'linewidth',wids{rr});
% $$$         if i==1
% $$$             leg{rr} = strrep(RUNS{rr}{1},'_',' ');
% $$$             legh(rr) = tmp;
% $$$         end
    end
    ylim([-1.5 1.5]);
    xlim([-3 31]);
    box on; 
    grid on;
    ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
    xlabel('Temperature $\Theta$ ($^\circ$C)');
% $$$ lg = legend(legh,leg);
% $$$ set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);

end

% $$$ % Print table at specific values:
% $$$ [tmp ind1] = min(abs(Te-22.5));
% $$$ [tmp ind2] = min(abs(Te-5));
% $$$ 
% $$$ for i=1:length(fields)
% $$$     var = mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3);
% $$$     RUNn = strrep(RUNS{rr}{1},'_',' ');
% $$$     varn = fields{i}{2};
% $$$     [RUNn '$5^\circ$C: ' sprintf('$%2.2f$PW',var(ind2)) ', $22.5^\circ$C: ' sprintf('$%2.2f$PW',var(ind1))]
% $$$ end

% $$$ end

% $$$ end

% $$$ % Plotting tendency terms for Jonathan:
% $$$ figure;
% $$$ plot(Te,dHdt(:,1,1)/1e15,'-k');
% $$$ hold on;
% $$$ plot(Te(1:4:end),dHdt(1:4:end,1,1)/1e15,'--k');
% $$$ plot(Te(1:16:end),dHdt(1:16:end,1,1)/1e15,'-xk','linewidth',3);
% $$$ legend('0.5$^\circ$C T-bins','2$^\circ$C T-bins','8$^\circ$C T-bins');
% $$$ plot(Te,N(:,1,1)/1e15,'-m');
% $$$ hold on;
% $$$ plot(Te(1:4:end),N(1:4:end,1,1)/1e15,'--m');
% $$$ plot(Te(1:16:end),N(1:16:end,1,1)/1e15,'-xm','linewidth',3);
% $$$ plot(Te,(dHdt(:,1,1)-N(:,1,1))/1e15,'-b');
% $$$ hold on;
% $$$ plot(Te(1:4:end),(dHdt(1:4:end,1,1)-N(1:4:end,1,1))/1e15,'--b');
% $$$ plot(Te(1:16:end),(dHdt(1:16:end,1,1)-N(1:16:end,1,1))/1e15,'-xb','linewidth',3);
% $$$ text(5,-4,'Total Heat Content Tendency d$\mathcal{H}$/dt');
% $$$ text(5,-6,'Internal Heat Content Tendency d$\mathcal{H}_I$/dt','color','m');
% $$$ text(5,-8,'External Heat Content Tendency d$\mathcal{H}_E$/dt','color','b');
% $$$ xlabel('Temperature ($^\circ$C)');
% $$$ ylabel('Net Heat Flux Tendency (PW)');
% $$$ title('January MOM025 Diathermal Heat Budget');

% $$$ 
% $$$ % Region averages for different runs:
% $$$ ff = mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3);
% $$$ [tmp ind1] = min(abs(Te-17.5));
% $$$ [tmp ind2] = min(abs(Te-25));
% $$$ RUNS{rr}{3} = mean(ff(ind1:ind2));
% $$$ ff = mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3);
% $$$ [tmp ind1] = min(abs(Te-0));
% $$$ [tmp ind2] = min(abs(Te-10));
% $$$ RUNS{rr}{4} = mean(ff(ind1:ind2));
% $$$ 
% $$$ end
% $$$ 
% $$$ % Plot bar graph:
% $$$ data = zeros(length(RUNS),2);
% $$$ labs = cell(length(RUNS),1);
% $$$ for rr=1:length(RUNS)
% $$$     data(rr,1) = RUNS{rr}{3};
% $$$     data(rr,2) = RUNS{rr}{4};
% $$$     labs{rr} = strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55','');
% $$$ end
% $$$ 
% $$$ % $$$ order = 1:length(RUNS);
% $$$ % $$$ order = [6 3 4 1 2 5 12 14 13 9 10 11];
% $$$ redi = [0 0 0 0 0 1 0 1 1 1 1 1];
% $$$ bakd = [0 0 1 1 1 0 0 0 0 0 0 0];
% $$$ horbar = [1 9];
% $$$ clf;
% $$$ b = barh(data(order,:),'BarWidth',1.5)
% $$$ xlabel('Numerical Mixing (PW)');
% $$$ cnt = 1;
% $$$ for rr=order
% $$$     text(-0.98,cnt,labs{rr});
% $$$     if (redi(cnt))
% $$$         text(-0.65,cnt,'REDI','color','r');
% $$$     end
% $$$     if (bakd(cnt))
% $$$         text(-0.65,cnt,'BAKD','color','c');
% $$$     end
% $$$     cnt = cnt+1;
% $$$ end
% $$$ hold on;
% $$$ for ii=1:length(horbar)
% $$$     plot([-1 0],(horbar(ii)+0.5)*[1 1],'--k','linewidth',2);
% $$$ end
% $$$ set(gca,'ytick',[]);
% $$$ title('Blue = Numerical Mixing $17.5^\circ$C-$25^\circ$C, Yellow = Numerical Mixing $0^\circ$C-$10^\circ$C');
% $$$ 
% $$$ %%%%WM Transformation / Volume Budget:
% $$$ % 05-12-17 Note: The volume budget of the Pacific now closes
% $$$ % satisfactorily, after fixing masks and including submeso transport.
% $$$ % I.e. the WMT term G (calculated via residual) is now 0.0046 Sv at
% $$$ % -3C, whereas before it was 0.5Sv (it should be zero).
% $$$ months = [1:12];
% $$$ 
% $$$ fields = { ...
% $$$           {dVdt(:,months,yrs), 'Tendency $\frac{\partial\mathcal{V}}{\partial t}$','m',2,'-'}, ...
% $$$           {JS(:,months,yrs), 'Surface Volume Flux $\mathcal{J}_S$',0.5*[1 1 1],2,':'}, ...
% $$$ % $$$           {G(:,months,yrs), 'WMT $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$           {WMT(:,months,yrs), 'Total WMT $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$           {WMTM(:,months,yrs), 'WMT $\mathcal{G}$ from Vertical Mixing','r',2,'-'}, ...
% $$$           {WMTF(:,months,yrs), 'WMT $\mathcal{G}$ from Surface Forcing','k',2,'-'}, ...
% $$$           {WMTI(:,months,yrs), 'WMT $\mathcal{G}$ from Implicit Mixing','b',2,'-'}, ...
% $$$ % $$$           {WMTR(:,months,yrs), 'Interior WMT $\mathcal{G}$ from Heat Fluxes','b',2,'--'}, ...
% $$$ % $$$           {-JI(:,months,yrs), 'ITF + SF + BS $-\mathcal{J_I}$',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {-JITF(:,months,yrs), 'ITF Volume Loss',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {-JSP(:,months,yrs), 'South Pacific Volume Loss',[0 0.5 0],2,'-.'}, ...
% $$$ % $$$           {-JBS(:,months,yrs), 'Bering Strait Volume Loss',[0 0.5 0],2,':'}, ...
% $$$           };
% $$$ 
% $$$ Mscale = 1/1e6;
% $$$ 
% $$$ %Fluxes only:
% $$$ figure;
% $$$ set(gcf,'Position',[207          97        1609         815]);
% $$$ leg = {};
% $$$ legh = [];
% $$$ for i=1:length(fields)
% $$$     hold on;
% $$$     if (length(fields{i}{1}(:,1)) == length(Te))
% $$$         x = Te;
% $$$     else
% $$$         x = T;
% $$$     end
% $$$     legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Mscale,3),fields{i}{5}, 'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$     leg{i} = fields{i}{2};
% $$$ end
% $$$ ylim([-50 80]);
% $$$ xlim([-3 31]);
% $$$ box on;
% $$$ grid on;
% $$$ ylabel('Water Mass Transformation (Sv)');
% $$$ xlabel('Temperature $\Theta$ ($^\circ$C)');
% $$$ lg = legend(legh,leg);
% $$$ set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);
% $$$ 
% $$$ %%% Temperature vs. time:
% $$$ months = [1:12];
% $$$ fields = { ...
% $$$           {N(:,months,yrs), 'Internal HC Tendency $\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
% $$$ % $$$           {dHdt(:,months,yrs), 'Total HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$ % $$$           {EHC(:,months,yrs), 'External HC Tendency $\partial\mathcal{H}/\partial t$','m',2,'-'}, ...
% $$$           {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$ % $$$           {F(:,months,yrs), 'Surface Heat Fluxes $\mathcal{F}$','k',2,'-'}, ...
% $$$ % $$$           {P(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}$',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$           {PI(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}_I$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$           {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {GM(:,months,yrs), 'GM $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {SUB(:,months,yrs), 'SUB $\mathcal{S}$',[0 0.5 0],2,':'}, ...
% $$$ % $$$           {HWMTI(:,months,yrs), 'Advective Implicit Mixing','b',2,'--'}, ...
% $$$ % $$$           {HWMTM(:,months,yrs), 'Advective Vertical Mixing','r',2,'--'}, ...
% $$$ % $$$           {HWMTF(:,months,yrs), 'Advective Surface Forcing','k',2,'--'}, ...
% $$$           {M(:,months,yrs)+I(:,months,yrs)+R(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}+\mathcal{R}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {Nmon(:,months,yrs), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$ % $$$           {SW(:,months,yrs), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$ % $$$           {dHdt(:,months,yrs), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',2,'--'}, ...
% $$$ % $$$           {CIA(:,months,yrs), 'Across-Isotherm Advection $\mathcal{G}\Theta\rho_0C_p$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {VDFkppiw(:,months,yrs), 'Background','b',2,'-'}, ...
% $$$ % $$$           {VDFkppish(:,months,yrs), 'Shear Instability','g',2,'-'}, ...
% $$$ % $$$           {M(:,months,yrs) - VDFsum(:,months,yrs) + VDFkppbl(:,months,yrs), 'Vertical Mixing SUM','r',2,'--'}, ...
% $$$ % $$$           {VDFkppbl(:,months,yrs), 'KPP Boundary Layer','c',2,'-'}, ...
% $$$ % $$$           {VDFwave(:,months,yrs), 'Vertical Diffusion WAVE','k',2,'-'}, ...
% $$$ % $$$           {VDFkppicon(:,months,yrs), 'Vertical Diffusion KPPICON','y',2,'-'}, ...
% $$$ % $$$           {VDFkppdd(:,months,yrs), 'Vertical Diffusion KPPDD','m',2,'-'}, ...
% $$$ % $$$           {VDFsum(:,months,yrs), 'Vertical Diffusion SUM','r',2,'--'}, ...
% $$$ % $$$           {KPPnloc(:,months,yrs), 'KPP Non-local','m',2,'--'}, ...
% $$$           };
% $$$ 
% $$$ % Fluxes:
% $$$ Tint = 0;
% $$$ scale = 1/1e15;label = '(PW)';x = Te;
% $$$ caxs = [-0.5 0.5];
% $$$ sp = 0.01;
% $$$ 
% $$$ Tint = 1;
% $$$ scale = 1/1e23;label = '($10^{23}$J)';x = Te;
% $$$ caxs = [-1 1];
% $$$ sp = 0.01;
% $$$ 
% $$$ % $$$ % Transformations:
% $$$ % $$$ scale = 1/1e6;label = '(Sv)';
% $$$ % $$$ caxs = [-250 250];
% $$$ % $$$ sp = 25;
% $$$ % $$$ caxs = [-70 70];
% $$$ % $$$ sp = 3.5;
% $$$ 
% $$$ cint = [-1e10 caxs(1):sp:caxs(2) 1e10];
% $$$ 
% $$$ % $$$ figure;
% $$$ %set(gcf,'Position',get(0,'ScreenSize'));
% $$$ % $$$ set(gcf,'Position',[3    40   956   963]);
% $$$ climean = [1:10];
% $$$ for ii=1:length(fields)
% $$$     subplot(1,length(fields),ii);
% $$$ % $$$     subplot_tight(3,4,rr);%1,length(fields),ii);
% $$$ % $$$     V = mean(fields{ii}{1},3)'*scale;
% $$$     V = squeeze(monmean(fields{ii}{1},2,ndays))'*scale;
% $$$     tL = length(V(:,1));
% $$$     V = V-repmat(mean(V(climean,:),1),[tL 1]);
% $$$     if (Tint == 1)
% $$$         V = cumsum(V*86400*365,1);
% $$$     end
% $$$     if (length(fields{ii}{1}(:,1)) == length(Te))
% $$$         x = Te;
% $$$     else
% $$$         x = T;
% $$$     end
% $$$ % $$$     [X,Y] = ndgrid(1:tL,x);
% $$$     [X,Y] = ndgrid(1976:2015,x);
% $$$     contourf(X,Y,V,cint,'linestyle','none');
% $$$     cb = colorbar('Location','NorthOutside','FontSize',15);    
% $$$     set(gca,'ytick',-5:5:35);
% $$$     ylabel(cb,label);
% $$$ % $$$     set(gca,'xtick',[1:tL]);
% $$$ % $$$     if (rr>=9)
% $$$ % $$$         set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun', ...
% $$$ % $$$                             'Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ % $$$     else
% $$$ % $$$         set(gca,'xticklabel',[]);
% $$$ % $$$     end
% $$$     ylim([-3 31]);
% $$$     grid on;
% $$$     caxis(caxs);
% $$$ % $$$     xlabel('Month');
% $$$ % $$$     ylabel('Temperature ($^\circ$C)');
% $$$ % $$$     xlabel(cb,strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
% $$$ % $$$     xlabel(cb,[model ' ' fields{ii}{2} ' ' ...
% $$$ % $$$                label],'FontSize',20);
% $$$ % $$$     title(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
% $$$     title(fields{ii}{2});%strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
% $$$     if (mod(rr,4) == 0)
% $$$         pos = get(gca,'Position');
% $$$         cb = colorbar;
% $$$         set(gca,'Position',pos);
% $$$     end
% $$$     set(gca,'FontSize',15);
% $$$ end
% $$$ cmap = redblue((length(cint)-3)*2);
% $$$ cmap = cmap(1:(length(cint)-3),:);
% $$$ colormap(cmap);
% $$$ colormap(redblue);
% $$$ % $$$ 
% $$$ % $$$ %% Global Seasonal Cycle TS
% $$$ % $$$ months = 1:12;
% $$$ % $$$ Ts = 21.5;
% $$$ % $$$ [tmp ind] = min(abs(Te-Ts));
% $$$ % $$$ 
% $$$ % $$$ Fscale = 1/1e15;
% $$$ % $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ %set(gcf,'Position',get(0,'ScreenSize'));
% $$$ % $$$ set(gcf,'Position',[34          40        1164         963]);
% $$$ % $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ subplot(2,1,1);
% $$$ % $$$ fields = {
% $$$ % $$$           {F(ind,months,yrs)+PI(ind,months,yrs), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$ % $$$           {N(ind,months,yrs), 'Total $\mathcal{N}$','m',2,'-'}, ...
% $$$ % $$$           {M(ind,months,yrs)+I(ind,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           };
% $$$ % $$$ 
% $$$ % $$$ leg = {};
% $$$ % $$$ legh = [];
% $$$ % $$$ for i=1:length(fields)
% $$$ % $$$     hold on;
% $$$ % $$$ % $$$     for j=1:length(P(1,1,:))
% $$$ % $$$ % $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$ % $$$ % $$$              'color',0.7*[1 1 1] ...
% $$$ % $$$ % $$$              ,'linewidth',0.5);
% $$$ % $$$ % $$$     end
% $$$ % $$$     legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
% $$$ % $$$          'color',fields{i}{3} ...
% $$$ % $$$          ,'linewidth',fields{i}{4});
% $$$ % $$$     hold on;
% $$$ % $$$     leg{i} = fields{i}{2};
% $$$ % $$$ end
% $$$ % $$$ xlabel('Month');
% $$$ % $$$ ylabel('PW');
% $$$ % $$$ lg = legend(legh,leg);
% $$$ % $$$ ylim([-5 5]);
% $$$ % $$$ xlim([1 tL]);
% $$$ % $$$ set(gca,'xtick',[1:tL]);
% $$$ % $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ % $$$ set(gca,'FontSize',25);
% $$$ % $$$ grid on;box on;
% $$$ % $$$ LabelAxes(gca,1,25,0.003,0.925);
% $$$ % $$$ 
% $$$ % $$$ subplot(2,1,2);
% $$$ % $$$ fields = {
% $$$ % $$$           {M(ind,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {VDFkppiw(ind,months,yrs), 'Background',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFkppish(ind,months,yrs), 'Shear Instability',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFkppbl(ind,months,yrs), 'KPP Boundary Layer',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFwave(ind,months,yrs), 'Topographic Internal Wave',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFwave(ind,months,yrs), 'Topographic Internal Wave',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFkppdd(ind,months,yrs)+VDFkppicon(ind,months,yrs)+KPPnloc(ind, ...
% $$$ % $$$                                                   months,yrs),'Other',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$ % $$$           {I(ind,months,yrs), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$ % $$$           {M(ind,months,yrs)+I(ind,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           };
% $$$ % $$$ leg = {};
% $$$ % $$$ legh = [];
% $$$ % $$$ for i=1:length(fields)
% $$$ % $$$     hold on;
% $$$ % $$$ % $$$     for j=1:length(P(1,1,:))
% $$$ % $$$ % $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$ % $$$ % $$$              'color',0.7*[1 1 1] ...
% $$$ % $$$ % $$$              ,'linewidth',0.5);
% $$$ % $$$ % $$$     end
% $$$ % $$$     legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
% $$$ % $$$          'color',fields{i}{3} ...
% $$$ % $$$          ,'linewidth',fields{i}{4});
% $$$ % $$$     hold on;
% $$$ % $$$     leg{i} = fields{i}{2};
% $$$ % $$$ end
% $$$ % $$$ xlabel('Month');
% $$$ % $$$ ylabel('PW');
% $$$ % $$$ lg = legend(legh,leg);
% $$$ % $$$ ylim([-2 0]);
% $$$ % $$$ xlim([1 tL]);
% $$$ % $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ % $$$ set(gca,'xtick',[1:tL]);
% $$$ % $$$ set(gca,'FontSize',25);
% $$$ % $$$ grid on;box on;
% $$$ % $$$ LabelAxes(gca,2,25,0.003,0.925);
% $$$ 
% $$$     %%% Meridional heat flux:
% $$$ 
% $$$ % Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = 94;
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em6/mat_data/';
% $$$ % $$$ model = 'MOM025_kb1em6';
% $$$ % $$$ outputs = 30;
% $$$ base = '/srv/ccrc/data03/z3500785/MOM01_HeatDiag/mat_data/';
% $$$ model = 'MOM01';
% $$$ outputs = [222];
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ % $$$ shfluxa = shflux;
% $$$ mhfluxa = mhflux;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$ % $$$     shfluxa = shfluxa+shflux;
% $$$     mhfluxa = mhfluxa+mhflux;
% $$$ end
% $$$ % $$$ shflux = shfluxa/length(outputs);
% $$$ % $$$ shflux = monmean(shflux,3,ndays);
% $$$ mhflux = mhfluxa/length(outputs);
% $$$ mhflux = monmean(mhflux,2,ndays);
% $$$ 
% $$$ % $$$ % Calculate meridional heat flux inferred:
% $$$ % $$$ latV = linspace(-90,90,181);
% $$$ % $$$ V = zeros(size(latV));
% $$$ % $$$ for i=1:length(latV)
% $$$ % $$$     inds = lat < latV(i);
% $$$ % $$$     V(i) = nansum(area(inds).*shflux(inds));
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ %Center the flux:
% $$$ % $$$ V = V + (V(1)-V(end))/2;
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[260         339        1055         586]);
% $$$ set(gcf,'defaulttextfontsize',25);
% $$$ set(gcf,'defaultaxesfontsize',25);
% $$$ % $$$ plot(latV,V/1e15,'-r','linewidth',2);
% $$$ hold on;
% $$$ plot(latv,mhflux/1e15,'-b','linewidth',2);
% $$$ 
% $$$ xlabel('Latitude ($^\circ$N)');
% $$$ ylabel('Meridional Heat Flux (PW)');
% $$$ grid on;
% $$$ box on;
% $$$ xlim([-90 90]);
% $$$ ylim([-1 2]);
% $$$ set(gca,'xtick',[-90:30:90]);
