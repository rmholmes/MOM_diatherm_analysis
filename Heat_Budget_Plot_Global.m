% This script makes plots of the heat budget in the MOM
% simulations.

% $$$ close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = { ...
% ACCESS-OM2 Gadi runs:
% $$$ % 1-degree
         {'ACCESS-OM2_1deg_jra55_ryf',[31]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_gfdl50',[31]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds75',[3135]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds100',[3135]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf_kds135',[3135]}, ...
% 1/4-degree
         {'ACCESS-OM2_025deg_jra55_ryf',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_norediGM',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_noGM',[7680]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5',[7781]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_rediGM_kbvar',[7781]}, ...
         {'ACCESS-OM2_025deg_jra55_ryf_kds75',[7680]}, ...
% 1/10-degree
         {'ACCESS-OM2_01deg_jra55_ryf',[636643]}, ...
         {'ACCESS-OM2_01deg_jra55_ryf_k_smag_iso7p9',[648:655]}, ...
       };
Inetstr = [];

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);

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
    Inetstr = [Inetstr sprintf(' Inet = %5.1f PWdegC, Mnet = %5.1f PWdegC',mean(monmean(Inet,2,ndays(1:length(Inet(1,:,1)))),3)/1e15,mean(monmean(Mnet,2,ndays(1:length(Mnet(1,:,1)))),3)/1e15) ' ' model ' \n '];

    %%%%Heat Flux: ---------------------------------------------------------------------------------------------
% Production fields:
    fields = { ...
        {N(:,months,yrs), 'Tendency $\partial\mathcal{H}_I/\partial t$','m',2,'-'}, ...
        {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
        {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
        {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
        {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
             };

% $$$     % Combination terms for Fig. 15:
% $$$     fields = { ...
% $$$         {F(:,months,yrs)+PI(:,months,yrs)-N(:,months,yrs), 'Forcing - Tendency (PW)','k',2,'-'}, ...
% $$$         {R(:,months,yrs)+M(:,months,yrs)+MD(:,months,yrs)+SG(:,months,yrs), 'Explicit Mixing (PW)','r',2,':'}, ...
% $$$         {I(:,months,yrs), 'Numerical Mixing (PW)','b',2,'--'}, ...
% $$$              };

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
