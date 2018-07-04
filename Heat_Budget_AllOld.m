% This script makes plots of the heat budget in the MOM
% simulations.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[222]}, ...
% $$$ % MOM025-SIS:
% $$$     {'MOM025',[8:12]}, ...
% $$$     {'MOM025',[15:19]}, ...
% $$$     {'MOM025_kb1em6',[30]}, ...
% $$$     {'MOM025_kb3seg',[80:84]}, ...
% $$$     {'MOM025_kb1em5',[94]}, ...
% $$$     {'MOM025_wombat',[1978]}, ...
% ACCESS-OM2 025-degree:
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485',[78]}, ...
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485_redi',[59]}, ...
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485_gmredi',[73]}, ...
% $$$ %     {'ACCESS-OM2_025deg_jra55_ryf8485_KDS75',[??]}, ...
% ACCESS-OM2 1-degree:
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_Tcen',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_TcenGMS',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_gfdl50_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds75_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds100_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds135_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_kb1em5',[0]}, ...
       };

rr = 1;
% $$$ for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    ndays = diff(time_snap);
    ndays = ndays(1:12);
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
% $$$ region = 'Pacific';
    nyrs = tL/12;szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
    yrs = 1:nyrs;
    months = 1:12;
    
    ycur = 1;

% $$$     %% Global Calculations:
% $$$     for i=1:length(outputs)
% $$$         
% $$$ % $$$     % Annual or Monthly offline Binning:
% $$$ % $$$     load([base model sprintf('_output%03d_',outputs(i)) 'GlobalHBud_MonAnBin.mat']);
% $$$ % $$$     GWB = GWBann;
% $$$ 
% $$$         load([base model sprintf('_output%03d_',outputs(i)) region 'HBud.mat']);
% $$$         
% $$$         % Fluxes:
% $$$         P(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.PME+GWB.RMX,szTe); % PME effective heat flux (W)
% $$$         F(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH+GWB.VDS+GWB.FRZ+GWB.ETS,szTe); % Surface heat flux (W)
% $$$         M(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDF+GWB.KNL,szTe); % Vertical mixing flux (W)
% $$$         if (isfield(GWB,'VDFkppiw')) % Vertical mixing components
% $$$             VDFkppiw(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw,szTe);
% $$$             VDFkppish(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppish,szTe);
% $$$             VDFkppicon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppicon,szTe);
% $$$             VDFkppbl(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppbl,szTe);
% $$$             VDFkppdd(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppdd,szTe);
% $$$             VDFwave(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFwave,szTe);
% $$$             KPPnloc(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.KNL,szTe);
% $$$             VDFsum(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.VDFkppiw+GWB.VDFkppish+GWB.VDFkppicon+ ...
% $$$                 GWB.VDFkppbl+GWB.VDFkppdd+GWB.VDFwave+GWB.KNL,szTe);
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
% $$$             MD(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.MDS,szTe);; % GM (W)
% $$$             M(:,:,ycur:(ycur+nyrs-1)) = M(:,:,ycur:(ycur+nyrs-1)) + reshape(GWB.MDS,szTe); %ADD TO VERTICAL MIXING, but it's small...
% $$$         else
% $$$             MD(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'NUM')) % Pre-calculated numerical mixing
% $$$             NUM(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUM,szTe); % NUM (W)
% $$$         else
% $$$             NUM(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         if (isfield(GWB,'NUMH')) % Pre-calculated numerical mixing
% $$$                                 % from heat budget
% $$$             NUMH(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.NUMH,szTe); % NUM (W)
% $$$         else
% $$$             NUMH(:,:,ycur:(ycur+nyrs-1)) = zeros(size(P(:,:,ycur:(ycur+nyrs-1))));
% $$$         end    
% $$$         D(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN-GWB.ADV-GWB.SUB,szTe)-GM(:,:,ycur:(ycur+nyrs-1)); % Material derivative of T (W)
% $$$         SW(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SWH,szTe); % Short-wave heat
% $$$         JS(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SFW,szTe); % Surface Volume Flux
% $$$         SUB(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.SUB,szTe);
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
% $$$         JSH(:,:,ycur:(ycur+nyrs-1)) = JS(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;
% $$$ 
% $$$         % Interior heat source P:
% $$$         PI(:,:,ycur:(ycur+nyrs-1)) = P(:,:,ycur:(ycur+nyrs-1)) - JSH(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Interior heat source Q:
% $$$         QII(:,:,ycur:(ycur+nyrs-1)) = QI(:,:,ycur:(ycur+nyrs-1)) - JI(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;
% $$$ 
% $$$         % Across-isotherm advective heat flux:
% $$$         CIA(:,:,ycur:(ycur+nyrs-1)) = G(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;
% $$$ 
% $$$         % External HC Tendency:
% $$$         EHC(:,:,ycur:(ycur+nyrs-1)) = dVdt(:,:,ycur:(ycur+nyrs-1)).*repmat(Te,[1 12 nyrs])*rho0*Cp;
% $$$ 
% $$$         % Internal HC Tendency:
% $$$         N(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1)) - EHC(:,:,ycur:(ycur+nyrs-1));
% $$$ % $$$ N(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TEN,szTe);
% $$$ 
% $$$ % $$$ % Alternative method 1 for the N calculation:
% $$$ % $$$ N(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*cumsum(dVdt*dT,1,'reverse');
% $$$ % However, this does not work well, potentially because the
% $$$ % integral should conceptually then be defined on Tcenters as
% $$$ % opposed to Tedges. In any case, this method gives a non-zero
% $$$ % total heat flux due to implicit mixing and is much noisier. 
% $$$ 
% $$$         % Implicit mixing:
% $$$         I(:,:,ycur:(ycur+nyrs-1)) = N(:,:,ycur:(ycur+nyrs-1)) - F(:,:,ycur:(ycur+nyrs-1)) - P(:,:,ycur:(ycur+nyrs-1)) - M(:,:,ycur:(ycur+nyrs-1)) - R(:,:,ycur:(ycur+nyrs-1)) + JSH(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Non-advective flux into volume:
% $$$         B(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+I(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Monthly binned Internal HC Tendency:
% $$$         if (isfield('GWB','TENMON'))
% $$$             Nmon(:,:,ycur:(ycur+nyrs-1)) = reshape(GWB.TENMON,szTe);
% $$$         end
% $$$ 
% $$$         % Alternative method 2 (which is what I was using for
% $$$         % the GRL submit - but the above method is simpler to explain). 
% $$$ % $$$ Ialt(:,:,ycur:(ycur+nyrs-1)) = dHdt(:,:,ycur:(ycur+nyrs-1))-D(:,:,ycur:(ycur+nyrs-1))-CIA(:,:,ycur:(ycur+nyrs-1)) + QI(:,:,ycur:(ycur+nyrs-1));
% $$$ % $$$ Balt(:,:,ycur:(ycur+nyrs-1)) = F(:,:,ycur:(ycur+nyrs-1))+M(:,:,ycur:(ycur+nyrs-1))+Ialt(:,:,ycur:(ycur+nyrs-1))+R(:,:,ycur:(ycur+nyrs-1));
% $$$ % $$$ Nalt(:,:,ycur:(ycur+nyrs-1)) = Balt(:,:,ycur:(ycur+nyrs-1)) + PI(:,:,ycur:(ycur+nyrs-1)) - QII(:,:,ycur:(ycur+nyrs-1));
% $$$ % Gives very close results. The difference between N and Nalt, and
% $$$ % I and Ialt (i.e. (Ialt-I)./I, (Nalt-N)./N) is less than 1e-5 in
% $$$ % all cases (except where I==0, which gives Inf).
% $$$ 
% $$$         % WMT from B:
% $$$         WMTM(:,:,ycur:(ycur+nyrs-1)) = -diff(M(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTF(:,:,ycur:(ycur+nyrs-1)) = -diff(F(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTI(:,:,ycur:(ycur+nyrs-1)) = -diff(I(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMTR(:,:,ycur:(ycur+nyrs-1)) = -diff(R(:,:,ycur:(ycur+nyrs-1)),[],1)/dT/rho0/Cp;
% $$$         WMT(:,:,ycur:(ycur+nyrs-1)) = WMTM(:,:,ycur:(ycur+nyrs-1))+WMTF(:,:,ycur:(ycur+nyrs-1))+WMTI(:,:,ycur:(ycur+nyrs-1))+WMTR(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % WMT HB from B:
% $$$         HWMTM(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTM(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
% $$$         HWMTF(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTF(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
% $$$         HWMTI(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTI(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
% $$$         HWMTR(:,:,ycur:(ycur+nyrs-1)) = rho0*Cp*WMTR(:,:,ycur:(ycur+nyrs-1)).*repmat(T,[1 12 nyrs]);
% $$$         HWMT(:,:,ycur:(ycur+nyrs-1)) = HWMTM(:,:,ycur:(ycur+nyrs-1))+HWMTF(:,:,ycur:(ycur+nyrs-1))+HWMTI(:,:,ycur:(ycur+nyrs-1))+HWMTR(:,:,ycur:(ycur+nyrs-1));
% $$$ 
% $$$         % Alternative method 3 I from Volume budget (as for spatial structure calc FlI):
% $$$ % $$$ WMTI(:,:,ycur:(ycur+nyrs-1)) = avg(dVdt(:,:,ycur:(ycur+nyrs-1)),1) - avg(JS(:,:,ycur:(ycur+nyrs-1)),1)-WMTM(:,:,ycur:(ycur+nyrs-1))-WMTF(:,:,ycur:(ycur+nyrs-1))-WMTR(:,:,ycur:(ycur+nyrs-1));
% $$$ % $$$ I(:,:,ycur:(ycur+nyrs-1)) = zeros(size(I(:,:,ycur:(ycur+nyrs-1))));
% $$$ % $$$ I(1:(end-1),:,ycur:(ycur+nyrs-1)) = rho0*Cp*cumsum(WMTI(:,:,ycur:(ycur+nyrs-1))*dT,1,'reverse');
% $$$ % This also gives pretty consistent results, although has similar
% $$$ % noisy problems and a non-zero total heat flux at the cooler
% $$$ % temperatures as alternative method 1 above. This method is
% $$$ % slightly better at the warmest temperatures.
% $$$         ycur = ycur+nyrs;
% $$$     end
% $$$     months = [1:length(P(1,:,1))];
% $$$     yrs = [1:length(P(1,1,:))];

% $$$     yrs = [1 5];
% $$$
% $$$
% $$$ %%%%Heat Flux: ---------------------------------------------------------------------------------------------
% $$$ % Production fields:
% $$$ fields = { ...
% $$$           {N(:,months,yrs), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
% $$$           {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$ % $$$           {M(:,months,yrs)+R(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{R}+\mathcal{I}$','r',2,'-'}, ...
% $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
% $$$           {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$           {NUM(:,months,yrs), 'Numerical Mixing Direct $\mathcal{I}$','c',2,'-'}, ...
% $$$ % $$$           {SW(:,months,yrs), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$ % $$$           {dHdt(:,months,yrs), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',2,'--'}, ...
% $$$ % $$$           {PI(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}_I$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {F(:,months,yrs), 'Surface Heat Fluxes $\mathcal{F}$','k',2,'-'}, ...
% $$$ % $$$           {P(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}$',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$           {GM(:,months,yrs), 'GM $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {SUB(:,months,yrs), 'SUB $\mathcal{S}$',[0 0.5 0],2,':'}, ...
% $$$ % $$$           {HWMTI(:,months,yrs), 'Advective Implicit Mixing','b',2,'--'}, ...
% $$$ % $$$           {HWMTM(:,months,yrs), 'Advective Vertical Mixing','r',2,'--'}, ...
% $$$ % $$$           {HWMTF(:,months,yrs), 'Advective Surface Forcing','k',2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {Nmon(:,months,yrs), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$ % $$$           {CIA(:,months,yrs), 'Across-Isotherm Advection $\mathcal{G}\Theta\rho_0C_p$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           };
% $$$ % $$$ fields = { ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {VDFkppiw(:,months,yrs), 'Vertical Diffusion KPPIW','b',2,'-'}, ...
% $$$ % $$$           {VDFkppish(:,months,yrs), 'Vertical Diffusion KPPISH',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {VDFkppbl(:,months,yrs), 'Vertical Diffusion KPPBL','m',2,'-'}, ...
% $$$ % $$$           {VDFwave(:,months,yrs), 'Vertical Diffusion WAVE','k',2,'-'}, ...
% $$$ % $$$ % $$$           {VDFkppicon(:,months,yrs), 'Vertical Diffusion KPPICON','y',2,'-'}, ...
% $$$ % $$$ % $$$           {KPPnloc(:,months,yrs), 'KPP Non-local',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$ % $$$           {VDFkppdd(:,months,yrs), 'Vertical Diffusion KPPDD','m',2,'-'}, ...
% $$$ % $$$           {VDFkppdd(:,months,yrs)+VDFkppicon(:,months,yrs)+KPPnloc(:, ...
% $$$ % $$$                                                   months,yrs),'DD + KPP non-local + Int. Convection',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$ % $$$           {VDFsum(:,months,yrs), 'Vertical Mixing SUM','r',2,'--'}, ...
% $$$ % $$$           };
% $$$ 
% $$$ Fscale = 1/1e15;
% $$$ 
% $$$ yrtyps = {'-','-','--','-.',':'}; % line-types for different years
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
% $$$     
% $$$     % Plot years from a single run separately:
% $$$ % $$$     for j=1:length(yrs) 
% $$$ % $$$         h = plot(Te,monmean(fields{i}{1}(:,:,yrs(j)),2,ndays(months))*Fscale,yrtyps{j}, 'color',fields{i}{3} ...
% $$$ % $$$              ,'linewidth',3);
% $$$ % $$$         if (j == 1)
% $$$ % $$$             legh(i) = h;
% $$$ % $$$         end
% $$$ % $$$     end
% $$$ % $$$     leg{i} = fields{i}{2};
% $$$ 
% $$$     % Average years together for a single run:
% $$$ % $$$     legh(i) = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),fields{i}{5}, 'color',fields{i}{3} ...
% $$$ % $$$          ,'linewidth',fields{i}{4});
% $$$ % $$$     leg{i} = fields{i}{2};
% $$$     
% $$$     % Average years together for multiple runs:
% $$$     tmp = plot(x,mean(monmean(fields{i}{1},2,ndays(months))*Fscale,3),yrtyps{rr}, 'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$ % $$$     if i==1
% $$$ % $$$         leg{rr} = strrep(RUNS{rr}{1},'_',' ');
% $$$ % $$$         legh(rr) = tmp;
% $$$ % $$$     end
% $$$ end
% $$$ ylim([-1.5 1.5]);
% $$$ xlim([-3 31]);
% $$$ box on; 
% $$$ grid on;
% $$$ ylabel('Heat flux into fluid warmer than $\Theta$ (PW)');
% $$$ xlabel('Temperature $\Theta$ ($^\circ$C)');
% $$$ % $$$ lg = legend(legh,leg);
% $$$ % $$$ set(lg,'Position',[0.5881    0.5500    0.2041    0.2588]);
% $$$ 
% $$$ end

% $$$ end

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
% $$$ % $$$           {N(:,months,yrs), 'Internal HC Tendency $\mathcal{N}$','m',2,'-'}, ...
% $$$ % $$$           {F(:,months,yrs)+PI(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}_I$','k',2,'-'}, ...
% $$$ % $$$           {F(:,months,yrs), 'Surface Heat Fluxes $\mathcal{F}$','k',2,'-'}, ...
% $$$ % $$$           {P(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}$',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$           {PI(:,months,yrs), 'Surface Volume Fluxes $\mathcal{P}_I$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(:,months,yrs), 'Numerical Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$           {R(:,months,yrs), 'Redi Mixing $\mathcal{R}$',[0 0.5 0],2,'-'}, ...
% $$$ % $$$           {GM(:,months,yrs), 'GM $\mathcal{G}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {SUB(:,months,yrs), 'SUB $\mathcal{S}$',[0 0.5 0],2,':'}, ...
% $$$ % $$$           {HWMTI(:,months,yrs), 'Advective Implicit Mixing','b',2,'--'}, ...
% $$$ % $$$           {HWMTM(:,months,yrs), 'Advective Vertical Mixing','r',2,'--'}, ...
% $$$ % $$$           {HWMTF(:,months,yrs), 'Advective Surface Forcing','k',2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {Nmon(:,months,yrs), 'Monthly-Binned Total','m',2,'--'}, ...
% $$$ % $$$           {SW(:,months,yrs), 'Shortwave Redistribution',0.5*[1 1 1],2,'--'}, ...
% $$$ % $$$           {dHdt(:,months,yrs), 'HC Tendency $\frac{\partial\mathcal{H}}{\partial t}$','m',2,'--'}, ...
% $$$ % $$$           {CIA(:,months,yrs), 'Across-Isotherm Advection $\mathcal{G}\Theta\rho_0C_p$',[0.49 0.18 0.56],2,'--'}, ...
% $$$ % $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$ % $$$           {VDFtotal(:,months,yrs), 'Vertical Diffusion','r',2,'-'}, ...
% $$$ % $$$           {VDFkppiw(:,months,yrs), 'Background','b',2,'-'}, ...
% $$$ % $$$           {VDFkppish(:,months,yrs), 'Shear Instability','g',2,'-'}, ...
% $$$ % $$$           {VDFkppbl(:,months,yrs), 'KPP Boundary Layer','c',2,'-'}, ...
% $$$ % $$$           {VDFwave(:,months,yrs), 'Vertical Diffusion WAVE','k',2,'-'}, ...
% $$$ % $$$           {VDFkppicon(:,months,yrs), 'Vertical Diffusion KPPICON','y',2,'-'}, ...
% $$$ % $$$           {VDFkppdd(:,months,yrs), 'Vertical Diffusion KPPDD','m',2,'-'}, ...
% $$$ % $$$           {VDFsum(:,months,yrs), 'Vertical Diffusion SUM','r',2,'--'}, ...
% $$$ % $$$           {KPPnloc(:,months,yrs), 'KPP Non-local','m',2,'--'}, ...
% $$$           };
% $$$ 
% $$$ % Fluxes:
% $$$ scale = 1/1e15;label = '(PW)';x = Te;
% $$$ caxs = [-1.5 0];
% $$$ sp = 0.05;
% $$$ caxs = [-5 5];
% $$$ sp = 0.5;
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
% $$$ for ii=1:length(fields)
% $$$ % $$$     subplot(1,length(fields),ii);
% $$$     subplot_tight(3,4,rr);%1,length(fields),ii);
% $$$     V = mean(fields{ii}{1},3)'*scale;
% $$$     if (length(fields{ii}{1}(:,1)) == length(Te))
% $$$         x = Te;
% $$$     else
% $$$         x = T;
% $$$     end
% $$$     [X,Y] = ndgrid(1:tL,x);
% $$$     contourf(X,Y,V,cint);%,'linestyle','none');
% $$$ % $$$     cb = colorbar('Location','NorthOutside','FontSize',25);    
% $$$     set(gca,'ytick',-5:5:35);
% $$$     set(gca,'xtick',[1:tL]);
% $$$     if (rr>=9)
% $$$         set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun', ...
% $$$                             'Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$     else
% $$$         set(gca,'xticklabel',[]);
% $$$     end
% $$$     ylim([-3 31]);
% $$$     grid on;
% $$$     caxis(caxs);
% $$$ % $$$     xlabel('Month');
% $$$ % $$$     ylabel('Temperature ($^\circ$C)');
% $$$ % $$$     xlabel(cb,strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
% $$$ % $$$     xlabel(cb,[model ' ' fields{ii}{2} ' ' ...
% $$$ % $$$                label],'FontSize',20);
% $$$     title(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''));
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
% $$$ 
% $$$ %% Global Seasonal Cycle TS
% $$$ months = 1:12;
% $$$ Ts = 21.5;
% $$$ [tmp ind] = min(abs(Te-Ts));
% $$$ 
% $$$ Fscale = 1/1e15;
% $$$ 
% $$$ figure;
% $$$ %set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'Position',[34          40        1164         963]);
% $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ 
% $$$ 
% $$$ subplot(2,1,1);
% $$$ fields = {
% $$$           {F(ind,months,yrs)+PI(ind,months,yrs), 'Surface Forcing $\mathcal{F}$','k',2,'-'}, ...
% $$$           {N(ind,months,yrs), 'Total $\mathcal{N}$','m',2,'-'}, ...
% $$$           {M(ind,months,yrs)+I(ind,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           };
% $$$ 
% $$$ leg = {};
% $$$ legh = [];
% $$$ for i=1:length(fields)
% $$$     hold on;
% $$$ % $$$     for j=1:length(P(1,1,:))
% $$$ % $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$ % $$$              'color',0.7*[1 1 1] ...
% $$$ % $$$              ,'linewidth',0.5);
% $$$ % $$$     end
% $$$     legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
% $$$          'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$     hold on;
% $$$     leg{i} = fields{i}{2};
% $$$ end
% $$$ xlabel('Month');
% $$$ ylabel('PW');
% $$$ lg = legend(legh,leg);
% $$$ ylim([-5 5]);
% $$$ xlim([1 tL]);
% $$$ set(gca,'xtick',[1:tL]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'FontSize',25);
% $$$ grid on;box on;
% $$$ LabelAxes(gca,1,25,0.003,0.925);
% $$$ 
% $$$ subplot(2,1,2);
% $$$ fields = {
% $$$           {M(ind,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {VDFkppiw(ind,months,yrs), 'Background',[0 0.5 0],2,'-'}, ...
% $$$           {VDFkppish(ind,months,yrs), 'Shear Instability',[0 0.5 0],2,'-'}, ...
% $$$           {VDFkppbl(ind,months,yrs), 'KPP Boundary Layer',[0 0.5 0],2,'-'}, ...
% $$$           {VDFwave(ind,months,yrs), 'Topographic Internal Wave',[0 0.5 0],2,'-'}, ...
% $$$           {VDFwave(ind,months,yrs), 'Topographic Internal Wave',[0 0.5 0],2,'-'}, ...
% $$$           {VDFkppdd(ind,months,yrs)+VDFkppicon(ind,months,yrs)+KPPnloc(ind, ...
% $$$                                                   months,yrs),'Other',[0.49 0.18 0.56],2,'-'}, ...
% $$$ % $$$           {I(ind,months,yrs), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$           {M(ind,months,yrs)+I(ind,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$           };
% $$$ leg = {};
% $$$ legh = [];
% $$$ for i=1:length(fields)
% $$$     hold on;
% $$$ % $$$     for j=1:length(P(1,1,:))
% $$$ % $$$         plot(1:tL,fields{i}{1}(:,:,j)*Fscale,fields{i}{5}, ...
% $$$ % $$$              'color',0.7*[1 1 1] ...
% $$$ % $$$              ,'linewidth',0.5);
% $$$ % $$$     end
% $$$     legh(i) = plot(1:tL,mean(fields{i}{1},3)*Fscale,fields{i}{5}, ...
% $$$          'color',fields{i}{3} ...
% $$$          ,'linewidth',fields{i}{4});
% $$$     hold on;
% $$$     leg{i} = fields{i}{2};
% $$$ end
% $$$ xlabel('Month');
% $$$ ylabel('PW');
% $$$ lg = legend(legh,leg);
% $$$ ylim([-2 0]);
% $$$ xlim([1 tL]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'xtick',[1:tL]);
% $$$ set(gca,'FontSize',25);
% $$$ grid on;box on;
% $$$ LabelAxes(gca,2,25,0.003,0.925);


    %%% Spatial Structure:

    VAR = 'FlM';
    TYPE = 'VertInt';
    Tl = 22.5;
    name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
    eval(['load(name,''' VAR ''');']);
    eval([VAR '(isnan(' VAR ')) = 0.0;']);
    if (length(outputs)==1)
        eval([VAR ' = reshape(' VAR ',[length(' VAR '(:,1,1)) length(' VAR '(1,:,1)) 12 nyrs]);']);
        eval([VAR ' = mean(' VAR '(:,:,:,yrs),4);']);
    else
        eval([VAR 'a = ' VAR ';']);
        for i=2:length(outputs)
            name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
            eval(['load(name,''' VAR ''');']);
            eval([VAR '(isnan(' VAR ')) = 0.0;']);
            eval([VAR 'a = ' VAR 'a + ' VAR ';']);
        end
        eval([VAR ' = ' VAR 'a/length(outputs);']);
        eval([VAR '(' VAR '==0) = NaN;']);
    end
    eval(['FlM = ' VAR ';']);
% $$$ 
% $$$ FlMblall = FlM - FlMwave - FlMkppish;
% $$$ FlM = FlMblall;
% $$$ 
% $$$ % $$$ % CHECK spatial structure sums to total:
% $$$ % $$$ % $$$ Tls = [14.75:2.5:27.25]+0.25;
% $$$ % $$$ Tls = [5 10 15:2.5:27.5]-0.25;
% $$$ % $$$ SUM = zeros(size(Tls));
% $$$ % $$$ for ii = 1:length(Tls)
% $$$ % $$$ 
% $$$ % $$$     Tl = Tls(ii)
% $$$ % $$$     load([base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']);
% $$$ % $$$     eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$ % $$$     eval([VAR 'a = ' VAR ';']);
% $$$ % $$$     for i=2:length(outputs)
% $$$ % $$$         load([base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']);
% $$$ % $$$         eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$ % $$$         eval([VAR 'a = ' VAR 'a + ' VAR ';']);
% $$$ % $$$     end
% $$$ % $$$     eval([VAR ' = ' VAR 'a/length(outputs);']);
% $$$ % $$$     eval([VAR '(' VAR '==0) = NaN;']);
% $$$ % $$$     eval(['FlM = ' VAR ';']);
% $$$ % $$$     tmp = FlM;
% $$$ % $$$     tmp(isnan(tmp)) = 0.0;
% $$$ % $$$     Z = monmean(tmp(:,:,months),3,ndays(months));
% $$$ % $$$     Z(Z == 0) = NaN;
% $$$ % $$$     SUM(ii) = nansum(nansum(area.*Z));
% $$$ % $$$ end
% $$$ %plot(Tls,SUM/1e6,'Xb','MarkerSize',12,'LineWidth',2); 
% $$$ %plot(Tls,SUM/1e15,'Xb','MarkerSize',12,'LineWidth',2); 
% $$$ 
% $$$ % $$$ %%% Regional time series 
% $$$ % $$$ % $$$ 
% $$$ % $$$ %[mask_t,~] = Heat_Budget_Mask('Pacific','','','',base,'MOM025');
% $$$ % $$$ mask_t = zeros(size(lon));
% $$$ % $$$ months = [1:12];
% $$$ % $$$ %Region choice:
% $$$ % $$$ regions = { ...
% $$$ % $$$     {'Global',lat==lat,'k'}, ...
% $$$ % $$$     {'Equatorial',abs(lat)<=10,'k'}, ...
% $$$ % $$$     {'Eastern Pacific',lat > -10 & lat < 10 & lon > -160 & lon<-70,[0.4667    0.6745    0.1882]}, ...
% $$$ % $$$     {'Kuroshio', lat < 50 & lat > 10 & lon > -260 & lon < -140,'r'}, ...
% $$$ % $$$     {'Gulf Stream', lat < 50 & lat > 10 & lon > -100 & lon < 0 & ~mask_t,'b'}, ...
% $$$ % $$$     {'Indian', lat < 25 & lat > -55 & (lon > 20 | lon < -260),[0.4941    0.1843    0.5569]}, ...
% $$$ % $$$     {'South Pacific/Atlantic', lat < -10 & lat > -55 & lon > -260 & lon < 20,'g'}, ...
% $$$ % $$$           };
% $$$ % $$$ % $$$     {'Eastern Pacific',lat > -20 & lat < 20 & lon > -150 & lon<-60 & mask_t,[0.4667    0.6745    0.1882]}, ...
% $$$ % $$$ % $$$     {'WWV',abs(lat)<=5 & lon>
% $$$ % $$$ % $$$     {'NH',lat>10,'k'}, ...
% $$$ % $$$ % $$$     {'SH',lat<-10,'k'}, ...
% $$$ % $$$ % $$$     {'Outside Equatorial',lat<-10 | lat > 10,'k'}, ...
% $$$ % $$$ % $$$     {'Nino 3',lat > -5 & lat < 5 & lon > -150 & lon < -90,'m'}, ...
% $$$ % $$$ % $$$     {'Eastern Pacific',lat > -10 & lat < 10 & lon > -160 & mask_t,'c'}, ...
% $$$ % $$$ % $$$     {'Equatorial Atlantic',lat > -10 & lat < 10 & lon > -65 & lon < 20,[0 0.5 0]}, ...
% $$$ % $$$ % $$$     {'Western Pacific',lat > -10 & lat < 10 & lon > -260 & lon < -160,[0.5 0.5 0.5]}, ...
% $$$ % $$$ % $$$     {'Gulf Stream', lat < 45 & lat > 10 & lon > -100 & lon < -40 & ~mask_t,'b'}, ...
% $$$ % $$$ 
% $$$ % $$$ Field = zeros(12,length(regions));
% $$$ % $$$ AREA = zeros(12,length(regions));
% $$$ % $$$ for i=1:12
% $$$ % $$$     Mtmp = FlM(:,:,i);
% $$$ % $$$     Atmp = area;
% $$$ % $$$     Atmp(isnan(Mtmp)) = NaN;
% $$$ % $$$     
% $$$ % $$$     for ii=1:length(regions)
% $$$ % $$$         Field(i,ii) = nansum(nansum(area(regions{ii}{2}).*Mtmp(regions{ii}{2}),1),2);
% $$$ % $$$         AREA(i,ii) = nansum(nansum(Atmp(regions{ii}{2}),1),2);
% $$$ % $$$     end
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ %Display output in terminal for table:
% $$$ % $$$ str = {['Area fractions model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$        ' '};
% $$$ % $$$ for ii = 1:length(regions)
% $$$ % $$$     str{ii+2} = sprintf([regions{ii}{1} ' area = %3.2f'], ...
% $$$ % $$$                         monmean(AREA(:,ii),1,ndays(months))/monmean(AREA(:,1),1,ndays(months)));
% $$$ % $$$ end
% $$$ % $$$ str
% $$$ % $$$ 
% $$$ % $$$ str = {['Annual totals model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$        ' '};
% $$$ % $$$ for ii = 1:length(regions)
% $$$ % $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.2fPW (%3.0f)', ...
% $$$ % $$$                         monmean(Field(:,ii),1,ndays(months))/1e15,monmean(Field(:,ii),1,ndays(months))/monmean(Field(:,1),1,ndays(months))*100)];
% $$$ % $$$ end
% $$$ % $$$ str
% $$$ % $$$ 
% $$$ % $$$ str = {['SC range model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$        ' '};
% $$$ % $$$ for ii = 1:length(regions)
% $$$ % $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.2fPW (%3.0f)', ...
% $$$ % $$$                         (max(Field(:,ii))-min(Field(:,ii)))/1e15,(max(Field(:,ii))-min(Field(:,ii)))/(max(Field(:,1))-min(Field(:,1)))*100)];
% $$$ % $$$ end
% $$$ % $$$ str
% $$$ % $$$ 
% $$$ % $$$ % $$$ mn1 = 4;mn2 = 7;
% $$$ % $$$ % $$$ str = {['SC range Apr-Jul model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$ % $$$        ' '};
% $$$ % $$$ mn1 = 5;mn2 = 8;
% $$$ % $$$ str = {['SC range May-Aug model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$        ' '};
% $$$ % $$$ for ii = 1:length(regions)
% $$$ % $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.5fPW (%3.0f)', ...
% $$$ % $$$                         (Field(mn1,ii)-Field(mn2,ii))/1e15,(Field(mn1,ii)-Field(mn2,ii))/(Field(mn1,1)-Field(mn2,1))*100)];
% $$$ % $$$ end
% $$$ % $$$ str
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ % $$$ % The 1% (within Nino 3) area count of the annual mean:
% $$$ % $$$ % $$$ inds = lat > -5 & lat < 5 & lon > -150 & lon < -90;
% $$$ % $$$ % $$$ M1p = zeros(length(find(inds)),12);
% $$$ % $$$ % $$$ A1p = area(inds);
% $$$ % $$$ % $$$ AREAtotal = zeros(12,1);
% $$$ % $$$ % $$$ for i=1:12
% $$$ % $$$ % $$$     Mtmp = FlM(:,:,i);
% $$$ % $$$ % $$$     M1p(:,i) = A1p.*Mtmp(inds);
% $$$ % $$$ % $$$     Atmp = area;
% $$$ % $$$ % $$$     Atmp(isnan(Mtmp)) = NaN;
% $$$ % $$$ % $$$     AREAtotal(i) = nansum(nansum(Atmp,1),2);
% $$$ % $$$ % $$$ end
% $$$ % $$$ % $$$ M1pA = monmean(M1p,2,ndays(months));
% $$$ % $$$ % $$$ [M1pA,I] = sort(M1pA);
% $$$ % $$$ % $$$ M1p = M1p(I,:);
% $$$ % $$$ % $$$ Atotal = monmean(AREAtotal,1,ndays(months));
% $$$ % $$$ % $$$ [tmp ind] = min(abs(cumsum(A1p(I))/Atotal - 0.005));
% $$$ % $$$ % $$$ %plot(cumsum(A1p(I))/Atotal*100,cumsum(M1pA)/monmean(MAll,2,ndays(months))*100)
% $$$ % $$$ % $$$ M1pA = sum(M1pA(1:ind));
% $$$ % $$$ % $$$ M1pSCR = (sum(M1p(1:ind,mn1))-sum(M1p(1:ind,mn2)));
% $$$ % $$$ % $$$ MAll = Field(:,1)';
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ str = {['1% area within Nino 3 model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$ % $$$ % $$$        sprintf(' 1p Annual-Mean = %3.2fPW (%3.0f)',M1pA/1e15,M1pA/monmean(MAll,2,ndays(months))*100); ...
% $$$ % $$$ % $$$        sprintf(' 1p SC Apr-Jul = %3.2fPW (%3.0f)',M1pSCR/1e15,M1pSCR/((MAll(mn1))-(MAll(mn2)))*100)}
% $$$ % $$$ 
% $$$ % $$$ % $$$ % Average fluxes, isotherm separationl, diffusivity:
% $$$ % $$$ % $$$ load([base 'MOM025_output002_Ziso_T23C.mat']);
% $$$ % $$$ % $$$ topiso = ziso;
% $$$ % $$$ % $$$ load([base 'MOM025_output002_Ziso_T22C.mat']);
% $$$ % $$$ % $$$ botiso = ziso;
% $$$ % $$$ % $$$ AvgFlux = MAll./AREAtotal;
% $$$ % $$$ % $$$ NaNs = isnan(topiso-botiso);
% $$$ % $$$ % $$$ AvgDZ = zeros(12,1);
% $$$ % $$$ % $$$ for ii=1:12
% $$$ % $$$ % $$$     AvgDZ(ii) = nansum(nansum(((topiso(:,:,ii)-botiso(:,:,ii)).*area),1),2)./ ...
% $$$ % $$$ % $$$                 nansum(area(~isnan((topiso(:,:,ii)-botiso(:,:,ii)))));
% $$$ % $$$ % $$$ end
% $$$ % $$$ % $$$ AvgDiff = monmean(AvgFlux/rho0/Cp.*AvgDZ',2,ndays(months))
% $$$ % $$$ % $$$ AvgFlux = monmean(AvgFlux,2,ndays(months))
% $$$ % $$$ % $$$ AvgDZ = monmean(AvgDZ',2,ndays(months))
% $$$ % $$$ % $$$ sprintf('Average flux across the isotherm = %3.1f Wm-2',AvgFlux)
% $$$ % $$$ % $$$ sprintf('Average diffusivity = %3.5f m2s-1',AvgDiff)
% $$$ % $$$ % $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ % $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ % $$$ 
% $$$ % $$$ legstr = {};
% $$$ % $$$ for ii=1:length(regions)
% $$$ % $$$     plot(1:12,Field(:,ii)/1e15,'-','color',regions{ii}{3},'LineWidth',2);
% $$$ % $$$     hold on;
% $$$ % $$$     legstr{ii} = regions{ii}{1};
% $$$ % $$$ end
% $$$ % $$$ % $$$ plot(1:12,MEEP/1e15,'--r','LineWidth',4);
% $$$ % $$$ % $$$ hold on;
% $$$ % $$$ % $$$ plot(1:12,MEq/1e15,'--k','LineWidth',4);
% $$$ % $$$ % $$$ plot(1:12,MN/1e15,':b','LineWidth',4);
% $$$ % $$$ % $$$ plot(1:12,MS/1e15,':','color',[0 0.5 0],'LineWidth',4);
% $$$ % $$$ % $$$ plot(1:12,(MN+MS)/1e15,':k','LineWidth',4);
% $$$ % $$$ % $$$ plot(1:12,MAll/1e15,'-k','LineWidth',4);
% $$$ % $$$ xlabel('Month');
% $$$ % $$$ ylabel(['PW']);
% $$$ % $$$ % $$$ title(['Vertical mixing heat flux through $' num2str(Tl) '^\circ$C isotherm']);
% $$$ % $$$ title(['Implicit mixing heat flux through $' num2str(Tl) '^\circ$C isotherm']);
% $$$ % $$$ 
% $$$ % $$$ leg = legend(legstr);
% $$$ % $$$ % $$$ 'Eastern Equatorial Pacific', ...
% $$$ % $$$ % $$$              'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
% $$$ % $$$ % $$$              ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
% $$$ % $$$ % $$$              ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
% $$$ % $$$ % $$$              ['Total']);
% $$$ % $$$ set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
% $$$ % $$$ ylim([-0.8 0]);
% $$$ % $$$ xlim([1 12]);
% $$$ % $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ % $$$ set(gca,'xtick',[1:12]);
% $$$ % $$$ set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
% $$$ % $$$ set(gca,'FontSize',25);
% $$$ % $$$ grid on;
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ axes('Position',[0.74 0.2451 0.13 0.6799]);
% $$$ % $$$ % $$$ plot([0 1],[min(MAll/1e15) min(MAll/1e15)],'-k','linewidth',2);
% $$$ % $$$ % $$$ hold on;
% $$$ % $$$ % $$$ plot([0 1],[max(MAll/1e15) max(MAll/1e15)],'-k','linewidth',2);
% $$$ % $$$ % $$$ text(0.5,mean([min(MAll/1e15) max(MAll/1e15)]),sprintf('%3.2f ',max(MAll/1e15)-min(MAll/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ % $$$ plot([1 2],[min(MEq/1e15) min(MEq/1e15)],'--k','linewidth',2);
% $$$ % $$$ % $$$ plot([1 2],[max(MEq/1e15) max(MEq/1e15)],'--k','linewidth',2);
% $$$ % $$$ % $$$ text(1.5,mean([min(MEq/1e15) max(MEq/1e15)]),sprintf('%3.2f ',max(MEq/1e15)-min(MEq/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ % $$$ plot([2 3],[min(MEEP/1e15) min(MEEP/1e15)],'--r','linewidth',2);
% $$$ % $$$ % $$$ plot([2 3],[max(MEEP/1e15) max(MEEP/1e15)],'--r','linewidth',2);
% $$$ % $$$ % $$$ text(2.5,mean([min(MEEP/1e15) max(MEEP/1e15)]),sprintf('%3.2f ',max(MEEP/1e15)-min(MEEP/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','r','FontSize',25);
% $$$ % $$$ % $$$ plot([3 4],[min(MN/1e15) min(MN/1e15)],':b','linewidth',2);
% $$$ % $$$ % $$$ plot([3 4],[max(MN/1e15) max(MN/1e15)],':b','linewidth',2);
% $$$ % $$$ % $$$ text(3.5,mean([min(MN/1e15) max(MN/1e15)]),sprintf('%3.2f ',max(MN/1e15)-min(MN/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','b','FontSize',25);
% $$$ % $$$ % $$$ plot([4 5],[min(MS/1e15) min(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$ % $$$ % $$$                     0.5 0]);
% $$$ % $$$ % $$$ plot([4 5],[max(MS/1e15) max(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$ % $$$ % $$$                     0.5 0]);
% $$$ % $$$ % $$$ text(4.5,mean([min(MS/1e15) max(MS/1e15)]),sprintf('%3.2f ',max(MS/1e15)- min(MS/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color',[0 0.5 0],'FontSize',25);
% $$$ % $$$ % $$$ plot([5 6],[min((MN+MS)/1e15) min((MN+MS)/1e15)],':k','linewidth',2);
% $$$ % $$$ % $$$ plot([5 6],[max((MN+MS)/1e15) max((MN+MS)/1e15)],':k','linewidth',2);
% $$$ % $$$ % $$$ text(5.5,mean([min((MN+MS)/1e15) max((MN+MS)/1e15)]),sprintf('%3.2f ',max((MN+MS)/1e15)-min((MN+MS)/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ % $$$ xlim([0 6]);
% $$$ % $$$ % $$$ ylim([-0.8 0]);
% $$$ % $$$ % $$$ set(gca,'xtick',[]);
% $$$ % $$$ % $$$ set(gca,'ytick',[]);
% $$$ % $$$ % $$$ box off;
% $$$ % $$$ % $$$ grid off;
% $$$ % $$$ % $$$ axis off;
% $$$ % $$$ % $$$ title('Seasonal Range','FontSize',25);
% $$$ % $$$ 
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
% $$$     if (length(FlM(:,1,1))>3000)
% $$$         xvec = 1:10:xL;
% $$$         yvec = 1:10:yL;
% $$$     elseif (length(FlM(:,1,1))>1000)
% $$$         xvec = 1:4:xL;
% $$$         yvec = 1:4:yL;
% $$$     else
        xvec = 1:1:xL;
        yvec = 1:1:yL;
% $$$         xvec = 1:1:xL;
% $$$         yvec = 1:1:yL;
% $$$     end        
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

    %Colormap and continents:
    sp = 5;
    clim = [-125 0];

    cCH = 1; % 0 = symmetric redblue
             % 1 = negative definite parula
             % 2 = negative parula with +ve's possible
    if (cCH==0)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = redblue(npts-3);
        for i=1:(npts-3)
            if (cmap(i,:) == 1.0)
                cmap(i,:) = [0.94 0.94 0.94];
            end
        end
    elseif (cCH==1)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = parula(npts-3);
        cmap(end,:) = [0.97 0.97 0.8];
        cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
    else
        buf = 1;
        clim = [clim(1) (buf+1)*sp];
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = parula(npts-3+(buf+1));
        cmap(end+2,:) = [1 0.7 0.9];
        cmap(end-(buf+1):end-1,:) = repmat([0.97 0.97 0.8],[buf+1 1]);
        cmap(end-(buf+2),:) = (cmap(end-(buf+2),:)+cmap(end-(buf+1),:))/2;
% $$$         cmap(end-1,:) = (cmap(end,:)+cmap(end-1,:))/2;
    end        

    tmp = LAND;
    tmp(isnan(LAND)) = clim(1)-sp/2;
    tmp(~isnan(LAND)) = NaN;
    LAND = tmp;
    cmap(2:(end+1),:) = cmap;
    cmap(1,:) = [0 0 0];

    climn = [clim(1)-sp clim(2)];

    doSEAS = 0; %Plot 6 plots 2-months apart for seasonality

    if (doSEAS)
        months = {1,3,5,7,9,11};
        labels = {'(a) Jan', ...
                  '(b) Mar', ...
                  '(c) May', ...
                  '(d) Jul', ...
                  '(e) Sep', ...
                  '(f) Nov'};
        months = {2,4,6,8,10,12};
        labels = {'(a) Feb', ...
                  '(b) Apr', ...
                  '(c) Jun', ...
                  '(d) Aug', ...
                  '(e) Oct', ...
                  '(f) Dec'};
    end
    
%Mean of all months:
figure;
set(gcf,'Position',[3          59        1916         914]);
set(gcf,'defaulttextfontsize',20);
set(gcf,'defaultaxesfontsize',20);

    if (doSEAS)
        poss = [0.0867    0.6910    0.4    0.26; ...
                0.5221    0.6910    0.4    0.26; ...
                0.0867    0.3926    0.4    0.26; ...
                0.5221    0.3926    0.4    0.26; ...
                0.0867    0.0899    0.4    0.26; ...
                0.5221    0.0899    0.4    0.26];
    else
        poss = [0.1300    0.4553    0.7693    0.4697; ...
                0.1300    0.1389    0.2343    0.2680; ...
                0.3951    0.1389    0.2343    0.2680; ...
                0.6681    0.1389    0.2343    0.2680];
    end
    i = 1;
% $$$ for i=1:length(months)
% $$$     if (doSEAS)
% $$$         subplot(3,2,i);
% $$$     else        
% $$$         if (i == 1)
            subplot(5,3,[1 9]);
% $$$         else
% $$$             subplot(5,3,[10 13]+(i-2));
% $$$         end
% $$$     end
% $$$     subplot_tight(2,1,rr);
% $$$     subplot_tight(2,1,1);
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
    
    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z,cpts,'linestyle','none');
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
        ylim(cb,clim);
    end
    hold on;
    % Plot regions:
% $$$     if (i==1)
    if (exist('regions'))
        xlims = get(gca,'xlim');
        for ii=2:length(regions)
            contour(lon,lat,regions{ii}{2},[0.5 0.5],'--','color',regions{ii}{3},'linewidth',1);
% $$$             contour(lon,lat,regions{ii}{2},[0.5 0.5],'--k','linewidth',1);
            % Add text label with total:
            flon = min(lon(regions{ii}{2}));
            flat = nanmean(lat(regions{ii}{2}));
            if (ii==2)
                flat = 7;
            elseif (flat>15)
                flat = 42;
            end
            text(flon,flat,sprintf('%3.2fPW (%3.0f%%)',monmean(Field(:,ii),1,ndays(1:12))/1e15,monmean(Field(:,ii),1,ndays(1:12))/monmean(Field(:,1),1,ndays(1:12))*100),'color',regions{ii}{3},'Interpreter','none','BackgroundColor','w','Margin',0.01);
        end
    end
% $$$     end
% $$$     if (i>1)
        xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (i<=2)
        ylabel('Latitude ($^\circ$N)');
% $$$     end
    
% $$$     if (i>1)
% $$$         text(77,37.5,labels{i},'BackgroundColor','w','Margin',0.5,'HorizontalAlignment','right');
% $$$         set(gca,'xtick',[-240:60:60]);
% $$$     else
% $$$         text(-279,41,labels{i},'BackgroundColor','w','Margin',0.5);
        set(gca,'xtick',[-270:30:60]);
% $$$     end        
% $$$     set(gca,'Position',[poss(i,:)]);
    ylim([-45 45]);
% $$$     ylim([-55 55]);
% $$$     set(gca,'ytick',[-45:15:45]);
% $$$     ylim([-75 75]);
    set(gca,'ytick',[-75:15:75]);
% $$$     set(gca,'xticklabel',[]);
    set(gca,'FontSize',20);
% $$$ end 
    colormap(cmap);
    title([strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55','') ...
           ' ' num2str(Tl) '$^\circ$C Vertical Mixing']);
end
%colormap(parula);%flipud(lbmap(50,'RedBlue')));

% $$$ % $$$ %%% Plot spatial pattern of net heat flux and SST:
% $$$ % $$$ 
% $$$ % $$$ % Load Variable and calculate mean:
% $$$ % $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ % $$$ shfluxa = shflux;
% $$$ % $$$ SSTa = SST;
% $$$ % $$$ for i=2:length(outputs)
% $$$ % $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$ % $$$     shfluxa = shfluxa+shflux;
% $$$ % $$$     SSTa = SSTa+SST;
% $$$ % $$$ end
% $$$ % $$$ shflux = shfluxa/length(outputs);
% $$$ % $$$ SST = SSTa/length(outputs);
% $$$ % $$$ 
% $$$ % $$$ %If MOM01, fix NaN's in grid:
% $$$ % $$$ if (strfind(model,'01'))
% $$$ % $$$     lon = repmat(lon(:,500),[1 yL]);
% $$$ % $$$     latv = nanmean(lat,1);
% $$$ % $$$     lat = repmat(latv,[xL 1]);
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ %Sum of all positives:
% $$$ % $$$ shfluxA = mean(shflux,3);
% $$$ % $$$ shfluxP = nansum(shfluxA(shfluxA>0).*area(shfluxA>0))/1e15
% $$$ % $$$ shfluxM = nansum(shfluxA(shfluxA<0).*area(shfluxA<0))/1e15
% $$$ % $$$ 
% $$$ % $$$ %Sum of mean flux above mean isotherm:
% $$$ % $$$ SSTA = mean(SST,3);
% $$$ % $$$ isot = 21.5;
% $$$ % $$$ shfluxP = nansum(shfluxA(SSTA>=isot).*area(SSTA>=isot))/1e15
% $$$ % $$$ shfluxM = nansum(shfluxA(SSTA<isot).*area(SSTA<isot))/1e15
% $$$ % $$$ 
% $$$ % $$$ [xL,yL] = size(lon);
% $$$ % $$$ xvec = 1:3:xL;
% $$$ % $$$ yvec = 1:3:yL;
% $$$ % $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ % $$$ 
% $$$ % $$$ months = {[1:12], ...
% $$$ % $$$                             [3], ...
% $$$ % $$$                             [7], ...
% $$$ % $$$                             [11]};
% $$$ % $$$ 
% $$$ % $$$ labels = {'(a) Annual', ...
% $$$ % $$$           '(b) March', ...
% $$$ % $$$           '(c) July', ...
% $$$ % $$$           '(d) November'};
% $$$ % $$$ 
% $$$ % $$$ clim = [-200 200];
% $$$ % $$$ sp = 20;
% $$$ % $$$ 
% $$$ % $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ % $$$ npts = length(cpts)
% $$$ % $$$ cmap = redblue(npts-3);
% $$$ % $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ % $$$ set(0,'defaulttextfontsize',20);
% $$$ % $$$ set(0,'defaultaxesfontsize',20);
% $$$ % $$$ poss = [0.1300    0.415  0.7693    0.56; ...
% $$$ % $$$         0.1300    0.1    0.2343    0.2680; ...
% $$$ % $$$         0.3951    0.1    0.2343    0.2680; ...
% $$$ % $$$         0.6681    0.1    0.2343    0.2680];
% $$$ % $$$ for i=1:length(months)
% $$$ % $$$     if (i == 1)
% $$$ % $$$         subplot(5,3,[1 9]);
% $$$ % $$$     else
% $$$ % $$$         subplot(5,3,[10 13]+(i-2));
% $$$ % $$$     end
% $$$ % $$$     X = lon(xvec,yvec);
% $$$ % $$$     Y = lat(xvec,yvec);
% $$$ % $$$     Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$     Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$     Z = Z(xvec,yvec);
% $$$ % $$$     Z2 = Z2(xvec,yvec);
% $$$ % $$$     Z(Z==0) = NaN;
% $$$ % $$$     
% $$$ % $$$     % Map projection:
% $$$ % $$$ % $$$     map = 1;
% $$$ % $$$ % $$$     ax = worldmap('World');
% $$$ % $$$ % $$$     setm(ax, 'Origin', [0 260 0])
% $$$ % $$$ % $$$     contourfm(Y,X,Z,[-1e10 -500:20:500 1e10],'linestyle','none');
% $$$ % $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$ % $$$     hold on;
% $$$ % $$$ % $$$     [c,h] = contourm(Y,X,Z2,[-3:2:21 25:2:35],'-k');
% $$$ % $$$     [c,h] = contour(X,Y,Z2,[-3:2:35],'-k');
% $$$ % $$$     clabel(c,h);        
% $$$ % $$$ % $$$     hand = clabelm(c,h);        
% $$$ % $$$ % $$$     set(hand,'FontSize',10,'BackgroundColor','none');
% $$$ % $$$     if (i == 1)
% $$$ % $$$         [c,h] = contour(X,Y,Z2,[21.5 21.5],'-k','linewidth',2);
% $$$ % $$$         clabel(c,h)
% $$$ % $$$ % $$$         [c,h] = contourm(Y,X,Z2,[23 23],'-k','linewidth',2);
% $$$ % $$$ % $$$         hand = clabelm(c,h);
% $$$ % $$$ % $$$         set(hand,'FontSize',10,'BackgroundColor','none');
% $$$ % $$$     end
% $$$ % $$$     caxis(clim);
% $$$ % $$$     if (i==1)
% $$$ % $$$         cb = colorbar;
% $$$ % $$$         ylabel(cb,'Wm$^{-2}$');
% $$$ % $$$     end
% $$$ % $$$     ylim([-75 75]);
% $$$ % $$$     if (i>1)
% $$$ % $$$         xlabel('Longitude ($^\circ$E)');
% $$$ % $$$     end
% $$$ % $$$     if (i<=2)
% $$$ % $$$         ylabel('Latitude ($^\circ$N)');
% $$$ % $$$     end
% $$$ % $$$     if (i>1)
% $$$ % $$$ % $$$         text(-277,53,labels{i},'BackgroundColor','w');
% $$$ % $$$         text(-277,64,labels{i},'BackgroundColor','w','margin',0.01);
% $$$ % $$$         set(gca,'xtick',[-240:60:60]);
% $$$ % $$$     else
% $$$ % $$$ % $$$         text(-279,55,labels{i},'BackgroundColor','w');
% $$$ % $$$         text(-279,69,labels{i},'BackgroundColor','w','margin',0.01);
% $$$ % $$$         set(gca,'xtick',[-270:30:60]);
% $$$ % $$$     end        
% $$$ % $$$     set(gca,'ytick',[-60:30:60]);
% $$$ % $$$     set(gca,'Position',[poss(i,:)]);
% $$$ % $$$     set(gca,'color','k');
% $$$ % $$$ 
% $$$ % $$$ % $$$     setm(gca, 'MlabelParallel', 'south');
% $$$ % $$$ % $$$     setm(gca, 'MlabelParallel', 'south');   
% $$$ % $$$ % $$$     land = shaperead('landareas.shp', 'UseGeoCoords', true);
% $$$ % $$$ % $$$     geoshow(land, 'FaceColor', [0 0 0])
% $$$ % $$$ end 
% $$$ % $$$ colormap(cmap);
% $$$ 
% $$$ %%%%%%%%%%%% Supp FIGURES
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ 
% $$$ % $$$ %% Plot SST difference and shflux difference plots:
% $$$ % $$$ 
% $$$ % $$$ % $$$ % Load MOM01 data:
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/';
% $$$ % $$$ % $$$ model = 'MOM01';
% $$$ % $$$ % $$$ outputs = [111 222];
% $$$ % $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ % $$$ ndays = diff(time_snap);
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ [xL,yL] = size(lon)
% $$$ % $$$ % $$$ %If MOM01, fix NaN's in grid:
% $$$ % $$$ % $$$ if (strfind(model,'01'))
% $$$ % $$$ % $$$     lonv = lon(:,500);
% $$$ % $$$ % $$$     latv = nanmean(lat,1);
% $$$ % $$$ % $$$     [tmp ind] = min(abs(latv - 60))
% $$$ % $$$ % $$$     [lon,lat] = ndgrid(lonv,latv(1:ind));
% $$$ % $$$ % $$$ end
% $$$ % $$$ % $$$ MOM01lon = lon;
% $$$ % $$$ % $$$ MOM01lat = lat;
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ load([base model sprintf('_output%03d_SurfaceVars.mat', ...
% $$$ % $$$ % $$$                          outputs(1))]);
% $$$ % $$$ % $$$ MOM01shflux = shflux(:,1:ind,:);
% $$$ % $$$ % $$$ MOM01SST = SST(:,1:ind,:);
% $$$ % $$$ 
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ % $$$ % $$$ model = 'MOM025_kb3seg';
% $$$ % $$$ % $$$ outputs = [75:79];
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ % $$$ % $$$ model = 'MOM025_kb1em5';
% $$$ % $$$ % $$$ outputs = 94;
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ % $$$ model = 'MOM025';
% $$$ % $$$ % $$$ outputs = [8:12];
% $$$ % $$$ % $$$ outputs = [14]
% $$$ % $$$ % $$$ 
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em6/mat_data/';
% $$$ % $$$ % $$$ model = 'MOM025_kb1em6';
% $$$ % $$$ % $$$ outputs = 30;
% $$$ % $$$ 
% $$$ % $$$ % $$$ % Load ACCESS-OM2:
% $$$ % $$$ % $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ % $$$ % $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ % $$$ % $$$ outputs = 36;
% $$$ % $$$ 
% $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ ndays = diff(time_snap);
% $$$ % $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ % $$$ 
% $$$ % $$$ % Average all years:
% $$$ % $$$ if (length(yrs)>1)
% $$$ % $$$     shflux = reshape(shflux,[length(shflux(:,1,1)) length(shflux(1,:,1)) 12 nyrs]);
% $$$ % $$$     shflux = mean(shflux(:,:,:,yrs),4);
% $$$ % $$$     SST = reshape(SST,[length(SST(:,1,1)) length(SST(1,:,1)) 12 nyrs]);
% $$$ % $$$     SST = mean(SST(:,:,:,yrs),4);
% $$$ % $$$     ndays = ndays(1:12);
% $$$ % $$$ else
% $$$ % $$$     shfluxa = shflux;
% $$$ % $$$     SSTa = SST;
% $$$ % $$$     SSTa = SSTa(:,:,13:24);
% $$$ % $$$     ndays = ndays(1:12);
% $$$ % $$$     for i=2:length(outputs)
% $$$ % $$$         load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$ % $$$         shfluxa = shfluxa+shflux;
% $$$ % $$$         SSTa = SSTa+SST;
% $$$ % $$$     end
% $$$ % $$$     shflux = shfluxa/length(outputs);
% $$$ % $$$     SST = SSTa/length(outputs);
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ if (max(max(max(SST)))>100)
% $$$ % $$$     SST = SST-273.15;
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ % Calculate bias from first run:
% $$$ % $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ % $$$ ndays = diff(time_snap);
% $$$ % $$$ if (rr >=2)
% $$$ % $$$     SSTcur = SST;
% $$$ % $$$     load([base RUNS{1}{1} sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ % $$$     SST = reshape(SST,[length(SST(:,1,1)) length(SST(1,:,1)) 12 nyrs]);
% $$$ % $$$     SST = mean(SST(:,:,:,yrs),4);
% $$$ % $$$     if (max(max(max(SST)))>100)
% $$$ % $$$     SST = SST-273.15;
% $$$ % $$$     end
% $$$ % $$$     SSTbias = monmean(SSTcur,3,ndays) - monmean(SST,3,ndays);
% $$$ % $$$ else
% $$$ % $$$ % WOA13 SST:
% $$$ % $$$ WOAname = '/srv/ccrc/data03/z3500785/WOA13/woa13_decav_t00_04v2.nc';
% $$$ % $$$ WOASST = ncread(WOAname,'t_an',[1 1 1 1],[1440 720 1 1]);
% $$$ % $$$ [WOAlon,WOAlat] = ndgrid(ncread(WOAname,'lon'),ncread(WOAname,'lat'));
% $$$ % $$$ 
% $$$ % $$$ %Shift longitudes:
% $$$ % $$$ [tmp ind] = min(abs(WOAlon(:,1)-80));
% $$$ % $$$ WOASST = cat(1,WOASST(ind+1:end,:),WOASST(1:ind,:));
% $$$ % $$$ WOAlon = cat(1,WOAlon(ind+1:end,:)-360,WOAlon(1:ind,:));
% $$$ % $$$ WOAlat = cat(1,WOAlat(ind+1:end,:),WOAlat(1:ind,:));
% $$$ % $$$ 
% $$$ % $$$ % Calculate bias from WOA:
% $$$ % $$$ SSTbias = monmean(SST,3,ndays)-interp2(WOAlon',WOAlat',WOASST',lon,lat,'linear');
% $$$ % $$$ % $$$ for i=1:12
% $$$ % $$$ % $$$     SST(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01SST(:,:,i)',lon,lat,'linear')-SST(:,:,i);
% $$$ % $$$ % $$$     shflux(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01shflux(:,:,i)',lon,lat,'linear')-shflux(:,:,i);
% $$$ % $$$ % $$$ end
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ [xL,yL] = size(lon);
% $$$ % $$$ xvec = 1:2:xL;
% $$$ % $$$ yvec = 1:2:yL;
% $$$ % $$$ 
% $$$ % $$$ % $$$ figure;
% $$$ % $$$ % $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ % $$$ % $$$ set(gcf,'defaulttextfontsize',15);
% $$$ % $$$ % $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ % $$$ % $$$ poss = [0.15  0.56 0.3 0.34; ...
% $$$ % $$$ % $$$         0.48  0.56 0.3 0.34; ...
% $$$ % $$$ % $$$         0.15  0.1 0.3 0.34; ...
% $$$ % $$$ % $$$         0.48  0.1 0.3 0.34];
% $$$ % $$$ clvls = [-1e10 -1:0.05:1 1e10];
% $$$ % $$$ subplot(2,3,rr);
% $$$ % $$$ contourf(lon(xvec,yvec),lat(xvec,yvec),SSTbias(xvec,yvec),clvls,'linestyle', ...
% $$$ % $$$          'none');
% $$$ % $$$ hold on;
% $$$ % $$$ % $$$ contour(lon,lat,SSTbias,[-3 -2 -1 1 2 3],'-k');
% $$$ % $$$ if (rr == 1)
% $$$ % $$$     contour(lon,lat,SSTbias,[-1:0.5:1],'-k');
% $$$ % $$$ end
% $$$ % $$$ set(gca,'color','k');
% $$$ % $$$ % $$$ title('$\kappa_B=10^{-5}$ - WOA13 SST Year 2 ($^\circ$C)');
% $$$ % $$$ if (rr == 1)
% $$$ % $$$     title([strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''),' may','') ' - WOA SST ($^\circ$C)']);
% $$$ % $$$ else
% $$$ % $$$     title([strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''),' may','') ' - kds50 SST ($^\circ$C)']);
% $$$ % $$$ end    
% $$$ % $$$ if (rr == 1)
% $$$ % $$$     caxis([-1 1]);
% $$$ % $$$ else
% $$$ % $$$     caxis([-0.5 0.5]);
% $$$ % $$$ end
% $$$ % $$$ colorbar;
% $$$ % $$$ colormap(redblue(length(clvls)-1));
% $$$ % $$$ set(gca,'FontSize',15);
% $$$ % $$$ xlabel('Longitude ($^\circ$E)');
% $$$ % $$$ ylabel('Latitude ($^\circ$N)');
% $$$ % $$$ % $$$ set(gca,'Position',poss(rr,:));
% $$$ % $$$ end
% $$$ 
% $$$ % $$$ 
% $$$ % $$$ %
% $$$ % $$$ for i=1:length(months)
% $$$ % $$$     if (i == 1)
% $$$ % $$$         subplot(5,3,[1 9]);
% $$$ % $$$     else
% $$$ % $$$         subplot(5,3,[10 13]+(i-2));
% $$$ % $$$     end
% $$$ % $$$     X = lon(xvec,yvec);
% $$$ % $$$     Y = lat(xvec,yvec);
% $$$ % $$$     Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$     Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$     Z = Z(xvec,yvec);
% $$$ % $$$     Z2 = Z2(xvec,yvec);
% $$$ % $$$ % $$$     contourf(X,Y,Z2,[-1e10 -5:0.25:5 1e10],'linestyle','none');
% $$$ % $$$     contourf(X,Y,Z,[-1e10 -500:10:500 1e10],'linestyle','none');
% $$$ % $$$     hold on;
% $$$ % $$$ % $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$ % $$$ % $$$     if (i==1)
% $$$ % $$$ % $$$         [c,h] = contour(X,Y,Z,[-100:10:-10],'--k');
% $$$ % $$$ % $$$         [c,h] = contour(X,Y,Z,[10:10:100],'-k');
% $$$ % $$$ % $$$     end
% $$$ % $$$ % $$$     else
% $$$ % $$$ % $$$         [c,h] = contour(X,Y,Z2,[-3:4:35],'-k');
% $$$ % $$$ % $$$     end
% $$$ % $$$ % $$$     clabel(c,h);
% $$$ % $$$     caxis([-100 100]);
% $$$ % $$$     if (i==1)
% $$$ % $$$         cb = colorbar;
% $$$ % $$$ % $$$         ylabel(cb,'$^\circ$C');
% $$$ % $$$         ylabel(cb,'Wm$^{-2}$');
% $$$ % $$$     end
% $$$ % $$$     ylim([-75 60]);
% $$$ % $$$     if (i>1)
% $$$ % $$$         xlabel('Longitude ($^\circ$E)');
% $$$ % $$$     end
% $$$ % $$$     if (i<=2)
% $$$ % $$$         ylabel('Latitude ($^\circ$N)');
% $$$ % $$$     end
% $$$ % $$$     if (i>1)
% $$$ % $$$         text(-276,53,labels{i},'BackgroundColor','w');
% $$$ % $$$     else
% $$$ % $$$         text(-278,55,labels{i},'BackgroundColor','w');
% $$$ % $$$     end        
% $$$ % $$$     set(gca,'Position',[poss(i,:)]);
% $$$ % $$$     set(gca,'color','k');
% $$$ % $$$ end 
% $$$ % $$$ colormap(redblue);
% $$$ % $$$ 
% $$$ % $$$ %%% Outcrop area plot:
% $$$ % $$$ 
% $$$ % $$$ Aout = zeros(size(Te));
% $$$ % $$$ for i=1:length(Te)
% $$$ % $$$     i
% $$$ % $$$     Aout(i) = nansum(area(mean(SST,3)>Te(i)));
% $$$ % $$$ end
% $$$ % $$$ figure;
% $$$ % $$$ plot(Aout/1e6,Te,'-k','linewidth',2);
% $$$ % $$$ xlabel('Outcrop area (km$^2$)');
% $$$ % $$$ ylabel('Temperature ($^\circ$C)');
% $$$ % $$$ plot(-diff(Aout,[],1)/dT/1e6,T,'-k','linewidth',2);
% $$$ % $$$ xlabel('Outcrop area (km$^2$/$^\circ$C)');
% $$$ % $$$ ylabel('Temperature ($^\circ$C)');
% $$$ 

% $$$ %%% Equatorial Temperature Bias Plot  -------------------------------------------------------------------------
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(1))]);
% $$$ vars = {'temp','mld','vdif','vnlc','ndif'};
% $$$ for i=1:length(vars)
% $$$     eval(['sz = size(' vars{i} ');']);
% $$$     sz(end) = 12;
% $$$     eval([vars{i} 'all = reshape(' vars{i} ',[sz nyrs]);']);
% $$$ end
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf(['_output%03d_varsat_Eq.mat'],outputs(i))]);
% $$$     for i=1:length(vars)
% $$$         eval(['sz = size(' vars{i} ');']);
% $$$         sz(end) = 12;
% $$$         eval([vars{i} 'all = cat(4,' vars{i} 'all,reshape(' vars{i} ',[sz nyrs]));']);
% $$$     end
% $$$ end
% $$$ for i=1:length(vars)
% $$$     eval(['sz = size(' vars{i} 'all);']);
% $$$     eval([vars{i} ' = mean(' vars{i} 'all,length(sz));']);
% $$$     eval(['clear ' vars{i} 'all;']);
% $$$ end
% $$$ 
% $$$ [xL,zL,tL] = size(temp);
% $$$ TL = length(T);
% $$$ 
% $$$ months = [1:12];
% $$$ temp = monmean(temp(:,:,months),3,ndays(months));
% $$$ 
% $$$ %ACCESS-OM2:
% $$$ temp = temp-273.15;
% $$$ 
% $$$ % WOA13:
% $$$ WOAname = '/srv/ccrc/data03/z3500785/WOA13/woa13_decav_t00_04v2.nc';
% $$$ WOAlat = ncread(WOAname,'lat');
% $$$ WOAlon = ncread(WOAname,'lon');
% $$$ WOAdep = ncread(WOAname,'depth');
% $$$ [tmp Eqind] = min(abs(WOAlat));
% $$$ 
% $$$ WOAT = squeeze(ncread(WOAname,'t_an',[1 Eqind 1 1],[1440 1 102 1]));
% $$$ 
% $$$ %Shift longitudes:
% $$$ [tmp ind] = min(abs(WOAlon-80));
% $$$ WOAT = cat(1,WOAT(ind+1:end,:),WOAT(1:ind,:));
% $$$ WOAlon = cat(1,WOAlon(ind+1:end)-360,WOAlon(1:ind));
% $$$ [WOAlon,WOAdep] = ndgrid(WOAlon,WOAdep);
% $$$ 
% $$$ % Calculate bias from WOA:
% $$$ Tbias = temp-interp2(WOAlon',-WOAdep',WOAT',Xt,-Zt,'linear');
% $$$ 
% $$$ %Colormap:
% $$$ clim = [-3 3];
% $$$ sp = 0.25;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$     
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ subplot(2,3,rr);
% $$$ contourf(Xt,-Zt,Tbias,cpts,'linestyle','none');
% $$$ hold on;
% $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[0:2:35],'-k');
% $$$ clabel(c,h,[0:2:35]);
% $$$ [c,h] = contour(WOAlon,-WOAdep,WOAT,[20 20],'-k','linewidth',2);
% $$$ hold on;
% $$$ [c,h] = contour(Xt,-Zt,temp,[20 20],'--k','linewidth',2);
% $$$ ylim([-250 0]);
% $$$ xlim([-220 -80]);
% $$$ if (rr == 3)
% $$$     cb = colorbar;
% $$$ end
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Depth (m)');
% $$$ caxis(clim);
% $$$ colormap(cmap);
% $$$ set(gca,'FontSize',15);
% $$$ title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ' ...
% $$$                     ,'AOM'),'deg jra55',''),' may',''),'ryf8485 ','') ' - WOA13 Equatorial T ($^\circ$C)']);
% $$$ 
% $$$ end

% $$$ 
% $$$ %%% M-L T-t plots: 
% $$$ months = 1:12;
% $$$ fields = {
% $$$ % $$$           {F(:,months,yrs)+P(:,months,yrs), 'Surface Forcing $\mathcal{F}+\mathcal{P}$','k',2,'-'}, ...
% $$$           {M(:,months,yrs), 'Vertical Mixing $\mathcal{M}$','r',2,'-'}, ...
% $$$           {I(:,months,yrs), 'Implicit Mixing $\mathcal{I}$','b',2,'-'}, ...
% $$$ % $$$           {N(:,months,yrs), 'Tendency $\mathcal{N}$','m',2,'-'}, ...
% $$$ % $$$           {M(:,months,yrs)+I(:,months,yrs), 'Total Mixing $\mathcal{M}+\mathcal{I}$',[0 0.5 0],2,'--'}, ...
% $$$ % $$$           {P(:,months,yrs), 'P-E+R Effective Heat Flux $\mathcal{P}$',0.5*[1 1 1],2,':'}, ...
% $$$ % $$$           {S(:,months,yrs), 'Submesoscale',0.5*[0 0 1],1,'--'}, ...
% $$$ % $$$           {F(:,months,yrs)+P(:,months,yrs), '$\mathcal{F}+\mathcal{P}$','k',2,'--'}, ...
% $$$ % $$$           {Nmon(:,months,yrs), 'Monthly-Binned Tendency','m',2,'--'}, ...
% $$$ % $$$           {Nsnap(:,months,yrs), 'Calculated Tendency',0.3*[1 1 1],2,':'}, ...
% $$$ % $$$          {N2(:,months,yrs), 'Calculated Tendency 2',0.3*[1 1 1],2,'--'}, ...
% $$$           };
% $$$ 
% $$$ % Fluxes:
% $$$ scale = 1/1e15;label = '(PW)';
% $$$ caxs = [-15 15];x = Te;
% $$$ 
% $$$ figure;
% $$$ %set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'Position',[272         107        1302         863]);
% $$$ for ii=1:length(fields)
% $$$     subplot(1,length(fields),ii);
% $$$     V = mean(fields{ii}{1},3)'*scale;
% $$$     [X,Y] = ndgrid(1:tL,Te);
% $$$     cint = [-1e10 -10:0.05:10 1e10];
% $$$     contourf(X,Y,V,cint);%,'linestyle','none');
% $$$     hold on;
% $$$     plot([1 tL],[22.5 22.5],'--k','linewidth',2);
% $$$     cb = colorbar('Location','NorthOutside','FontSize',25);    
% $$$ % $$$     set(gca,'xtick',0.5:1:11.5);
% $$$ % $$$     set(gca,'xticklabel',[1:12]);
% $$$     set(gca,'ytick',-5:5:35);
% $$$     set(gca,'xtick',[1:tL]);
% $$$     ylim([-3 31]);
% $$$     grid on;
% $$$     %    caxis(caxs);
% $$$     caxis([-1 0]);
% $$$     xlabel('Month');
% $$$     ylabel('Temperature ($^\circ$C)');
% $$$     %    xlabel(cb,['MOM025 Global ' fields{ii}{2} ' '
% $$$     %    label],'FontSize',20);
% $$$     xlabel(cb,[fields{ii}{2} ' (PW)'],'FontSize',25);
% $$$     set(gca,'FontSize',25);
% $$$ end
% $$$ % $$$ colormap(redblue);
% $$$ % $$$ caxis([-2 2]);
% $$$ cmap = redblue((length(cint)-3)*2);
% $$$ cmap = cmap(1:(length(cint)-3),:);
% $$$ colormap(cmap);
% $$$ % $$$ caxis([-5 5]);
% $$$ hold on;
% $$$ plot([1 tL],[22.5 22.5],'--k','linewidth',2);
% $$$ 
% $$$ %%% Spatial Integrated Mixing
% $$$ 
% $$$ % Load Base Variables:
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = [8:12];
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ %Region choice:
% $$$ LATsplit = 10;
% $$$ reg = [-150 -90 -5 5];
% $$$ EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);
% $$$ 
% $$$ obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ LAND = obj.SST(:,:,1);
% $$$ 
% $$$ %If MOM01, fix NaN's in grid:
% $$$ if (strfind(model,'01'))
% $$$     lon = repmat(lon(:,500),[1 yL]);
% $$$     latv = nanmean(lat,1);
% $$$     lat = repmat(latv,[xL 1]);
% $$$ end
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:1:xL;
% $$$ yvec = 1:1:yL;
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ 
% $$$ months = [1:12];
% $$$ 
% $$$ %Colormap and continents:
% $$$ clim = [-50 0];
% $$$ sp = 2.5;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ 
% $$$ % $$$ cmap = flipud(lbmap(2*(npts-3),'RedBlue'));
% $$$ % $$$ cmap = cmap(1:(npts-3),:);
% $$$ % $$$ cmap = redblue(2*(npts-3));
% $$$ % $$$ cmap = cmap(1:(npts-3),:);
% $$$ %cmap = summer(npts-3);
% $$$ cmap = parula(npts-3);
% $$$ 
% $$$ %Custom parula:
% $$$ cmap = parula(npts-3);
% $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ % $$$ mid = round((npts-3)/2);
% $$$ % $$$ cmap(2:(mid-1),1) = interp1([1 mid],[cmap(1,1) cmap(mid,1)],2:(mid-1),'linear');
% $$$ % $$$ cmap(2:(mid-1),2) = interp1([1 mid],[cmap(1,2) cmap(mid,2)],2:(mid-1),'linear');
% $$$ % $$$ cmap(2:(mid-1),3) = interp1([1 mid],[cmap(1,3) cmap(mid,3)],2:(mid-1),'linear');
% $$$ % $$$ cmap((mid+1):(end-1),1) = interp1([mid npts-3],[cmap(mid,1) cmap(end,1)],(mid+1):(npts-4),'linear');
% $$$ % $$$ cmap((mid+1):(end-1),2) = interp1([mid npts-3],[cmap(mid,2) cmap(end,2)],(mid+1):(npts-4),'linear');
% $$$ % $$$ cmap((mid+1):(end-1),3) = interp1([mid npts-3],[cmap(mid,3) cmap(end,3)],(mid+1):(npts-4),'linear');
% $$$ 
% $$$ 
% $$$ % $$$ %Custom:
% $$$ % $$$ cmap = zeros(npts-3,3);
% $$$ % $$$ cmap(1,:) = [0 0 1];
% $$$ % $$$ cmap(end,:) = [0.9359 0.9323 0.7947];
% $$$ % $$$ %cmap(end,:) = [1 1 0.85];
% $$$ % $$$ %cmap(end,:) = [0.97 0.97 0.8];
% $$$ % $$$ cmap(2:(end-1),1) = interp1([1 npts-3],[cmap(1,1) cmap(end,1)],2:(npts-4),'linear');
% $$$ % $$$ cmap(2:(end-1),2) = interp1([1 npts-3],[cmap(1,2) cmap(end,2)],2:(npts-4),'linear');
% $$$ % $$$ cmap(2:(end-1),3) = interp1([1 npts-3],[cmap(1,3) cmap(end,3)],2:(npts-4),'linear');
% $$$ 
% $$$ tmp = LAND;
% $$$ tmp(isnan(LAND)) = clim(1)-sp/2;
% $$$ tmp(~isnan(LAND)) = NaN;
% $$$ LAND = tmp;
% $$$ cmap(2:(end+1),:) = cmap;
% $$$ cmap(1,:) = [0 0 0];
% $$$ 
% $$$ climn = [clim(1)-sp clim(2)];
% $$$ 
% $$$ %Mean of all months:
% $$$ figure;
% $$$ set(gcf,'Position',[3    40   956   963]);
% $$$ set(0,'defaulttextfontsize',15);
% $$$ set(0,'defaultaxesfontsize',15);
% $$$ poss = [0.100    0.6993    0.7750    0.2837; ...
% $$$         0.100    0.3996    0.7750    0.2837; ...
% $$$         0.100    0.100    0.7750    0.2837];
% $$$ 
% $$$ Tls = [27.5 17.5 10]
% $$$ labels = {'(a) 27.5$^\circ$C','(b) 17.5$^\circ$C','(c) 10$^\circ$C'};
% $$$ 
% $$$ % $$$ Tls = [-1.5 -1 2.5]
% $$$ % $$$ labels = {'-1.5$^\circ$C','-1$^\circ$C','2.5$^\circ$C'};
% $$$ 
% $$$ for i=1:length(Tls)
% $$$     Tl = Tls(i)
% $$$     % Load Variable and calculate mean:
% $$$     load([base model sprintf('_output%03d_VertInt_T',outputs(1)) strrep(num2str(Tl),'.','p') 'C.mat']);
% $$$     if (strfind(model,'025'))
% $$$         FlM(isnan(FlM)) = 0.0;
% $$$         FlMa = FlM;
% $$$         for ii=2:length(outputs)
% $$$             load([base model sprintf('_output%03d_VertInt_T',outputs(ii)) strrep(num2str(Tl),'.','p') 'C.mat']);
% $$$             FlM(isnan(FlM)) = 0.0;
% $$$             FlMa = FlMa+FlM;
% $$$         end
% $$$         FlM = FlMa/length(outputs);
% $$$         FlM(FlM==0) = NaN;
% $$$     end
% $$$ 
% $$$     subplot(3,1,i);
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     if (length(months)>1)
% $$$         tmp = FlM;
% $$$         tmp(isnan(tmp)) = 0.0;
% $$$         Z = monmean(tmp(:,:,months),3,ndays(months));
% $$$         Z(Z == 0) = NaN;
% $$$     else
% $$$         Z = FlM(:,:,months);
% $$$     end
% $$$     Z = Z(xvec,yvec);
% $$$     
% $$$     Z(Z<clim(1)) = clim(1);
% $$$ % $$$     contourf(X,Y,Z.*cos(Y/180*pi),cpts,'linestyle','none');
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$     hold on;    
% $$$     contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
% $$$     caxis(climn);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Wm$^{-2}$');
% $$$     ylim(cb,clim);
% $$$     hold on;
% $$$     xlims = get(gca,'xlim');
% $$$ % $$$     plot(xlims,LATsplit*[1 1],'--k');
% $$$ % $$$     plot(xlims,-LATsplit*[1 1],'--k');
% $$$ % $$$     plot([reg(1:2) reg(2:-1:1) reg(1)],[reg(3) reg(3:4) reg(4:-1:3)],'--k');
% $$$     ylim([-75 75]);
% $$$     if (i==3)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     ylabel('Latitude ($^\circ$N)');
% $$$     text(-276,65,labels{i},'BackgroundColor','w','Margin',0.5);
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$ % $$$     set(gca,'ytick',[-45:15:45]);
% $$$     set(gca,'xtick',[-360:60:360]);
% $$$     if (i~=3)
% $$$         set(gca,'xticklabel',[]);
% $$$     end
% $$$ end 
% $$$ colormap(cmap);
% $$$ %colormap(parula);%flipud(lbmap(50,'RedBlue')));
% $$$ 
% $$$ %%% Meridional heat flux:
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
% $$$ 
% $$$ %Calculate fractions of surface heat flux regionally:
% $$$ 
% $$$ %Region choice:
% $$$ LATsplit = 10;
% $$$ reg = [-150 -90 -5 5];
% $$$ EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);
% $$$ 
% $$$ Eqinds = abs(lat)<=LATsplit;
% $$$ Ninds = lat>LATsplit;
% $$$ Sinds = lat<-LATsplit;
% $$$ 
% $$$ for i=1:12
% $$$     SHFtmp = shflux(:,:,i);
% $$$     SHFEq(i) = nansum(nansum(area(Eqinds).*SHFtmp(Eqinds),1),2);
% $$$     SHFN(i) = nansum(nansum(area(Ninds).*SHFtmp(Ninds),1),2);
% $$$     SHFS(i) = nansum(nansum(area(Sinds).*SHFtmp(Sinds),1),2);
% $$$     SHFEEP(i) = nansum(nansum(area(EEPinds).*SHFtmp(EEPinds),1),2);
% $$$ 
% $$$     Atmp = area;
% $$$     Atmp(isnan(SHFtmp)) = NaN;
% $$$     AREAtotal(i) = nansum(nansum(Atmp));
% $$$     AREAEEP(i) = nansum(nansum(Atmp(EEPinds)));
% $$$     AREAEq(i) = nansum(nansum(Atmp(Eqinds)));
% $$$     AREAN(i) = nansum(nansum(Atmp(Ninds)));
% $$$     AREAS(i) = nansum(nansum(Atmp(Sinds)));
% $$$ end
% $$$ SHFAll = SHFEq+SHFN+SHFS;
% $$$ 
% $$$ %Display output in terminal for table:
% $$$ str = {['Area fractions model ' model] ; ...
% $$$ sprintf(' NH area = %3.2f',monmean(AREAN,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' SH area = %3.2f',monmean(AREAS,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' Eq area = %3.2f',monmean(AREAEq,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' EEP area = %3.2f',monmean(AREAEEP,2,ndays)/monmean(AREAtotal,2,ndays))}
% $$$ str = {['Annual totals model ' model]  ; ...
% $$$ sprintf(' Total = %3.2f',monmean(SHFAll,2,ndays)/1e15) ; ...
% $$$ sprintf(' NH = %3.2fPW (%3.0f)',monmean(SHFN,2,ndays)/1e15,monmean(SHFN,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' SH = %3.2fPW (%3.0f)',monmean(SHFS,2,ndays)/1e15,monmean(SHFS,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' Eq = %3.2fPW (%3.0f)',monmean(SHFEq,2,ndays)/1e15,monmean(SHFEq,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' EEP = %3.2fPW (%3.0f)',monmean(SHFEEP,2,ndays)/1e15,monmean(SHFEEP,2,ndays)/monmean(SHFAll,2,ndays)*100)}
% $$$ 
% $$$ % Plot regional heat fluxes as a function of season:
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ 
% $$$ % $$$ subplot(3,1,1);
% $$$ plot(1:12,SHFEEP/1e15,'--r','LineWidth',4);
% $$$ hold on;
% $$$ plot(1:12,SHFEq/1e15,'--k','LineWidth',4);
% $$$ plot(1:12,SHFN/1e15,':b','LineWidth',4);
% $$$ plot(1:12,SHFS/1e15,':','color',[0 0.5 0],'LineWidth',4);
% $$$ plot(1:12,(SHFN+SHFS)/1e15,':k','LineWidth',4);
% $$$ plot(1:12,SHFAll/1e15,'-k','LineWidth',4);
% $$$ xlabel('Month');
% $$$ ylabel(['PW']);
% $$$ title(['Surface Heat Flux']);
% $$$ leg = legend('Eastern Equatorial Pacific', ...
% $$$              'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
% $$$              ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
% $$$              ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
% $$$              ['Total']);
% $$$ set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
% $$$ ylim([-20 20]);
% $$$ xlim([1 12]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'xtick',[1:12]);
% $$$ set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
% $$$ set(gca,'FontSize',25);
% $$$ grid on;
% $$$ 
% $$$ 
% $$$ %%% Latitudinal Slices:
% $$$ 
% $$$ % Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/access-om2/1deg_jra55_ryf8485_kds50_s13p8_mushy/mat_data/';
% $$$ % $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_s13p8_mushy';
% $$$ % $$$ outputs = 57;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_ENSO/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = 4;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM01_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM01';
% $$$ % $$$ outputs = [111 222];
% $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ lonsl = 140;
% $$$ load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(1))]);
% $$$ vars = {'temp','u','v','kappa','taux','tauy','mld','vdif','vnlc','pmer','sufc','swrd'};
% $$$ for i=1:length(vars)
% $$$     eval([vars{i} 'a = ' vars{i} ';']);
% $$$ end
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf(['_output%03d_varsat_' num2str(lonsl) 'W.mat'],outputs(i))]);
% $$$     for i=1:length(vars)
% $$$         eval([vars{i} 'a = ' vars{i} 'a + ' vars{i} ';']);
% $$$     end
% $$$ end
% $$$ for i=1:length(vars)
% $$$     eval([vars{i} ' = ' vars{i} 'a/length(outputs);']);
% $$$     eval(['clear ' vars{i} 'a;']);
% $$$ end
% $$$ if(max(max(max(temp)))>100)
% $$$     temp = temp-273.15;
% $$$ end
% $$$ 
% $$$ [yL,zL,tL] = size(temp);
% $$$ TL = length(T);
% $$$ 
% $$$ % Depth of isotherms:
% $$$ Zi = zeros(yL,TL,tL);
% $$$ for ti=1:tL
% $$$     for yi=1:yL
% $$$         tvec = squeeze(temp(yi,:,ti));
% $$$         zvec = -Zt(yi,:);
% $$$         tvec(isnan(tvec)) = -1000;
% $$$         tvec = tvec - 0.01*(1:zL);
% $$$         Zi(yi,:,ti) = interp1(tvec,zvec,T,'linear');
% $$$         ind = find(~isnan(Zi(yi,:,ti)),1,'last');
% $$$         Zi(yi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
% $$$     end
% $$$ end
% $$$ Yi = repmat(Yt(:,1),[1 TL]);
% $$$ 
% $$$ var = cumsum(vdif,2,'reverse'); % Vertical Mixing Flux
% $$$ clim = [-250 0];
% $$$ sp = 10;
% $$$ doWMT = 0;
% $$$ 
% $$$ % $$$ var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
% $$$ % $$$ var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
% $$$ % $$$ clim = [-1 1]*1e-5*86400; % FOR WMT
% $$$ % $$$ sp = 0.1*1e-5*86400;
% $$$ % $$$ doWMT = 1;
% $$$ 
% $$$ months = {[1:12],[3],[7],[11]};
% $$$ monthsu01 = {[1:4],[1],[3],[4]};
% $$$ labels = {'Annual','March','July','November'};
% $$$ % $$$ 
% $$$ % $$$ months = {[6],[7],[8],[9]};
% $$$ % $$$ monthsu01 = {[2],[3],[3],[3]};
% $$$ % $$$ labels = {'June','July','August','September'};
% $$$ % $$$ 
% $$$ % $$$ months = {[7:9],[10:12],[8],[11]};
% $$$ % $$$ labels = {'JAS','OND','August','November'};
% $$$ 
% $$$ %Colormap:
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ 
% $$$ if (doWMT)
% $$$     cmap = redblue(npts-3);
% $$$ else
% $$$     cmap = parula(npts-3);
% $$$     cmap(end,:) = [1 1 1];
% $$$     cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ end
% $$$ 
% $$$ % $$$ %Save for schematic:
% $$$ % $$$ Xl = Yi;
% $$$ % $$$ Yl = nanmonmean(Zi(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$ Zl = nanmonmean(vdif(:,:,months{i}),3,ndays(months{i}));
% $$$ % $$$ XlC = Yt;
% $$$ % $$$ YlC = -Zt;
% $$$ % $$$ ZlC = monmean(temp(:,:,months{i}),3,ndays(months{i}));
% $$$ 
% $$$ % $$$ [tmp Eqind] = min(abs(Yu(:,1)));
% $$$ % $$$ tauweight = abs(taux(Eqind,:))*200;
% $$$ % $$$ % Wind stress vectors:
% $$$ % $$$ sp  =5;
% $$$ % $$$ yvec = Yu(:,1);
% $$$ figure;
% $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ set(gcf,'defaulttextfontsize',20);
% $$$ set(gcf,'defaultaxesfontsize',20);
% $$$ for i=1:length(months)
% $$$ subplot(2,2,i);
% $$$ contourf(Yi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
% $$$ hold on;
% $$$ [c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
% $$$ clabel(c,h,[0:2:35]);
% $$$ [c,h] = contour(Yt,-Zt,monmean(temp(:,:,months{i}),3, ...
% $$$                                ndays(months{i})),[21.5 21.5],'-k','linewidth',2);
% $$$ if (strcmp(model,'MOM01'))
% $$$     mnu = monthsu01{i};
% $$$ else
% $$$     mnu = months{i};
% $$$ end
% $$$ ucol = [0.8706    0.4902         0];
% $$$ [c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
% $$$                 'color',ucol);
% $$$ [c,h] = contour(Yu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
% $$$                 'color',ucol);
% $$$ plot(Yu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
% $$$ ylim([-200 0]);
% $$$ %ylim([-250 0]);
% $$$ xlim([-10 10]);
% $$$ cb = colorbar;
% $$$ if (doWMT)
% $$$     ylabel(cb,'m/day');
% $$$ else
% $$$     ylabel(cb,'Wm$^{-2}$');
% $$$ end
% $$$ xlabel('Latitude ($^\circ$N)');
% $$$ ylabel('Depth (m)');
% $$$ caxis(clim);
% $$$ text(-9.6,-188,labels{i},'Backgroundcolor','w','FontSize',20);
% $$$ text(9.6,-188,[num2str(lonsl) '$^\circ$W'],'Backgroundcolor','w','FontSize',20,'HorizontalAlignment','Right');
% $$$ 
% $$$ % $$$ %Add wind-stress vectors:
% $$$ % $$$ pos = get(gca,'Position')
% $$$ % $$$ wsh = axes('Position',[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]);
% $$$ % $$$ % $$$ for ii=1:sp:length(yvec)
% $$$ % $$$     plot(0,0,'o','MarkerSize',abs(mean(mean(taux(yvec>=-10 & yvec<=10,months{i}),2),1))*200);
% $$$ % $$$     hold on;
% $$$ % $$$ % $$$ end
% $$$ % $$$ %quiver(yvec,zeros(size(yvec)),mean(taux(1:sp:end,months{i}),2),mean(tauy(1:sp:end,months{i}),2));
% $$$ % $$$ xlim([-10 10]);
% $$$ % $$$ ylim([-1 1]);
% $$$ % $$$ box off;axis off;
% $$$ % $$$ set(wsh,'Position',[[pos(1) pos(2)+pos(4)+0.005 pos(3) 0.03]]);
% $$$ 
% $$$ LabelAxes(gca,i,20,0.008,0.95);
% $$$ end
% $$$ colormap(cmap);
% $$$ 

%%% Equatorial Slices:

% Load Variable and calculate mean:
load([base model sprintf(['_output%03d_varsat_EqPM5.mat'],outputs(1))]);
vars = {'temp','mld','vdif','vnlc','ndif'};
for i=1:length(vars)
    eval(['sz = size(' vars{i} ');']);
    sz(end) = 12;
    eval([vars{i} 'all = reshape(' vars{i} ',[sz nyrs]);']);
end
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_varsat_EqPM5.mat'],outputs(i))]);
    for i=1:length(vars)
        eval(['sz = size(' vars{i} ');']);
        sz(end) = 12;
        eval([vars{i} 'all = cat(4,' vars{i} 'all,reshape(' vars{i} ',[sz nyrs]));']);
    end
end
for i=1:length(vars)
    eval(['sz = size(' vars{i} 'all);']);
    eval([vars{i} ' = mean(' vars{i} 'all,length(sz));']);
    eval(['clear ' vars{i} 'all;']);
end

[xL,zL,tL] = size(temp);
TL = length(T);

% $$$ %ACCESS-OM2:
% $$$ temp = temp-273.15;

% Depth of isotherms:
Zi = zeros(xL,TL,tL);
for ti=1:tL
    for xi=1:xL
        tvec = squeeze(temp(xi,:,ti));
        zvec = -Zt(xi,:);
        tvec(isnan(tvec)) = -1000;
        tvec = tvec - 0.01*(1:zL);
        Zi(xi,:,ti) = interp1(tvec,zvec,T,'linear');
        ind = find(~isnan(Zi(xi,:,ti)),1,'last');
        Zi(xi,(ind+1):end,ti) = max(zvec);%linspace(Zi(yi,ind;
    end
end
Xi = repmat(Xt(:,1),[1 TL]);

% $$$ var = cumsum(vdif+vnlc,2,'reverse'); % Vertical Mixing Flux
var = ndif; % Vertical Mixing Flux
clim = [-35 0];
sp = 0.25;
doWMT = 0;

% $$$ var = (vdif+vnlc)/rho0/Cp*86400; % Vertical Mixing Transformation
% $$$ var = (pmer+sufc)/rho0/Cp*86400; % Surface Forcing Transformation
% $$$ clim = [-1 1]*1e-5*86400; % FOR WMT
% $$$ sp = 0.1*1e-5*86400;
% $$$ doWMT = 0;

months = {[1:12]};
% $$$ months = {[1:12],[3],[7],[11]};
% $$$ monthsu01 = {[1:4],[1],[3],[4]};
% $$$ labels = {'Annual','March','July','November'};
% $$$ 
% $$$ months = {[6],[7],[8],[9]};
% $$$ monthsu01 = {[2],[3],[3],[3]};
% $$$ labels = {'June','July','August','September'};
% $$$ 
% $$$ months = {[7:9],[10:12],[8],[11]};
% $$$ labels = {'JAS','OND','August','November'};

%Colormap:
cpts = [-1e10 clim(1):sp:clim(2) 1e10];
npts = length(cpts)

if (doWMT)
    cmap = redblue(npts-3);
else
    cmap = parula(npts-3);
    cmap(end,:) = [1 1 1];
    cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
end

% $$$ %Save for schematic:
% $$$ Xe = Xi;
% $$$ Ye = nanmonmean(Zi(:,:,months{i}),3,ndays(months{i}));
% $$$ Ze = nanmonmean(vdif(:,:,months{i}),3,ndays(months{i}));
% $$$ XeC = X;
% $$$ YeC = -Z;
% $$$ ZeC = monmean(temp(:,:,months{i}),3,ndays(months{i}));

% $$$ [tmp Eqind] = min(abs(Yu(:,1)));
% $$$ tauweight = abs(taux(Eqind,:))*200;
% $$$ % Wind stress vectors:
% $$$ sp  =5;
% $$$ yvec = Yu(:,1);

% $$$ figure;
% $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ for i=1:length(months)
i= 1;
% $$$ subplot(2,3,rr);
contourf(Xi,nanmonmean(Zi(:,:,months{i}),3,ndays(months{i})),nanmonmean(var(:,:,months{i}),3,ndays(months{i})),cpts,'linestyle','none');
hold on;
[c,h] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[0:1:35],'-k');
clabel(c,h,[0:2:35]);
[c,h] = contour(Xt,-Zt,nanmonmean(temp(:,:,months{i}),3,ndays(months{i})),[23 23],'-k','linewidth',2);
% $$$ if (strcmp(model,'MOM01'))
% $$$     mnu = monthsu01{i};
% $$$ else
% $$$     mnu = months{i};
% $$$ end
% $$$ ucol = [0.8706    0.4902         0];
% $$$ [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[-2:0.2:-0.2],'--', ...
% $$$                 'color',ucol);
% $$$ [c,h] = contour(Xu,-Zu,mean(u(:,:,mnu),3),[0.2:0.2:2],'-', ...
% $$$                 'color',ucol);
plot(Xu(:,1),-monmean(mld(:,months{i}),2,ndays(months{i})),'--','color',[0 0.5 0],'linewidth',3);
% $$$ clabel(c,h,'color','w');
ylim([-250 0]);
xlim([-200 -80]);
% $$$ if (rr == 3 | rr == 5)
cb = colorbar;
if (doWMT)
    ylabel(cb,'m/day');
else
    ylabel(cb,'Wm$^{-2}$');
end
% $$$ end
xlabel('Longitude ($^\circ$E)');
ylabel('Depth (m)');
set(gca,'FontSize',15);
caxis(clim);
% $$$ text(-218,-288,labels{i},'Backgroundcolor','w','FontSize',20);
title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ' ...
                    ,'AOM'),'deg jra55',''),' may',''),'ryf8485 ','') ...
       ' Numerical Mixing']);

% $$$ LabelAxes(gca,i,20,0.008,0.95);
% $$$ end
colormap(cmap);

end

% $$$ 
% $$$ 
% $$$ %%% Temperature-latitude heat function:
% $$$ 
% $$$ % Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ % $$$ outputs = 75;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ % $$$ model = 'MOM025_kb1em5';
% $$$ % $$$ outputs = 94;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = [8:12];
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em6/mat_data/';
% $$$ % $$$ model = 'MOM025_kb1em6';
% $$$ % $$$ outputs = 30;
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load streamfunction:
% $$$ load([base model sprintf('_output%03d_Tpsi.mat',outputs(1))]);
% $$$ NaNs = PSI == 0;
% $$$ % $$$ PSI = cumsum(-PSI/rho0*1e9,2); %cumsum from low temps and correct units
% $$$ PSI = cumsum(PSI/rho0*1e9,2,'reverse');
% $$$ 
% $$$ % Load SST:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ SST(SST==0) = NaN;
% $$$ 
% $$$ maxSST = squeeze(max(max(SST,[],1),[],3));
% $$$ meanSST = squeeze(nanmean(nanmean(SST,1),3));
% $$$ minSST = squeeze(min(min(SST,[],1),[],3));
% $$$ 
% $$$ %calculate heat function:
% $$$ H = cumsum(rho0*Cp*PSI*dT,2);
% $$$ 
% $$$ PSI(NaNs) = NaN;
% $$$ H(NaNs) = NaN;
% $$$ 
% $$$ [YY,TT] = ndgrid(latv,T);
% $$$ 
% $$$ %monthsc = {[1:12],[3],[7],[11]};
% $$$ monthsc = {[1:12]};%,[3],[7],[11]};
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[157         359        1742         586]);
% $$$ %set(gcf,'Position',[3    40   956   963]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ for i=1:length(monthsc)
% $$$     months = monthsc{i}
% $$$ subplot(length(monthsc),2,2*i-1);
% $$$ contourf(YY,TT,monmean(PSI(:,:,months)/1e6,3,ndays(months)),[-1e10 -100:5:100 1e10],'-k');
% $$$ hold on;
% $$$ contour(YY,TT,monmean(PSI(:,:,months)/1e6,3,ndays(months)),[0 0],'-k','linewidth',2);
% $$$ % $$$ plot(latv,maxSST,'--r','linewidth',2);
% $$$ plot(latv,meanSST,'--','color',[0 0.5 0],'linewidth',2);
% $$$ % $$$ plot(latv,minSST,'--b','linewidth',2);
% $$$ caxis([-30 30]);
% $$$ colormap(redblue);
% $$$ xlabel('Latitude ($^\circ$N)');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ text(-80,32,'(a) Transport Streamfunction (Sv)');
% $$$ colorbar;
% $$$ set(gca,'Position',[0.1300    0.1313    0.2963    0.7937]);
% $$$ 
% $$$ subplot(length(monthsc),2,2*i);
% $$$ contourf(YY,TT,monmean(H(:,:,months)/1e15,3,ndays(months)),[-1e10 -5:0.15:5 1e10],'-k');
% $$$ hold on;
% $$$ contour(YY,TT,monmean(H(:,:,months)/1e15,3,ndays(months)),[0 0],'-k','linewidth',2);
% $$$ % $$$ plot(latv,maxSST,'--r','linewidth',2);
% $$$ plot(latv,meanSST,'--','color',[0 0.5 0],'linewidth',2);
% $$$ % $$$ plot(latv,minSST,'--b','linewidth',2);
% $$$ caxis([-1.5 1.5]);
% $$$ colormap(redblue);
% $$$ xlabel('Latitude ($^\circ$N)');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ text(-80,32,'(b) Heat Function (PW)');
% $$$ colorbar;
% $$$ set(gca,'Position',[0.5306    0.1261    0.2963    0.7937]);
% $$$ end
% $$$ 
% $$$ % $$$ load([base model sprintf('_output%03d_HFunc.mat',outputs(1))]);
% $$$ 
% $$$ % $$$ TENa = HFETS+HFFRZ+HFKNL+HFPME+HFRMX+HFSUB+HFSWH+HFVDF+HFVDS;
% $$$ % $$$ latv = max(lat,[],1);
% $$$ 
% $$$ 
% $$$ %%% Plot Seasonal cycle of wind stress and SST
% $$$ 
% $$$ % Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = [8:12];
% $$$ % $$$ model = 'MOM01';
% $$$ % $$$ outputs = [333];
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ SSTa = SST;
% $$$ tauxa = taux;
% $$$ tauya = tauy;
% $$$ for i=2:length(outputs)
% $$$     i
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     SSTa = SSTa+SST;
% $$$     tauxa = tauxa+taux;
% $$$     tauya = tauya+tauy;
% $$$ end
% $$$ SST = SSTa/length(outputs);
% $$$ taux = tauxa/length(outputs);
% $$$ tauy = tauya/length(outputs);
% $$$ 
% $$$ ylims = [-30 30];
% $$$ xlims = [-260 -60];
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ [tmp ln1] = min(abs(lon(:,1)-xlims(1)));
% $$$ [tmp ln2] = min(abs(lon(:,1)-xlims(2)));
% $$$ [tmp lt1] = min(abs(lat(1,:)-ylims(1)));
% $$$ [tmp lt2] = min(abs(lat(1,:)-ylims(2)));
% $$$ xvec = (ln1-2):1:(ln2-2);
% $$$ yvec = (lt1-2):1:(lt2-2);
% $$$ [xL,yL] = size(lonu);
% $$$ [tmp ln1] = min(abs(lonu(:,1)-xlims(1)));
% $$$ [tmp ln2] = min(abs(lonu(:,1)-xlims(2)));
% $$$ [tmp lt1] = min(abs(latu(1,:)-ylims(1)));
% $$$ [tmp lt2] = min(abs(latu(1,:)-ylims(2)));
% $$$ xvec2 = (ln1-2):15:(ln2-2);
% $$$ yvec2 = (lt1-2):15:(lt2-2);
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ 
% $$$ months = {[1],[3],[5],[7],[9],[11]};
% $$$ % $$$           [1:12], ...
% $$$ % $$$           [3], ...
% $$$ % $$$           [7], ...
% $$$ % $$$           [11]};
% $$$ 
% $$$ labels = {'Jan', ...
% $$$           'Mar', ...
% $$$           'May', ...
% $$$           'Jul', ...
% $$$           'Sep', ...
% $$$           'Nov'};
% $$$ % $$$ 'Annual', ...
% $$$ % $$$           'March', ...
% $$$ % $$$           'July', ...
% $$$ % $$$           'November'};
% $$$ 
% $$$ % $$$ %Colormap:
% $$$ % $$$ clim = [0 0.1];
% $$$ % $$$ sp = 0.01;
% $$$ % $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ % $$$ npts = length(cpts)
% $$$ % $$$ cmap = parula(npts-3);
% $$$ % $$$ cmap = parula(npts-3);
% $$$ % $$$ % $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ % $$$ cmap(end,:) = [1 1 1];
% $$$ % $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ % $$$ cmap = flipud(cmap);
% $$$ 
% $$$ 
% $$$ clim = [20 30];
% $$$ sp = 1;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = colormap(redblue(npts-3));
% $$$ cmap = colormap(flipud(lbmap(npts-3,'RedBlue')));
% $$$ 
% $$$ figure;
% $$$ %set(gcf,'Position',[3          59        1916         914]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ set(gcf,'Position',[322          58        1247         945]);
% $$$ % $$$ poss = [0.1300    0.4553    0.7693    0.4697; ...
% $$$ % $$$         0.1300    0.1389    0.2343    0.2680; ...
% $$$ % $$$         0.3951    0.1389    0.2343    0.2680; ...
% $$$ % $$$         0.6681    0.1389    0.2343    0.2680];
% $$$ poss = [0.0867    0.6910    0.4    0.26; ...
% $$$         0.5221    0.6910    0.4    0.26; ...
% $$$         0.0867    0.3926    0.4    0.26; ...
% $$$         0.5221    0.3926    0.4    0.26; ...
% $$$         0.0867    0.0899    0.4    0.26; ...
% $$$         0.5221    0.0899    0.4    0.26];
% $$$ for i=1:length(months)
% $$$ % $$$     if (i == 1)
% $$$ % $$$         subplot(5,3,[1 9]);
% $$$ % $$$     else
% $$$ % $$$         subplot(5,3,[10 13]+(i-2));
% $$$ % $$$     end
% $$$     subplot(3,2,i);
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     Z = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$     Z2 = monmean(taux(:,:,months{i}),3,ndays(months{i}));
% $$$     Z3 = monmean(tauy(:,:,months{i}),3,ndays(months{i}));
% $$$     Z = Z(xvec,yvec);
% $$$     Z(Z==0) = NaN;
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$ % $$$     contourf(lonu(xvec,yvec),latu(xvec,yvec),sqrt(Z2(xvec,yvec).^2+ ...
% $$$ % $$$                                                   Z3(xvec,yvec).^2), ...
% $$$ % $$$              [-1e10 0:0.01:0.5 1e10],'linestyle','none');
% $$$     hold on;
% $$$ % $$$     [c,h] = contour(X,Y,Z,[-3:2:35],'-k');
% $$$ % $$$     clabel(c,h);
% $$$     quiver(lonu(xvec2,yvec2),latu(xvec2,yvec2),Z2(xvec2,yvec2),Z3(xvec2,yvec2),0.75,'-k');
% $$$     caxis(clim);
% $$$     xlim(xlims);
% $$$     ylim(ylims);
% $$$     hold on;
% $$$     plot([-150 -90 -90 -150 -150],[-5 -5 5 5 -5],'-','color',[0 0.5 ...
% $$$                         0],'linewidth',2);
% $$$     if (i>=5)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     else
% $$$         set(gca,'xticklabel',[]);
% $$$     end
% $$$     if (i==1 | i == 3 | i ==5 )
% $$$         ylabel('Latitude ($^\circ$N)');
% $$$     else
% $$$         set(gca,'yticklabel',[]);
% $$$     end
% $$$ % $$$     if (i>1)
% $$$         text(-257,-25,labels{i},'BackgroundColor','w');
% $$$ % $$$     else
% $$$ % $$$         text(-278,35,labels{i},'BackgroundColor','w');
% $$$ % $$$     end        
% $$$     if (i==2 | i == 4 | i ==6)
% $$$         cb = colorbar;
% $$$         ylabel(cb,'$^\circ$C');
% $$$     end
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'color','k');
% $$$     LabelAxes(gca,i,20,0.008,0.93);
% $$$ end 
% $$$ colormap(cmap);
% $$$ % $$$ colormap(parula);
% $$$ 
% $$$ 
% $$$ 
% $$$ %%% Heat convergence into layer plots:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [7];
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ wname = sprintf('/srv/ccrc/data03/z3500785/MOM_HeatDiag/output%03d/ocean_wmass.nc',outputs(1));
% $$$ 
% $$$ % $$$ Tw = 22.5; %Warm temperature
% $$$ % $$$ Tc = 15;   %Cool temperature
% $$$ Tw = 30; %Warm temperature
% $$$ Tc = 22.5;   %Cool temperature
% $$$ 
% $$$ [tmp iw] = min(abs(Te-Tw));
% $$$ [tmp ic] = min(abs(Te-Tc));
% $$$ 
% $$$ FlM = NaN*zeros(xL,yL,tL); % vdiffuse and nonlocal_KPP
% $$$ FlF = NaN*zeros(xL,yL,tL); % surface forcing
% $$$ FlP = NaN*zeros(xL,yL,tL); % P-E+R
% $$$ ty  = NaN*zeros(xL,yL,tL); % meridional mass flux
% $$$ tx  = NaN*zeros(xL,yL,tL); % zonal mass flux
% $$$ hy  = NaN*zeros(xL,yL,tL); % meridional heat flux
% $$$ hx  = NaN*zeros(xL,yL,tL); % zonal heat flux
% $$$ FlA = NaN*zeros(xL,yL,tL); % advection + submeso
% $$$ FlT = NaN*zeros(xL,yL,tL); % tendency
% $$$ 
% $$$ for ti=1:tL
% $$$     ii = iw;
% $$$     sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,iw-ii+1,iw-ic+1)
% $$$     FlT(:,:,ti) = ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlA(:,:,ti) = ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlP(:,:,ti) = ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlM(:,:,ti) = ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlF(:,:,ti) = ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     ty(:,:,ti)  = ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     tx(:,:,ti)  = ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     hy(:,:,ti)  = Te(ii)*Cp*ty(:,:,ti);
% $$$     hx(:,:,ti)  = Te(ii)*Cp*tx(:,:,ti);
% $$$ 
% $$$     for ii=iw-1:-1:ic
% $$$         sprintf('Calculating water-mass heat budget time %03d of %03d, temp %03d of %03d',ti,tL,iw-ii+1,iw-ic+1)
% $$$     FlT(:,:,ti) = FlT(:,:,ti)+ncread(wname,'temp_tendency_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlA(:,:,ti) = FlA(:,:,ti)+ncread(wname,'temp_advection_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_submeso_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlP(:,:,ti) = FlP(:,:,ti)+ncread(wname,'sfc_hflux_pme_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_rivermix_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlM(:,:,ti) = FlM(:,:,ti)+ncread(wname,'temp_vdiffuse_diff_cbt_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_nonlocal_KPP_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     FlF(:,:,ti) = FlF(:,:,ti)+ncread(wname,'temp_vdiffuse_sbc_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'sw_heat_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'frazil_on_nrho',[1 1 ii ti],[xL yL 1 1])+...
% $$$                   ncread(wname,'temp_eta_smooth_on_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     ty(:,:,ti)  = ty(:,:,ti) + ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     tx(:,:,ti)  = tx(:,:,ti) + ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     hy(:,:,ti)  = hy(:,:,ti) + Te(ii)*Cp*ncread(wname,'ty_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$     hx(:,:,ti)  = hx(:,:,ti) + Te(ii)*Cp*ncread(wname,'tx_trans_nrho',[1 1 ii ti],[xL yL 1 1]);
% $$$ 
% $$$     end
% $$$ end
% $$$ 
% $$$ save([base model sprintf('_output%03d',outputs(1)) '_VertInt_T' ...
% $$$       strrep(num2str(Tc),'.','p') 'C_T' strrep(num2str(Tw),'.','p') ...
% $$$       'C.mat'],'FlM','FlT','FlA','FlP','FlF','ty','tx','hy','hx','Tw','Tc');
% $$$ 
% $$$ load([base model sprintf('_output%03d',outputs(1)) '_VertInt_T' ...
% $$$       strrep(num2str(Tc),'.','p') 'C_T' strrep(num2str(Tw),'.','p') ...
% $$$       'C.mat']);
% $$$ FlA = nansum(FlA,3)/12;FlA(FlA==0) = NaN;
% $$$ FlT = nansum(FlT,3)/12;FlT(FlT==0) = NaN;
% $$$ FlM = nansum(FlM,3)/12;FlM(FlM==0) = NaN;
% $$$ FlF = nansum(FlF,3)/12;FlF(isnan(FlA)) = NaN;
% $$$ FlP = nansum(FlP,3)/12;FlP(isnan(FlA)) = NaN;
% $$$ ty = nansum(ty*1e9,3)/12;ty(ty==0) = NaN;
% $$$ tx = nansum(tx*1e9,3)/12;tx(tx==0) = NaN;
% $$$ hy = nansum(hy*1e9,3)/12;hy(hy==0) = NaN;
% $$$ hx = nansum(hx*1e9,3)/12;hx(hx==0) = NaN;
% $$$ 
% $$$ xvec = 1:1:xL;yvec = 1:1:yL;
% $$$ xvec2 = 1:5:xL-2;yvec2 = 1:5:yL-2;
% $$$ xvec3 = 1:8:xL;yvec3 = 1:4:yL;
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3    40   956   963]);
% $$$ subplot(3,1,1);
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlF(xvec,yvec)+FlP(xvec,yvec));
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title(['Surface Forcing heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
% $$$ caxis([-150 150]);
% $$$ set(gca,'color','k');
% $$$ ylim([-50 50]);
% $$$ subplot(3,1,2);
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlM(xvec,yvec));
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title(['Vertical Mixing heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
% $$$ caxis([-150 150]);
% $$$ set(gca,'color','k');
% $$$ ylim([-50 50]);
% $$$ subplot(3,1,3);
% $$$ pcolPlot(lon(xvec,yvec),lat(xvec,yvec),FlA(xvec,yvec)-FlT(xvec,yvec));
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title(['TEN-ADV heating between ' sprintf('%3.1f',Tc) '$^\circ$C and ' sprintf('%3.1f',Tw) '$^\circ$C']);
% $$$ caxis([-150 150]);
% $$$ set(gca,'color','k');
% $$$ ylim([-50 50]);
% $$$ colormap(redblue);
% $$$ subplot(4,1,4);
% $$$ X = avg(avg(lon,1),2);
% $$$ Y = avg(avg(lat,1),2);
% $$$ CF = pi*6371000/180; %conversion factor from lat/lon to distance.
% $$$ % $$$ hx(isnan(hx)) = 0;
% $$$ % $$$ hy(isnan(hy)) = 0;
% $$$ % $$$ DIV = avg(diff(hx,[],1)./diff(CF*lon.*cos(pi/180*lat),[],1),2) + ...
% $$$ % $$$       avg(diff(hy,[],2)./diff(CF*lat,[],2),1);
% $$$ DIV = avg(diff(hx,[],1)./avg(area,1),2) + ...
% $$$       avg(diff(hy,[],2)./avg(area,2),1);
% $$$ pcolPlot(X(xvec2,yvec2),Y(xvec2,yvec2),DIV(xvec2,yvec2));
% $$$ hold on;
% $$$ hand = quiver(lon(xvec3,yvec3),lat(xvec3,yvec3),hx(xvec3,yvec3),hy(xvec3,yvec3),20,'-k');
% $$$ % $$$ hand = quiver(lon(xvec3,yvec3),lat(xvec3,yvec3),hx(xvec3,yvec3)./sqrt(hx(xvec3,yvec3).^2+hy(xvec3,yvec3).^2),hy(xvec3,yvec3)./sqrt(hx(xvec3,yvec3).^2+hy(xvec3,yvec3).^2),'-k');
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title(['Heat Fluxes and divergence between %3.1fC and %3.1fC',Tc,Tw));
% $$$ caxis([-150 150]);
% $$$ % $$$ ylim([-50 50]);
% $$$ % $$$ ylim([-50 50]);
% $$$ xlim([-260 -70]);
% $$$ ylim([-40 40]);
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%% Components of mixing plot:
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(gcf,'defaulttextfontsize',17);
% $$$ set(gcf,'defaultaxesfontsize',17);
% $$$ 
% $$$ poss = [0.15  0.56 0.3 0.34; ...
% $$$         0.48  0.56 0.3 0.34; ...
% $$$         0.15  0.15 0.3 0.34; ...
% $$$         0.48  0.15 0.3 0.34];
% $$$ 
% $$$ obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ LAND = obj.SST(:,:,1);
% $$$ 
% $$$ clim = [-40 0];
% $$$ sp = 2;
% $$$ doWMT = 0; % plot WMT instead of flux
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ if (doWMT)
% $$$     cmap = flipud(lbmap(npts-3,'RedBlue'));
% $$$     cmap = redblue(npts-3);
% $$$     for i=1:(npts-3)
% $$$         if (cmap(i,:) == 1.0)
% $$$             cmap(i,:) = [0.94 0.94 0.94];
% $$$         end
% $$$     end
% $$$ else
% $$$     cmap = parula(npts-3);
% $$$     cmap = parula(npts-3);
% $$$     cmap(end,:) = [0.97 0.97 0.8];
% $$$     cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ end
% $$$ tmp = LAND;
% $$$ tmp(isnan(LAND)) = clim(1)-sp/2;
% $$$ tmp(~isnan(LAND)) = NaN;
% $$$ LAND = tmp;
% $$$ cmap(2:(end+1),:) = cmap;
% $$$ cmap(1,:) = [0 0 0];
% $$$ climn = [clim(1)-sp clim(2)];
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:1:xL;
% $$$ yvec = 1:1:yL;
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ months = [1:12];
% $$$ labels = {'(a) Background', ...
% $$$           '(b) Shear Instability', ...
% $$$           '(c) KPP Boundary Layer', ...
% $$$           '(d) Internal Wave'};
% $$$ VARS = {'FlMkppiw','FlMkppish','FlMkppbl','FlMwave'};
% $$$ 
% $$$ for ii=1:4
% $$$     subplot(2,2,ii);
% $$$     
% $$$ VAR = VARS{ii};
% $$$ TYPE = 'VertInt';
% $$$ Tl = 21.5;
% $$$ name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$ eval(['load(name,''' VAR ''');']);
% $$$ eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$ eval([VAR 'a = ' VAR ';']);
% $$$ for i=2:length(outputs)
% $$$     name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$     eval(['load(name,''' VAR ''');']);
% $$$     eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$     eval([VAR 'a = ' VAR 'a + ' VAR ';']);
% $$$ end
% $$$ eval([VAR ' = ' VAR 'a/length(outputs);']);
% $$$ eval([VAR '(' VAR '==0) = NaN;']);
% $$$ eval(['FlM = ' VAR ';']);
% $$$ 
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     if (length(months)>1)
% $$$         tmp = FlM;
% $$$         tmp(isnan(tmp)) = 0.0;
% $$$         Z = monmean(tmp(:,:,months),3,ndays(months));
% $$$         Z(Z == 0) = NaN;
% $$$     else
% $$$         Z = FlM(:,:,months);
% $$$     end
% $$$     Z = Z(xvec,yvec);
% $$$     if (doWMT)
% $$$         Z = Z*86400;
% $$$     end
% $$$     
% $$$     Z(Z<clim(1)) = clim(1);
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$     hold on;    
% $$$     contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
% $$$     caxis(climn);
% $$$     hold on;
% $$$ 
% $$$     if (ii>=3)
% $$$     xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (ii==1 | ii==3)
% $$$     ylabel('Latitude ($^\circ$N)');
% $$$     end
% $$$     if (ii==2 | ii==4)
% $$$         cb = colorbar;
% $$$         if (~doWMT)
% $$$             ylabel(cb,'Wm$^{-2}$');
% $$$         else
% $$$             ylabel(cb,'m/day');
% $$$         end            
% $$$         ylim(cb,clim);
% $$$     end
% $$$     text(-278,-40.5,labels{ii},'BackgroundColor','w','Margin',0.5);
% $$$ % $$$     title(labels{ii});
% $$$     ylim([-45 45]);
% $$$     set(gca,'ytick',[-45:15:45]);
% $$$     set(gca,'xtick',[-240:60:60]);
% $$$     set(gca,'Position',poss(ii,:));
% $$$ end 
% $$$ colormap(cmap);
% $$$ 
% $$$ %%% Plot shflux difference between runs with SST contours:
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = 94;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ % $$$ model = 'MOM025_kb3seg';
% $$$ % $$$ outputs = [75:79];
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux1 = monmean(shfluxa/length(outputs),3,ndays);
% $$$ SST1 = monmean(SSTa/length(outputs),3,ndays);
% $$$ clear shfluxa SSTa;
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ % $$$ outputs = [12]
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux2 = monmean(shfluxa/length(outputs),3,ndays);
% $$$ SST2 = monmean(SSTa/length(outputs),3,ndays);
% $$$ clear shfluxa SSTa;
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:1:xL;
% $$$ yvec = 291:1:706;
% $$$ 
% $$$ clim = [-25 25];
% $$$ sp = 2.5;
% $$$ 
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ poss = [0.1300    0.415  0.7693    0.56; ...
% $$$         0.1300    0.1    0.2343    0.2680; ...
% $$$         0.3951    0.1    0.2343    0.2680; ...
% $$$         0.6681    0.1    0.2343    0.2680];
% $$$ subplot(5,3,[1 9]);
% $$$ X = lon(xvec,yvec);
% $$$ Y = lat(xvec,yvec);
% $$$ Z = shflux1(xvec,yvec)-shflux2(xvec,yvec);
% $$$ contourf(X,Y,Z,cpts,'linestyle','none');
% $$$ hold on;
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec)-SST2(xvec,yvec),[-1.5:0.125:-0.125],'--k');
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec)-SST2(xvec,yvec),[0.125:0.125:1.5],'-k');
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec),[21.5 21.5],'-k','linewidth',3);
% $$$ % $$$ [c,h] = contour(X,Y,SST2(xvec,yvec),[21.5 21.5],'--','color',[0 0.5 ...
% $$$ % $$$                     0],'linewidth',2);
% $$$ caxis(clim);
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ ylim([-45 45]);
% $$$ set(gca,'ytick',[-45:15:45]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ set(gca,'xtick',[-240:60:60]);
% $$$ set(gca,'Position',[poss(1,:)]);
% $$$ set(gca,'color','k');
% $$$ colormap(cmap);
