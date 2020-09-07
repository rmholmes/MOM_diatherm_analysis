
% $$$ % Plot heat budget:
% $$$ load('ACCESS-CM2_PIcontrol_Closure.mat');
% $$$ 
% $$$ time = time/365.25;
% $$$ tL = length(time);
% $$$ 
% $$$ hfds_prior = cumsum(hfds_prior);
% $$$ hfds_new = cumsum(hfds_new);
% $$$ temp_tend = cumsum(temp_tend);
% $$$ OHC = OHC-OHC(1);
% $$$ hfds_prior = hfds_prior - hfds_prior(1);
% $$$ hfds_new = hfds_new - hfds_new(1);
% $$$ temp_tend = temp_tend - temp_tend(1);
% $$$ 
% $$$ % Annual average:
% $$$ tL_a = floor(tL/12);
% $$$ OHC_a = zeros(tL_a,1);
% $$$ hfds_prior_a = OHC_a;
% $$$ hfds_new_a = OHC_a;
% $$$ temp_tend_a = OHC_a;
% $$$ for ti=1:tL_a
% $$$     OHC_a(ti) = mean(OHC(((ti-1)*12+1):(ti*12)));
% $$$     hfds_prior_a(ti) = mean(hfds_prior(((ti-1)*12+1):(ti*12)));
% $$$     hfds_new_a(ti) = mean(hfds_new(((ti-1)*12+1):(ti*12)));
% $$$     temp_tend_a(ti) = mean(temp_tend(((ti-1)*12+1):(ti*12)));
% $$$ end
% $$$ time_a = 1:tL_a;
% $$$ 
% $$$ figure;
% $$$ plot(time_a,OHC_a,'-k','linewidth',2);
% $$$ hold on;
% $$$ plot(time_a,hfds_prior_a,'-r','linewidth',2);
% $$$ plot(time_a,hfds_new_a,'-b','linewidth',2);
% $$$ plot(time_a,temp_tend_a,'--m','linewidth',1);
% $$$ xlim([0 610]);
% $$$ ylim([-0.1 2]*1e24);
% $$$ xlabel('Year');
% $$$ ylabel('Heat Content [J]');
% $$$ legend('Total OHC (monthly averages)',...
% $$$        'Time-integrated hfds-prior',...
% $$$        'Time-integrated hfds-new',...
% $$$        'Time-integrated temp-tendency');
% $$$ 
% $$$        

% Plot mass budget:
load('ACCESS-CM2_PIcontrol_Mass_Closure.mat');

time = time/365.25;
tL = length(time);

pme_river = cumsum(pme_river);
pme_mass = cumsum(pme_mass);
pme_net = cumsum(pme_net);
wfimelt = cumsum(wfimelt);
wfiform = cumsum(wfiform);

OM = OM-OM(1);
OMeta = OMeta-OMeta(1);
pme_river = pme_river-pme_river(1);
pme_mass = pme_mass-pme_mass(1);
pme_net = pme_net-pme_net(1);
wfimelt = wfimelt-wfimelt(1);
wfiform = wfiform-wfiform(1);

% Annual average:
tL_a = floor(tL/12);
OM_a = zeros(tL_a,1);
OMeta_a = OM_a;
pme_river_a = OM_a;
pme_mass_a = OM_a;
pme_net_a = OM_a;
wfimelt_a = OM_a;
wfiform_a = OM_a;
for ti=1:tL_a
    OM_a(ti) = mean(OM(((ti-1)*12+1):(ti*12)));
    OMeta_a(ti) = mean(OMeta(((ti-1)*12+1):(ti*12)));
    pme_river_a(ti) = mean(pme_river(((ti-1)*12+1):(ti*12)));
    pme_mass_a(ti) = mean(pme_mass(((ti-1)*12+1):(ti*12)));
    pme_net_a(ti) = mean(pme_net(((ti-1)*12+1):(ti*12)));
    wfimelt_a(ti) = mean(wfimelt(((ti-1)*12+1):(ti*12)));
    wfiform_a(ti) = mean(wfiform(((ti-1)*12+1):(ti*12)));
end
time_a = 1:tL_a;

figure;
plot(time_a,OM_a,'-k','linewidth',2);
hold on;
plot(time_a,OMeta_a,'--r','linewidth',2);
plot(time_a,pme_river_a,'-b','linewidth',2);
% $$$ plot(time_a,pme_mass_a,'-g','linewidth',2);
% $$$ plot(time_a,pme_net_a,'-y','linewidth',2);
% $$$ plot(time_a,wfimelt_a,'-c','linewidth',2);
% $$$ plot(time_a,wfiform_a,'-m','linewidth',2);
xlim([0 610]);
ylim([-50 10]*1e15);
xlabel('Year');
ylabel('Mass [kg]');
legend('Total OM (monthly averages)',...
       'Total OM from sea-level (monthly averages)',...
       'Time-integrated pme,river,ice etc.');
% $$$        'Time-integrated pme-mass',...
% $$$        'Time-integrated pme-net',...
% $$$        'Time-integrated wfimelt',...
% $$$        'Time-integrated wfiform');

       

