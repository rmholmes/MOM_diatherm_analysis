% This script makes line plots of specific regions

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';

RUNS = {'ACCESS-OM2_025deg_jra55_ryf_norediGM', ...
        'ACCESS-OM2_025deg_jra55_ryf_noGM', ...
        'ACCESS-OM2_025deg_jra55_ryf', ...
       };
colors = {'r','r','r'};
styles = {'-',':','--'};
names = {'ACCESS-OM2-025','ACCESS-OM2-025-N','ACCESS-OM2-025-NG'};
region = 'Kuroshio';
region = 'GulfStream';
Tcaxs = [-2 30];

% $$$ RUNS = {'ACCESS-OM2_025deg_jra55_ryf', ...
% $$$         'ACCESS-OM2_025deg_jra55_ryf_rediGM_kbvar', ...
% $$$         'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5', ...
% $$$        };
% $$$ colors = {'r',[0 0.5 0],'k'};
% $$$ styles = {'-',':','--'};
% $$$ names = {'ACCESS-OM2-025-NG','ACCESS-OM2-025-NG-kbv','ACCESS-OM2-025-NG-kb5'};
% $$$ region = 'Kuroshio';
% $$$ region = 'EastTropPac';
% $$$ Tcaxs = [10 30];

RUNS = {'ACCESS-OM2_025deg_jra55_ryf_norediGM', ...
        'ACCESS-OM2_025deg_jra55_ryf_norediGM_smoothkppbl'
       };
colors = {'r',[0 0.5 0],'k'};
styles = {'-',':','--'};
names = {'ACCESS-OM2-025','ACCESS-OM2-025-KPPBL-smooth'};
region = 'Kuroshio';
region = 'EastTropPac';
Tcaxs = [10 30];

rr = 1;
NUMs = {};
udhs = {};
vdhs = {};
Tdzs = {};
for rr = 1:length(RUNS)
    model = RUNS{rr}
    load([base model '_' region '_nummixprofiles.mat']);

    % time-averaging:
    NUMs{rr} = mean(mean(NUM,2),3);
    vars = {'EKE','Tdxsq','Tdysq','Tdzsq','aiso','rey_bih', ...
            'udxsq','vdysq'};
    vars = {'Tdxsq','Tdysq','Tdzsq','udxsq','vdysq'};
    for vi = 1:length(vars)
        eval([vars{vi} ' = monmean(' vars{vi} ',1,ndays);']);
    end
    
    udhs{rr} = sqrt(0.5*(udxsq+vdysq));
    Tdhs{rr} = sqrt(0.5*(Tdxsq+Tdysq));
    Tdzs{rr} = sqrt(Tdzsq);
    EKEs{rr} = EKE;
end

poss = [0.05    0.0900    0.16    0.8150; ...
        0.23    0.0900    0.16    0.8150; ...
        0.41    0.0900    0.16    0.8150; ...
        0.59    0.0900    0.16    0.8150; ...
        0.77    0.0900    0.16    0.8150];
    
    figure;
    set(gcf,'Position',[1          36        1920         970]);
    subplot(1,5,1);
    for rr = 1:length(RUNS)
        plot(NUMs{rr}/1e12,Te,styles{rr},'color',colors{rr},'linewidth',2);
        hold on;
    end
    ylabel('Temperature $\Theta^*$ ($^\circ$C)');
    title('Numerical Mixing $\mathcal{I}(\Theta^*)$');
    xlabel('TW');
    xlim([-80 0]);
% $$$     xlim([-60 0]);
    ylim(Tcaxs);
    legend(names);
    set(gca,'Position',poss(1,:));
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    text(xlims(2)-0.12*diff(xlims),ylims(2)-0.03*diff(ylims),'(a)');

    subplot(1,5,2);
    for rr = 1:length(RUNS)
        plot(EKEs{rr},T,styles{rr},'color',colors{rr},'linewidth',2);
        hold on;
    end
    title('EKE');
    xlabel('m$^2$s$^{-2}$');
    set(gca,'yticklabel',[]);
    xlim([0 4]*1e-3);
% $$$     xlim([0 10]*1e-3);
    ylim(Tcaxs);
    set(gca,'Position',poss(2,:));
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    text(xlims(2)-0.12*diff(xlims),ylims(2)-0.03*diff(ylims),'(b)');
    
    subplot(1,5,3);
    for rr = 1:length(RUNS)
        plot(udhs{rr},T,styles{rr},'color',colors{rr},'linewidth',2);;
        hold on;
    end
% $$$     ylabel('Temperature ($^\circ$C)');
    set(gca,'yticklabel',[]);
    title('$\sqrt{\frac{1}{2}\left(\overline{|\Delta_x u|^2}+\overline{|\Delta_y v|^2}\right)}$');
    xlabel('ms$^{-1}$');
    ylim(Tcaxs);
    set(gca,'Position',poss(3,:));
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    text(xlims(2)-0.12*diff(xlims),ylims(2)-0.03*diff(ylims),'(c)');

    subplot(1,5,4);
    for rr = 1:length(RUNS)
        plot(Tdhs{rr},T,styles{rr},'color',colors{rr},'linewidth',2);;
        hold on;
    end
% $$$     ylabel('Temperature ($^\circ$C)');
    set(gca,'yticklabel',[]);
    title('$\sqrt{\frac{1}{2}\left(\overline{|\Delta_x \Theta|^2}+\overline{|\Delta_y \Theta|^2}\right)}$');
    xlabel('$^\circ$C');
    ylim(Tcaxs);
    set(gca,'Position',poss(4,:));
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    text(xlims(2)-0.12*diff(xlims),ylims(2)-0.03*diff(ylims),'(d)');
    
    subplot(1,5,5);
    for rr = 1:length(RUNS)
        plot(Tdzs{rr},T,styles{rr},'color',colors{rr},'linewidth',2);;
        hold on;
    end
% $$$     ylabel('Temperature ($^\circ$C)');
    set(gca,'yticklabel',[]);
    title('$\sqrt{\overline{|\Delta_z \Theta|^2}}$');
    xlabel('$^\circ$C');
    xlim([0 15]);
% $$$     xlim([0 8]);
    ylim(Tcaxs);
    set(gca,'Position',poss(5,:));
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    text(xlims(2)-0.12*diff(xlims),ylims(2)-0.03*diff(ylims),'(e)');

