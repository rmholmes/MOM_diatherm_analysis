
base = '/short/e14/rmh561/access-om2/archive/';
pre = '1deg_jra55_ryf8485_';
post = '_may/';
runs = {'gfdl50','kds50','kds75','kds100','kds135'};
labs = {'GFDL50','KDS50','KDS75','KDS100','KDS135'};
color = {'k','b','r','m',[0 0.5 0]};

z = {};
dz = {};

for i=1:length(runs)
    name = [base pre runs{i} post 'output036/ocean/ocean.nc'];
    ztmp = ncread(name,'st_ocean');
    z{i} = (ztmp(1:end-1)+ztmp(2:end))/2;
    dz{i} = diff(ztmp);
end


figure;
set(gcf,'Color','w');
for i=1:length(runs)
    plot(dz{i},z{i},'o-','color',color{i},'linewidth',2, ...
         'MarkerSize',5);
    hold on;
end
set(gca,'ydir','reverse');
xlabel('\Delta z (m)');
ylabel('z (m)');
title('Vertical Resolution Comparison');
legend(labs);
ylim([0 300]);
xlim([0 50]);


    
