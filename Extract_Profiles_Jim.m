
% Extract profiles of the vertical heat flux and the temperature
% from the MOM025 Control simulations of Holmes et al. 2019 JPO at
% specific locations.

outputs = [101:105]; % outputs to use

locations = [-10 0; % 10W, 0N
             -23 0; % 23W, 0N
             -140 0; % 140W, 0N
             -110 0; % 110W, 0N
             ]

out = outputs(1);
outb = sprintf('output%03d/',out);

% Load some grid-variables:
lon = ncread([outb 'ocean.nc'],'xt_ocean');
lat = ncread([outb 'ocean.nc'],'yt_ocean');
T = ncread([outb 'ocean_wmass.nc'],'neutralrho_edges');
Z = ncread([outb 'ocean.nc'],'st_ocean');

zL = 50;
tL = 12;
TL = length(T);

for li = 1:length(locations)
    [tmp lnind] = min(abs(lon-locations(li,1)));
    [tmp ltind] = min(abs(lat-locations(li,2)));

    temp = zeros(zL,tL);
    t_vdiff = zeros(TL-1,tL);
    
    for oi = 1:length(outputs)
        out = outputs(oi)
        [num2str(li) ' - ' num2str(out)]
        
        outb = sprintf('output%03d/',out);

        temp = temp + squeeze(ncread([outb 'ocean.nc'],'temp',[lnind ...
                            ltind 1 1],[1 1 zL tL]));

        t_vdiff = t_vdiff + squeeze(ncread([outb 'ocean_wmass.nc'],'temp_vdiffuse_diff_cbt_on_nrho', ...
                                           [lnind ltind 1 1],[1 1 TL-1 tL]));
    end
    t_vdiff = t_vdiff/length(outputs);
    temp = temp/length(outputs);

    JqT = cat(1,zeros(1,tL),cumsum(t_vdiff,1));

    name = ['Holmesetal2019_' num2str(-locations(li,1)) 'W_JqT.mat'];
    save(name,'Z','T','temp','JqT');
end

