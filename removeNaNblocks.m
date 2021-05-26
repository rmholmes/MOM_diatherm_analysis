%% Remove NaN's in the continents for display by copying in from
%% another BaseVar file:

% Original:
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
load([base model sprintf('_output%03d_BaseVars.mat',90)]);
region = 'Global';

% new run:
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'ACCESS-OM2_025deg_jra55_ryf_rediGM_kb1em5';
outputs = [7781];

for output = outputs
    save([base model sprintf('_output%03d_BaseVars.mat',output)], ...
         'lon','lat','lonu','latu','area','-append');
end

