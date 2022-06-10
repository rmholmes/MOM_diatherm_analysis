% Calculate dV/dt from snapshots:
baseL = '/g/data/e14/rmh561/access-om2/archive/';
baseD = [baseL '1deg_jra55_ryf_fatITF/'];

region = 'Pacific';
post = 'ocean/';

output = 31;
restart = output-1;
base = [baseD sprintf('output%03d/',output) post];
basem1 = [baseD sprintf('output%03d/',output-1) post];

sname = [base 'ocean_snap.nc'];
wname = [base 'ocean_wmass.nc'];
gname = [base 'ocean_grid.nc'];

area = ncread(gname,'area_t');[xL,yL] = size(area);

zL = 50;
Cp = 3992.10322329649; % J kg-1 degC-1
rho0 = 1035; % kgm-3
T = ncread(wname,'neutral');
Te = ncread(wname,'neutralrho_edges');
TL = length(T);dT = T(2)-T(1);

mask_t = 'Your Pacific mask...';

Vsnap  = zeros(TL+1,1); % Volume from snapshots (m3)

for zi = 1:zL
    %Temperature snapshot:
    tempsnap1 = ncread(sname,'temp',[1 1 zi 1],[xL yL 1 1]);
    if (max(max(tempsnap))>120);tempsnap = tempsnap-273.15;end;

    % Volume snapshot:
    Volsnap = ncread(sname,'dzt',[1 1 zi 1],[xL yL 1 1]).*area;
    
    % Mask in Pacific-only:
    Volsnap(~mask_t) = 0;
    
    %Accumulate sums:
    for Ti=1:TL
        inds = find(tempsnap>=Te(Ti) & tempsnap<Te(Ti+1));
        Vsnap(Ti) = nansum(Volsnap(inds));
    end
    inds = find(tempsnap>=Te(TL+1));
    Vsnap(TL+1) = nansum(Volsnap(inds));
end
%Integrate to get to T'>T:
Vsnap(:,1) = flipud(cumsum(flipud(Vsnap)));

