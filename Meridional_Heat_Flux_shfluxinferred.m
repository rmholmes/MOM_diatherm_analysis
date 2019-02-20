% This script makes a plot of the meridional heat flux from MOM-SIS
% simulations inferred from the surface heat flux

close all;
clear all;

baseL = '/short/e14/rmh561/mom/archive/';
model = 'MOM025_kb3seg';
baseD = [baseL 'MOM_HeatDiag_kb3seg/']; %Data Directory.
output = 86;
base = [baseD sprintf('output%03d/',output)];
fname = [base 'ocean.nc'];
gname = [base 'ocean_grid.nc'];

lat = ncread(gname,'geolat_t');
area = ncread(gname,'area_t');

ndays = ncread(fname,'average_DT');
shflux = ncread(fname,'net_sfc_heating');
[xL,yL,tL] = size(shflux);

shflux = sum(shflux.*repmat(permute(ndays(:),[3 2 1]),[xL yL 1]),3)/sum(ndays);

% Calculate meridional heat flux inferred:
latV = linspace(-90,90,181);
V = zeros(size(latV));
for i=1:length(latV)
    inds = lat < latV(i);
    V(i) = nansum(area(inds).*shflux(inds));
end

%Center the flux:
V = V + (V(1)-V(end))/2;

figure;

plot(latV,V/1e15,'-r','linewidth',2);
xlabel('Latitude ($^\circ$N)');
ylabel('Meridional Heat Flux (PW)');
grid on;
box on;
xlim([-90 90]);
ylim([-1 2]);
set(gca,'xtick',[-90:30:90]);
