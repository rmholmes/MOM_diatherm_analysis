
%This script combines the MOM01 3-month files into yearly files.
base = '/srv/ccrc/data03/z3500785/MOM01_HeatDiag/mat_data/';
model = 'MOM01';
outputs = [266 267 268 269];

fname1 = [base model sprintf('_output%03d_BaseVars.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output333_BaseVars.mat')];
copyfile(fname1,oname);

load(fname1);
timea = time;
time_snapa = time_snap;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
    timea = [timea; time];
    time_snapa = [time_snapa; time_snap];
end
time = timea;
time_snap = time_snapa([1:4 6:8 10:12 14:16]);
tL = length(time);
save(oname,'time','time_snap','tL','-append');


fname1 = [base model sprintf('_output%03d_GlobalHBud.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output333_GlobalHBud.mat')];
copyfile(fname1,oname);

vars = {'ADV','ETS','FRZ','H','Hsnap','KNL','PME','RMX','SFW', ...
        'SUB','SWH','TEN','TENMON','Temp','V','VDF','VDS','Vsnap'};
load(fname1);
for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
end
for i=2:length(outputs)
    load([base model sprintf('_output%03d_GlobalHBud.mat',outputs(i))]);
    for ii=1:length(vars)
        eval([vars{ii} 'a = [' vars{ii} 'a ' vars{ii} '];']);
    end
end
%Fix Vsnap and Hsnap:
Vsnapa = Vsnapa(:,[1:4 6:8 10:12 14:16]);
Hsnapa = Hsnapa(:,[1:4 6:8 10:12 14:16]);

for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a;']);
    eval(['save(oname,''' vars{i} ''',''-append'');']);
end

fname1 = [base model sprintf('_output%03d_VertInt_T22p5C.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output333_VertInt_T22p5C.mat')];
copyfile(fname1,oname);

load(fname1);
FlMa = FlM;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_VertInt_T22p5C.mat',outputs(i))]);
    FlMa = cat(3,FlMa,FlM);
end
FlM = FlMa;
save(oname,'FlM','-append');

fname1 = [base model sprintf('_output%03d_SurfaceVars.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output333_SurfaceVars.mat')];
copyfile(fname1,oname);

load(fname1);
SSTa = SST;
shfluxa = shflux;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    SSTa = cat(3,SSTa,SST);
    shfluxa = cat(3,shfluxa,shflux);
end
SST = SSTa;
save(oname,'SST','-append');
shflux = shfluxa;
save(oname,'shflux','-append');



fname1 = [base model sprintf('_output%03d_varsat_140W.mat', ...
 outputs(1))];
oname = [base model sprintf('_output333_varsat_140W.mat')];
copyfile(fname1,oname);

vars = {'temp','kappa','taux','tauy','mld','vdif','vnlc'};
load(fname1);
for i=1:length(vars)
    eval([vars{i} 'a = ' vars{i} ';']);
    ua = u;
    va = v;
end
for i=2:length(outputs)
    load([base model sprintf('_output%03d_varsat_140W.mat',outputs(i))]);
    for ii=1:length(vars)
        eval([vars{ii} 'a = cat(length(size(' vars{ii} ')),' vars{ii} 'a,' vars{ii} ');']);
    end
    ua = cat(3,ua,u);
    va = cat(3,va,v);
end

for i=1:length(vars)
    eval([vars{i} ' = ' vars{i} 'a;']);
    eval(['save(oname,''' vars{i} ''',''-append'');']);
end
u = ua;
v = va;
save(oname,'u','-append');
save(oname,'v','-append');

