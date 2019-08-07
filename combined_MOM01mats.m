
%This script combines the MOM01 3-month files into yearly files.
base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = '/g/data/e14/rmh561/MOM01_HeatDiag/mat_data/';
model = 'MOM01';

% $$$ for tttt=[1 2]
% $$$     if (tttt==1)
% $$$         outputs = [266 267 268 269];
% $$$         onum = 111;
% $$$     else
% $$$         outputs = [0 1 2 3];
% $$$         onum = 222;
% $$$     end
outputs = [4 5 6 7]
onum = 4567;

% BaseVars: ---------------------------
fname1 = [base model sprintf('_output%03d_BaseVars.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output%03d_BaseVars.mat',onum)]
copyfile(fname1,oname);

load(fname1);
timea = time;
ndaysa = ndays;
for i=2:length(outputs)
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
    timea = [timea; time];
    ndaysa = [ndaysa; ndays];
end
time = timea;
ndays = ndaysa;
tL = length(time);
save(oname,'time','ndays','tL','-append');

regions = {'Global','IndoPacific2BAS','Atlantic2BAS'};
for reg = 1:length(regions)
    region = regions{reg};

% GlobalHB: ---------------------------
fname1 = [base model sprintf(['_output%03d_' region '_HBud.mat'], ...
                             outputs(1))];
oname = [base model sprintf(['_output%03d_' region '_HBud.mat'],onum)]
copyfile(fname1,oname);

load(fname1);
GWBa = GWB;
names = fieldnames(GWBa);
for i=2:length(outputs)
    load([base model sprintf(['_output%03d_' region '_HBud.mat'],outputs(i))]);
    for ii=1:length(names)
        eval(['GWBa.' names{ii} ' = [GWBa.' names{ii} ' GWB.' names{ii} '];']);
    end
end
GWB = GWBa;
save(oname,'GWB');

% ZA: ---------------------------
fname1 = [base model sprintf('_output%03d_', ...
                             outputs(1)) region '_ZAHBud.mat'];
oname = [base model sprintf('_output%03d_',onum) region '_ZAHBud.mat']
copyfile(fname1,oname);

load(fname1);
ZAa = ZA;
names = fieldnames(ZAa);
for i=2:length(outputs)
    load([base model sprintf('_output%03d_', ...
                             outputs(i)) region '_ZAHBud.mat']);
    for ii=1:length(names)
        eval(['ZAa.' names{ii} ' = cat(3,ZAa.' names{ii} ',ZA.' names{ii} ');']);
    end
end
ZA = ZAa;
save(oname,'ZA','yu','yt');
end

% VertInt and WMT: ----------------------
type = {'VertInt'};
for typ = 1:length(type)
    for Ti=Te(1):(dT/2):Te(end)
        
        fname1 = [base model sprintf(['_output%03d_' type{typ} '_T' strrep(num2str(Ti),'.','p') 'C.mat'], ...
                             outputs(1))];
        if (exist(fname1))
            fname1
            oname = [base model sprintf(['_output%03d_' type{typ} '_T' strrep(num2str(Ti),'.','p') 'C.mat'],onum)]
            copyfile(fname1,oname);

            mat = load(fname1);
            mata = mat;
            names = fieldnames(mata);
            for i=2:length(outputs)
                mat = load([base model sprintf(['_output%03d_' type{typ} '_T' strrep(num2str(Ti),'.','p') 'C.mat'] ...
                             ,outputs(i))]);
                for ii=1:length(names)
                    eval(['sz = size(mat.' names{ii} ');']);
                    vec = find(sz==3);
                    if (length(vec)==1)
                        eval(['mata.' names{ii} ' = cat(vec,mata.' names{ii} ',mat.' names{ii} ');']);
                    end
                end
            end
            mat = mata;
            clear mata;
            for ii=1:length(names)
                eval([names{ii} ' = mat.' names{ii} ';']);
                save(oname,names{ii},'-append');
            end
        end
    end
end

% Any others that don't have temp -----------------------

type = {'SurfaceVars'}%,'varsat_110W','varsat_Eq'};
for typ = 1:length(type)

    fname1 = [base model sprintf(['_output%03d_' type{typ} '.mat'], ...
                                 outputs(1))];
    oname = [base model sprintf(['_output%03d_' type{typ} '.mat'],onum)]
    copyfile(fname1,oname);

    mat = load(fname1);
    mata = mat;
    names = fieldnames(mata);
    for i=2:length(outputs)
        mat = load([base model sprintf(['_output%03d_' type{typ} '.mat'] ...
                                       ,outputs(i))]);
        for ii=1:length(names)
            eval(['sz = size(mat.' names{ii} ');']);
            vec = find(sz==3);
            if (strcmp(names{ii},'u') | strcmp(names{ii},'v'))
                vec = 3;
            end
            if (length(vec)==1)
                eval(['mata.' names{ii} ' = cat(vec,mata.' names{ii} ',mat.' names{ii} ');']);
            end
        end
    end
    mat = mata;
    clear mata;
    for ii=1:length(names)
        eval([names{ii} ' = mat.' names{ii} ';']);
        if (ii==1)
            save(oname,names{ii},'-v7.3');
        else
            save(oname,names{ii},'-append');
        end
    end
end

end
