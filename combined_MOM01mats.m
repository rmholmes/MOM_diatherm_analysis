
%This script combines the MOM01 3-month files into yearly files.
% $$$ base = '/srv/ccrc/data03/z3500785/MOM01_HeatDiag/mat_data/';
base = '/g/data/e14/rmh561/MOM01_HeatDiag/mat_data/';
model = 'MOM01';

for tttt=[1 2]
    if (tttt==1)
        outputs = [266 267 268 269];
        onum = 111;
    else
        outputs = [0 1 2 3];
        onum = 222;
    end

% BaseVars: ---------------------------
fname1 = [base model sprintf('_output%03d_BaseVars.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output%03d_BaseVars.mat',onum)]
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

% GlobalHB: ---------------------------
fname1 = [base model sprintf('_output%03d_GlobalHBud.mat', ...
                             outputs(1))];
oname = [base model sprintf('_output%03d_GlobalHBud.mat',onum)]
copyfile(fname1,oname);

load(fname1);
GWBa = GWB;
names = fieldnames(GWBa);
for i=2:length(outputs)
    load([base model sprintf('_output%03d_GlobalHBud.mat',outputs(i))]);
    for ii=1:length(names)
        eval(['GWBa.' names{ii} ' = [GWBa.' names{ii} ' GWB.' names{ii} '];']);
    end
end
GWB = GWBa;
save(oname,'GWB');

% VertInt and WMT: ----------------------
type = {'VertInt','WMT'};
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

type = {'SurfaceVars','varsat_110W','varsat_Eq'};
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
