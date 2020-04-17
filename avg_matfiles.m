
% This script takes a mean across multiple outputs of .mat files
% and saves back into a single .mat file.

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
outputs = [101:120];
onum = 101120;

anavg = [ones(10,1); zeros(10,1)]; % Take annual average (choose these
                                   % manually for flexibility)

%%%% BaseVars (take from first output:
copyfile([base model sprintf('_output%03d_',outputs(1)) 'BaseVars.mat'], ...
         [base model sprintf('_output%03d_',onum) 'BaseVars.mat']);
load([base model sprintf('_output%03d_',onum) 'BaseVars.mat']);
if (find(anavg))         
    ndays = 1;
    tL = 1;
    time = 1;
    save([base model sprintf('_output%03d_',onum) 'BaseVars.mat'],'time','ndays','tL','-append');
end


%%%% Global budget

load([base model sprintf('_output%03d_',outputs(1)) 'Global_HBud.mat']);

GWBsave = GWB;
names = fieldnames(GWBsave);
if (anavg(1))
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    for ii=1:length(names)
        eval(['GWBsave.' names{ii} ' = monmean(GWBsave.' names{ii} ',2,ndays);']);
    end
end
    
for i=2:length(outputs)
    i
    load([base model sprintf('_output%03d_',outputs(i)) 'Global_HBud.mat']);
    
    if (anavg(i))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        for ii=1:length(names)
            eval(['GWB.' names{ii} ' = monmean(GWB.' names{ii} ',2,ndays);']);
        end
    end
        
    for fi=1:length(names)
        eval(['GWBsave.' names{fi} ' = GWBsave.' names{fi} ' + GWB.' names{fi} ';']);
    end
    
end

for fi=1:length(names)
    eval(['GWBsave.' names{fi} ' = GWBsave.' names{fi} '/length(outputs);']);
end
    
GWB = GWBsave;
save([base model sprintf('_output%03d_',onum) 'Global_HBud.mat'],'GWB');

%%%% Spatial fluxes

Tls = [0 5 10 15 20 22.5 25 27.5];
for ind=1:length(Tls)
    name = [base model sprintf('_output%03d',outputs(1)) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']
    str = load(name);
    str = rmfield(str,'Tl');
    strsave = str;
    names = fieldnames(strsave);
    if (anavg(1))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        for ii=1:length(names)
            eval(['strsave.' names{ii} ' = monmean(strsave.' names{ii} ',3,ndays);']);
        end
    end
    
    for i=2:length(outputs)
        name = [base model sprintf('_output%03d',outputs(i)) '_VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat']
        str = load(name);
        str = rmfield(str,'Tl');
    
        if (anavg(i))
            load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
            if (~exist('ndays'))
                ndays = diff(time_snap);
            end
            for ii=1:length(names)
                eval(['str.' names{ii} ' = monmean(str.' names{ii} ',3,ndays);']);
            end
        end
        
        for fi=1:length(names)
            eval(['strsave.' names{fi} ' = strsave.' names{fi} ' + str.' names{fi} ';']);
        end
    
    end

for fi=1:length(names)
    eval(['strsave.' names{fi} ' = strsave.' names{fi} '/length(outputs);']);
end
    
str = strsave;
str.Tl = Tls(ind);
save([base model sprintf('_output%03d_',onum) 'VertInt_T' strrep(num2str(Tls(ind)),'.','p') 'C.mat'],'-struct', 'str')

end

%%%% EqPM2
name = [base model sprintf('_output%03d',outputs(1)) '_varsat_EqPM2.mat']
str = load(name);
coords= {'Xt','Xu','Xw','Zt','Zu','Zw'};
for i=1:length(coords)
    eval([coords{i} ' = str.' coords{i} ';']);
    eval(['str=rmfield(str,''' coords{i} ''');']);
end
strsave = str;
names = fieldnames(strsave);
if (anavg(1))
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    for ii=1:length(names)
        eval(['lsz = length(size(strsave.' names{ii} '));']);
        if (lsz == 2)
            eval(['strsave.' names{ii} ' = monmean(strsave.' names{ii} ',2,ndays);']);
        else 
            eval(['strsave.' names{ii} ' = monmean(strsave.' names{ii} ',3,ndays);']);
        end            
    end
end
    
for i=2:length(outputs)
    name = [base model sprintf('_output%03d',outputs(i)) '_varsat_EqPM2.mat']
    str = load(name);
    coords= {'Xt','Xu','Xw','Zt','Zu','Zw'};
    for ii=1:length(coords)
        eval(['str=rmfield(str,''' coords{ii} ''');']);
    end
    names = fieldnames(str);
    if (anavg(i))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        for ii=1:length(names)
            eval(['lsz = length(size(str.' names{ii} '));']);
            if (lsz == 2)
                eval(['str.' names{ii} ' = monmean(str.' names{ii} ',2,ndays);']);
            else 
                eval(['str.' names{ii} ' = monmean(str.' names{ii} ',3,ndays);']);
            end            
        end
    end
    
    for fi=1:length(names)
        eval(['strsave.' names{fi} ' = strsave.' names{fi} ' + str.' names{fi} ';']);
    end
    
end

for fi=1:length(names)
    eval(['strsave.' names{fi} ' = strsave.' names{fi} '/length(outputs);']);
end
    
str = strsave;
for i=1:length(coords)
    eval(['str.' coords{i} ' = ' coords{i} ';']);
end
save([base model sprintf('_output%03d_',onum) 'varsat_EqPM2.mat'],'-struct', 'str')

%%%% ZA files
         
% $$$ regions = {'Atlantic2BAS','IndoPacific2BAS','Global'};
% $$$ regLets = {'A','P','G'};
regions = {'SO_Atlantic','SO_IndoPacific'};
regLets = {'SA','SP'};

for reg = 1:length(regions)
    region = regions{reg}
    regLet = regLets{reg};

    %% Make Vars
    type = 'ZA';
    load([base model sprintf('_output%03d_',outputs(1)) region '_' type 'HBud.mat']);

    ZAsave = ZA;
    names = fieldnames(ZAsave);
    if (anavg(1))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        for ii=1:length(names)
            eval(['ZAsave.' names{ii} ' = monmean(ZAsave.' names{ii} ',3,ndays);']);
        end
    end
    
    for i=2:length(outputs)
        i
        load([base model sprintf('_output%03d_',outputs(i)) region '_' type 'HBud.mat']);

        if (anavg(i))
            load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
            if (~exist('ndays'))
                ndays = diff(time_snap);
            end
            for ii=1:length(names)
                eval(['ZA.' names{ii} ' = monmean(ZA.' names{ii} ',3,ndays);']);
            end
        end
        
        for fi=1:length(names)
            eval(['ZAsave.' names{fi} ' = ZAsave.' names{fi} ' + ZA.' names{fi} ';']);
        end
    
    end

    for fi=1:length(names)
        eval(['ZAsave.' names{fi} ' = ZAsave.' names{fi} '/length(outputs);']);
    end
    
    ZA = ZAsave;
    save([base model sprintf('_output%03d_',onum) region '_' type 'HBud.mat'],'ZA','yt','yu');
end

%%% Surface vars:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
SSTsave = SST;
shfluxsave = shflux;
if (anavg(1))
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    SSTsave = monmean(SSTsave,3,ndays);
    shfluxsave = monmean(shfluxsave,3,ndays);
end

for i=2:length(outputs)
    i
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
    if (anavg(i))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        SST = monmean(SST,3,ndays);
        shflux = monmean(shflux,3,ndays);
    end

    SSTsave = SSTsave + SST;
    shfluxsave = shfluxsave + shflux;
end
SST = SSTsave/length(outputs);
shflux = shfluxsave/length(outputs);
save([base model sprintf('_output%03d_SurfaceVars.mat',onum)],'SST','shflux');


%%% XYtrans
% Do annual average XYtrans:
Tls = {'10','12p5','15','20','34'};
Tlsn = [10,12.5,15,20,34];
for ii=1:length(Tls)
    load([base model sprintf(['_output%03d_XYtrans_T' Tls{ii} 'C.mat'],outputs(1))]);
    xfluxT = xflux;
    yfluxT = yflux;
    for i=2:length(outputs)
        i
        load([base model sprintf(['_output%03d_XYtrans_T' Tls{ii} 'C.mat'],outputs(i))]);
        xfluxT = xfluxT+xflux;
        yfluxT = yfluxT+yflux;
    end
    xflux = xfluxT/length(outputs);
    yflux = yfluxT/length(outputs);
    Tl = Tlsn(ii)
    save([base model sprintf(['_output%03d_XYtrans_T' Tls{ii} 'C.mat'],onum)],'xflux','yflux','Tl');
end


    
    
