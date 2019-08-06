
% This script takes a mean across multiple outputs of the ZA-processed
% .mats and saves back into a .mat file

close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
model = 'MOM025_kb3seg';
outputs = [101:110];
onum = 101110;

anavg = zeros(length(outputs)); % Take annual average

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

% Also do surface vars for min SST:
load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))],'SST');
SSTsave = SST;
if (anavg(1))
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
    if (~exist('ndays'))
        ndays = diff(time_snap);
    end
    SSTsave = monmean(SSTsave,3,ndays);
end

for i=2:length(outputs)
    i
    load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))],'SST');
    if (anavg(i))
        load([base model sprintf('_output%03d_BaseVars.mat',outputs(i))]);
        if (~exist('ndays'))
            ndays = diff(time_snap);
        end
        SST = monmean(SST,3,ndays);
    end

    SSTsave = SSTsave + SST;
end
SST = SSTsave/length(outputs);
save([base model sprintf('_output%03d_SurfaceVars.mat',onum)],'SST');

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


    
    
