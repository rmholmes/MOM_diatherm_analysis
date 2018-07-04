% This script makes plots of the spatial structure of the
% diathermal fluxes in the MOM simulations.

% $$$ close all;
clear all;

base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ base = 'archive/mat_data/';

RUNS = { ...
% MOM01-SIS:
% $$$     {'MOM01',[222]}, ...
% $$$ % MOM025-SIS:
% $$$     {'MOM025',[8:12]}, ...
% $$$     {'MOM025',[15:19]}, ...
% $$$ % $$$     {'MOM025_kb1em6',[30]}, ...
% $$$     {'MOM025_kb3seg',[80:84]}, ...
% $$$     {'MOM025_kb1em5',[94]}, ...
% $$$     {'MOM025_wombat',[1978]}, ...
% ACCESS-OM2 025-degree:
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485',[78]}, ...
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485_redi',[59]}, ...
% $$$     {'ACCESS-OM2_025deg_jra55_ryf8485_gmredi',[73]}, ...
% $$$ %     {'ACCESS-OM2_025deg_jra55_ryf8485_KDS75',[??]}, ...
% ACCESS-OM2 1-degree:
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_Tcen',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_TcenGMS',[36]}, ...
         {'ACCESS-OM2_1deg_jra55_ryf8485_gfdl50_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds75_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds100_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds135_may',[36]}, ...
% $$$          {'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may_kb1em5',[0]}, ...
       };

rr = 1;
for rr = 1:length(RUNS);
    outputs = RUNS{rr}{2};
    model = RUNS{rr}{1};

    clearvars -except base RUNS rr outputs model;
    
    load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
    ndays = diff(time_snap);
    ndays = ndays(1:12);
    if (ndays(end) <= 0); ndays(end) = 365-ndays(end);end;
    region = 'Global';
% $$$ region = 'Pacific';
    nyrs = tL/12;szTe = [TL+1 12 nyrs];szT  = [TL 12 nyrs];
    yrs = 1:nyrs;
    months = 1:12;
    
    ycur = 1;


    %%% Spatial Structure:
    VAR = 'FlI';
    TYPE = 'VertInt';
    Tl = 27.5;
    name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
    eval(['load(name,''' VAR ''');']);
    eval([VAR '(isnan(' VAR ')) = 0.0;']);
    if (length(outputs)==1)
        eval([VAR ' = reshape(' VAR ',[length(' VAR '(:,1,1)) length(' VAR '(1,:,1)) 12 nyrs]);']);
        eval([VAR ' = mean(' VAR '(:,:,:,yrs),4);']);
    else
        eval([VAR 'a = ' VAR ';']);
        for i=2:length(outputs)
            name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
            eval(['load(name,''' VAR ''');']);
            eval([VAR '(isnan(' VAR ')) = 0.0;']);
            eval([VAR 'a = ' VAR 'a + ' VAR ';']);
        end
        eval([VAR ' = ' VAR 'a/length(outputs);']);
        eval([VAR '(' VAR '==0) = NaN;']);
    end
    eval(['FlM = ' VAR ';']);

% $$$ % CHECK spatial structure sums to total:
% $$$ Tls = [0:2.5:27.5];
% $$$ SUM = zeros(size(Tls));
% $$$ for ii = 1:length(Tls)
% $$$ 
% $$$     Tl = Tls(ii)
% $$$     name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$     eval(['load(name,''' VAR ''');']);
% $$$     eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$     if (length(outputs)==1)
% $$$         eval([VAR ' = reshape(' VAR ',[length(' VAR '(:,1,1)) length(' VAR '(1,:,1)) 12 nyrs]);']);
% $$$         eval([VAR ' = mean(' VAR '(:,:,:,yrs),4);']);
% $$$     else
% $$$         eval([VAR 'a = ' VAR ';']);
% $$$         for i=2:length(outputs)
% $$$             name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$             eval(['load(name,''' VAR ''');']);
% $$$             eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$             eval([VAR 'a = ' VAR 'a + ' VAR ';']);
% $$$         end
% $$$         eval([VAR ' = ' VAR 'a/length(outputs);']);
% $$$         eval([VAR '(' VAR '==0) = NaN;']);
% $$$     end
% $$$     eval(['FlM = ' VAR ';']);
% $$$     tmp = FlM;
% $$$     tmp(isnan(tmp)) = 0.0;
% $$$     Z = monmean(tmp(:,:,months),3,ndays(months));
% $$$     Z(Z == 0) = NaN;
% $$$     SUM(ii) = nansum(nansum(area.*Z));
% $$$ end
% $$$ plot(Tls,SUM/1e15,'Xb','MarkerSize',12,'LineWidth',2); 
% $$$ %plot(Tls,SUM/1e15,'Xb','MarkerSize',12,'LineWidth',2); 
% $$$
% $$$ %%% Regional time series 
% $$$ % $$$ 
% $$$ %[mask_t,~] = Heat_Budget_Mask('Pacific','','','',base,'MOM025');
% $$$ mask_t = zeros(size(lon));
% $$$ months = [1:12];
% $$$ %Region choice:
% $$$ regions = { ...
% $$$     {'Global',lat==lat,'k'}, ...
% $$$     {'Equatorial',abs(lat)<=10,'k'}, ...
% $$$     {'Eastern Pacific',lat > -10 & lat < 10 & lon > -160 & lon<-70,[0.4667    0.6745    0.1882]}, ...
% $$$     {'Kuroshio', lat < 50 & lat > 10 & lon > -260 & lon < -140,'r'}, ...
% $$$     {'Gulf Stream', lat < 50 & lat > 10 & lon > -100 & lon < 0 & ~mask_t,'b'}, ...
% $$$     {'Indian', lat < 25 & lat > -55 & (lon > 20 | lon < -260),[0.4941    0.1843    0.5569]}, ...
% $$$     {'South Pacific/Atlantic', lat < -10 & lat > -55 & lon > -260 & lon < 20,'g'}, ...
% $$$           };
% $$$ % $$$     {'Eastern Pacific',lat > -20 & lat < 20 & lon > -150 & lon<-60 & mask_t,[0.4667    0.6745    0.1882]}, ...
% $$$ % $$$     {'WWV',abs(lat)<=5 & lon>
% $$$ % $$$     {'NH',lat>10,'k'}, ...
% $$$ % $$$     {'SH',lat<-10,'k'}, ...
% $$$ % $$$     {'Outside Equatorial',lat<-10 | lat > 10,'k'}, ...
% $$$ % $$$     {'Nino 3',lat > -5 & lat < 5 & lon > -150 & lon < -90,'m'}, ...
% $$$ % $$$     {'Eastern Pacific',lat > -10 & lat < 10 & lon > -160 & mask_t,'c'}, ...
% $$$ % $$$     {'Equatorial Atlantic',lat > -10 & lat < 10 & lon > -65 & lon < 20,[0 0.5 0]}, ...
% $$$ % $$$     {'Western Pacific',lat > -10 & lat < 10 & lon > -260 & lon < -160,[0.5 0.5 0.5]}, ...
% $$$ % $$$     {'Gulf Stream', lat < 45 & lat > 10 & lon > -100 & lon < -40 & ~mask_t,'b'}, ...
% $$$ 
% $$$ Field = zeros(12,length(regions));
% $$$ AREA = zeros(12,length(regions));
% $$$ for i=1:12
% $$$     Mtmp = FlM(:,:,i);
% $$$     Atmp = area;
% $$$     Atmp(isnan(Mtmp)) = NaN;
% $$$     
% $$$     for ii=1:length(regions)
% $$$         Field(i,ii) = nansum(nansum(area(regions{ii}{2}).*Mtmp(regions{ii}{2}),1),2);
% $$$         AREA(i,ii) = nansum(nansum(Atmp(regions{ii}{2}),1),2);
% $$$     end
% $$$ end
% $$$ 
% $$$ %Display output in terminal for table:
% $$$ str = {['Area fractions model ' model ' Temp ' num2str(Tl)] ;
% $$$        ' '};
% $$$ for ii = 1:length(regions)
% $$$     str{ii+2} = sprintf([regions{ii}{1} ' area = %3.2f'], ...
% $$$                         monmean(AREA(:,ii),1,ndays(months))/monmean(AREA(:,1),1,ndays(months)));
% $$$ end
% $$$ str
% $$$ 
% $$$ str = {['Annual totals model ' model ' Temp ' num2str(Tl)] ;
% $$$        ' '};
% $$$ for ii = 1:length(regions)
% $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.2fPW (%3.0f)', ...
% $$$                         monmean(Field(:,ii),1,ndays(months))/1e15,monmean(Field(:,ii),1,ndays(months))/monmean(Field(:,1),1,ndays(months))*100)];
% $$$ end
% $$$ str
% $$$ 
% $$$ str = {['SC range model ' model ' Temp ' num2str(Tl)] ;
% $$$        ' '};
% $$$ for ii = 1:length(regions)
% $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.2fPW (%3.0f)', ...
% $$$                         (max(Field(:,ii))-min(Field(:,ii)))/1e15,(max(Field(:,ii))-min(Field(:,ii)))/(max(Field(:,1))-min(Field(:,1)))*100)];
% $$$ end
% $$$ str
% $$$ 
% $$$ % $$$ mn1 = 4;mn2 = 7;
% $$$ % $$$ str = {['SC range Apr-Jul model ' model ' Temp ' num2str(Tl)] ;
% $$$ % $$$        ' '};
% $$$ mn1 = 5;mn2 = 8;
% $$$ str = {['SC range May-Aug model ' model ' Temp ' num2str(Tl)] ;
% $$$        ' '};
% $$$ for ii = 1:length(regions)
% $$$     str{ii+2} = [regions{ii}{1} sprintf(' = %3.5fPW (%3.0f)', ...
% $$$                         (Field(mn1,ii)-Field(mn2,ii))/1e15,(Field(mn1,ii)-Field(mn2,ii))/(Field(mn1,1)-Field(mn2,1))*100)];
% $$$ end
% $$$ str
% $$$ 
% $$$ 
% $$$ % The 1% (within Nino 3) area count of the annual mean:
% $$$ inds = lat > -5 & lat < 5 & lon > -150 & lon < -90;
% $$$ M1p = zeros(length(find(inds)),12);
% $$$ A1p = area(inds);
% $$$ AREAtotal = zeros(12,1);
% $$$ for i=1:12
% $$$     Mtmp = FlM(:,:,i);
% $$$     M1p(:,i) = A1p.*Mtmp(inds);
% $$$     Atmp = area;
% $$$     Atmp(isnan(Mtmp)) = NaN;
% $$$     AREAtotal(i) = nansum(nansum(Atmp,1),2);
% $$$ end
% $$$ M1pA = monmean(M1p,2,ndays(months));
% $$$ [M1pA,I] = sort(M1pA);
% $$$ M1p = M1p(I,:);
% $$$ Atotal = monmean(AREAtotal,1,ndays(months));
% $$$ [tmp ind] = min(abs(cumsum(A1p(I))/Atotal - 0.005));
% $$$ %plot(cumsum(A1p(I))/Atotal*100,cumsum(M1pA)/monmean(MAll,2,ndays(months))*100)
% $$$ M1pA = sum(M1pA(1:ind));
% $$$ M1pSCR = (sum(M1p(1:ind,mn1))-sum(M1p(1:ind,mn2)));
% $$$ MAll = Field(:,1)';
% $$$ 
% $$$ str = {['1% area within Nino 3 model ' model ' Temp ' num2str(Tl)]  ; ...
% $$$        sprintf(' 1p Annual-Mean = %3.2fPW (%3.0f)',M1pA/1e15,M1pA/monmean(MAll,2,ndays(months))*100); ...
% $$$        sprintf(' 1p SC Apr-Jul = %3.2fPW (%3.0f)',M1pSCR/1e15,M1pSCR/((MAll(mn1))-(MAll(mn2)))*100)}
% $$$ 
% $$$ % $$$ % Average fluxes, isotherm separationl, diffusivity:
% $$$ % $$$ load([base 'MOM025_output002_Ziso_T23C.mat']);
% $$$ % $$$ topiso = ziso;
% $$$ % $$$ load([base 'MOM025_output002_Ziso_T22C.mat']);
% $$$ % $$$ botiso = ziso;
% $$$ % $$$ AvgFlux = MAll./AREAtotal;
% $$$ % $$$ NaNs = isnan(topiso-botiso);
% $$$ % $$$ AvgDZ = zeros(12,1);
% $$$ % $$$ for ii=1:12
% $$$ % $$$     AvgDZ(ii) = nansum(nansum(((topiso(:,:,ii)-botiso(:,:,ii)).*area),1),2)./ ...
% $$$ % $$$                 nansum(area(~isnan((topiso(:,:,ii)-botiso(:,:,ii)))));
% $$$ % $$$ end
% $$$ % $$$ AvgDiff = monmean(AvgFlux/rho0/Cp.*AvgDZ',2,ndays(months))
% $$$ % $$$ AvgFlux = monmean(AvgFlux,2,ndays(months))
% $$$ % $$$ AvgDZ = monmean(AvgDZ',2,ndays(months))
% $$$ % $$$ sprintf('Average flux across the isotherm = %3.1f Wm-2',AvgFlux)
% $$$ % $$$ sprintf('Average diffusivity = %3.5f m2s-1',AvgDiff)
% $$$ % $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ 
% $$$ legstr = {};
% $$$ for ii=1:length(regions)
% $$$     plot(1:12,Field(:,ii)/1e15,'-','color',regions{ii}{3},'LineWidth',2);
% $$$     hold on;
% $$$     legstr{ii} = regions{ii}{1};
% $$$ end
% $$$ % $$$ plot(1:12,MEEP/1e15,'--r','LineWidth',4);
% $$$ % $$$ hold on;
% $$$ % $$$ plot(1:12,MEq/1e15,'--k','LineWidth',4);
% $$$ % $$$ plot(1:12,MN/1e15,':b','LineWidth',4);
% $$$ % $$$ plot(1:12,MS/1e15,':','color',[0 0.5 0],'LineWidth',4);
% $$$ % $$$ plot(1:12,(MN+MS)/1e15,':k','LineWidth',4);
% $$$ % $$$ plot(1:12,MAll/1e15,'-k','LineWidth',4);
% $$$ xlabel('Month');
% $$$ ylabel(['PW']);
% $$$ % $$$ title(['Vertical mixing heat flux through $' num2str(Tl) '^\circ$C isotherm']);
% $$$ title(['Implicit mixing heat flux through $' num2str(Tl) '^\circ$C isotherm']);
% $$$ 
% $$$ leg = legend(legstr);
% $$$ % $$$ 'Eastern Equatorial Pacific', ...
% $$$ % $$$              'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
% $$$ % $$$              ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
% $$$ % $$$              ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
% $$$ % $$$              ['Total']);
% $$$ set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
% $$$ ylim([-0.8 0]);
% $$$ xlim([1 12]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'xtick',[1:12]);
% $$$ set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
% $$$ set(gca,'FontSize',25);
% $$$ grid on;
% $$$ % $$$ 
% $$$ % $$$ axes('Position',[0.74 0.2451 0.13 0.6799]);
% $$$ % $$$ plot([0 1],[min(MAll/1e15) min(MAll/1e15)],'-k','linewidth',2);
% $$$ % $$$ hold on;
% $$$ % $$$ plot([0 1],[max(MAll/1e15) max(MAll/1e15)],'-k','linewidth',2);
% $$$ % $$$ text(0.5,mean([min(MAll/1e15) max(MAll/1e15)]),sprintf('%3.2f ',max(MAll/1e15)-min(MAll/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ plot([1 2],[min(MEq/1e15) min(MEq/1e15)],'--k','linewidth',2);
% $$$ % $$$ plot([1 2],[max(MEq/1e15) max(MEq/1e15)],'--k','linewidth',2);
% $$$ % $$$ text(1.5,mean([min(MEq/1e15) max(MEq/1e15)]),sprintf('%3.2f ',max(MEq/1e15)-min(MEq/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ plot([2 3],[min(MEEP/1e15) min(MEEP/1e15)],'--r','linewidth',2);
% $$$ % $$$ plot([2 3],[max(MEEP/1e15) max(MEEP/1e15)],'--r','linewidth',2);
% $$$ % $$$ text(2.5,mean([min(MEEP/1e15) max(MEEP/1e15)]),sprintf('%3.2f ',max(MEEP/1e15)-min(MEEP/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','r','FontSize',25);
% $$$ % $$$ plot([3 4],[min(MN/1e15) min(MN/1e15)],':b','linewidth',2);
% $$$ % $$$ plot([3 4],[max(MN/1e15) max(MN/1e15)],':b','linewidth',2);
% $$$ % $$$ text(3.5,mean([min(MN/1e15) max(MN/1e15)]),sprintf('%3.2f ',max(MN/1e15)-min(MN/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color','b','FontSize',25);
% $$$ % $$$ plot([4 5],[min(MS/1e15) min(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$ % $$$                     0.5 0]);
% $$$ % $$$ plot([4 5],[max(MS/1e15) max(MS/1e15)],':','linewidth',2,'color',[0 ...
% $$$ % $$$                     0.5 0]);
% $$$ % $$$ text(4.5,mean([min(MS/1e15) max(MS/1e15)]),sprintf('%3.2f ',max(MS/1e15)- min(MS/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'color',[0 0.5 0],'FontSize',25);
% $$$ % $$$ plot([5 6],[min((MN+MS)/1e15) min((MN+MS)/1e15)],':k','linewidth',2);
% $$$ % $$$ plot([5 6],[max((MN+MS)/1e15) max((MN+MS)/1e15)],':k','linewidth',2);
% $$$ % $$$ text(5.5,mean([min((MN+MS)/1e15) max((MN+MS)/1e15)]),sprintf('%3.2f ',max((MN+MS)/1e15)-min((MN+MS)/1e15)),'VerticalAlignment','middle','HorizontalAlignment','Center','Rotation',270,'FontSize',25);
% $$$ % $$$ xlim([0 6]);
% $$$ % $$$ ylim([-0.8 0]);
% $$$ % $$$ set(gca,'xtick',[]);
% $$$ % $$$ set(gca,'ytick',[]);
% $$$ % $$$ box off;
% $$$ % $$$ grid off;
% $$$ % $$$ axis off;
% $$$ % $$$ title('Seasonal Range','FontSize',25);

    %%% Plot spatial pattern:

    try
        obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
        LAND = obj.SST(:,:,1);
    catch
        LAND = zeros(size(FlM(:,:,1)));
    end

    %If MOM01, fix NaN's in grid:
    if (strfind(model,'01'))
        lon = repmat(lon(:,500),[1 yL]);
        latv = nanmean(lat,1);
        lat = repmat(latv,[xL 1]);
    end

    [xL,yL] = size(lon);
    xvec = 1:1:xL;
    yvec = 1:1:yL;
    txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% $$$     months = {[1:12], ...
% $$$               [3], ...
% $$$               [7], ...
% $$$               [11]};
% $$$     labels = {'(a) Annual', ...
% $$$               '(b) March', ...
% $$$               '(c) July', ...
% $$$               '(d) November'};
    months = {[1:12]}; labels = {'Annual'};

    %Colormap and continents:
    sp = 1;
    clim = [-20 20];

    cCH = 0; % 0 = symmetric redblue
             % 1 = negative definite parula
             % 2 = negative parula with +ve's possible
    if (cCH==0)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = redblue(npts-3);
        for i=1:(npts-3)
            if (cmap(i,:) == 1.0)
                cmap(i,:) = [0.94 0.94 0.94];
            end
        end
    elseif (cCH>=1)
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        npts = length(cpts);
        cmap = parula(npts-3);
        cmap(end,:) = [0.97 0.97 0.8];
        cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
    end
    if (cCH == 2)
        buf = 6;
        clim = [clim(1) buf*sp];
        cpts = [-1e10 clim(1):sp:clim(2) 1e10];
        cmap(end+1,:) = cmap(end,:); % 1st positive bin
        cmap(end+buf-1,:) = [1 0.7 0.9]; % last pink bin
        for ii = 1:(buf-2)
            cmap(end-buf+1+ii,:) = cmap(end,:)*ii/(buf-1) + ...
                         cmap(end-buf+1,:)*(buf-1-ii)/(buf-1);
        end
    end        

    tmp = LAND;
    tmp(isnan(LAND)) = clim(1)-sp/2;
    tmp(~isnan(LAND)) = NaN;
    LAND = tmp;
    cmap(2:(end+1),:) = cmap;
    cmap(1,:) = [0 0 0];

    climn = [clim(1)-sp clim(2)];
    
%Mean of all months:
figure;
set(gcf,'Position',[3          59        1916         914]);
set(gcf,'defaulttextfontsize',20);
set(gcf,'defaultaxesfontsize',20);

poss = [0.1300    0.4553    0.7693    0.4697; ...
        0.1300    0.1389    0.2343    0.2680; ...
        0.3951    0.1389    0.2343    0.2680; ...
        0.6681    0.1389    0.2343    0.2680];
for i=1:length(months)
    if (i == 1)
        subplot(5,3,[1 9]);
    else
        subplot(5,3,[10 13]+(i-2));
    end
    X = lon(xvec,yvec);
    Y = lat(xvec,yvec);
    if (length(months{i})>1)
        tmp = FlM;
        tmp(isnan(tmp)) = 0.0;
        Z = monmean(tmp(:,:,months{i}),3,ndays(months{i}));
        Z(Z == 0) = NaN;
    else
        Z = FlM(:,:,months{i});
    end
    Z = Z(xvec,yvec);
    
    Z(Z<clim(1)) = clim(1);
    contourf(X,Y,Z,cpts,'linestyle','none');
    hold on;    
    contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
    caxis(climn);
    if (i==1)
        cb = colorbar;
        ylabel(cb,'Wm$^{-2}$');
        ylim(cb,clim);
    end
    hold on;
    % Plot regions:
    if (i==1)
    if (exist('regions'))
        xlims = get(gca,'xlim');
        for ii=2:length(regions)
            contour(lon,lat,regions{ii}{2},[0.5 0.5],'--','color',regions{ii}{3},'linewidth',1);
            % Add text label with total:
            flon = min(lon(regions{ii}{2}));
            flat = nanmean(lat(regions{ii}{2}));
            if (ii==2)
                flat = 7;
            elseif (flat>15)
                flat = 42;
            end
            text(flon,flat,sprintf('%3.2fPW (%3.0f%%)',monmean(Field(:,ii),1,ndays(1:12))/1e15,monmean(Field(:,ii),1,ndays(1:12))/monmean(Field(:,1),1,ndays(1:12))*100),'color',regions{ii}{3},'Interpreter','none','BackgroundColor','w','Margin',0.01);
        end
    end
    end
    xlabel('Longitude ($^\circ$E)');
    ylabel('Latitude ($^\circ$N)');
    set(gca,'xtick',[-270:30:60]);
    set(gca,'ytick',[-75:15:75]);
    set(gca,'Position',[poss(i,:)]);
    ylim([-45 45]);
% $$$     ylim([-55 55]);
% $$$     ylim([-75 75]);
    set(gca,'FontSize',15);
    colormap(cmap);
    title([strrep(strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),' jra55',''),'ryf8485',''),' may','') ...
           ' ' num2str(Tl) '$^\circ$C Numerical Mixing']);
end
end

% $$$ %%% Plot spatial pattern of net heat flux and SST:
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux = shfluxa/length(outputs);
% $$$ SST = SSTa/length(outputs);
% $$$ 
% $$$ %If MOM01, fix NaN's in grid:
% $$$ if (strfind(model,'01'))
% $$$     lon = repmat(lon(:,500),[1 yL]);
% $$$     latv = nanmean(lat,1);
% $$$     lat = repmat(latv,[xL 1]);
% $$$ end
% $$$ 
% $$$ %Sum of all positives:
% $$$ shfluxA = mean(shflux,3);
% $$$ shfluxP = nansum(shfluxA(shfluxA>0).*area(shfluxA>0))/1e15
% $$$ shfluxM = nansum(shfluxA(shfluxA<0).*area(shfluxA<0))/1e15
% $$$ 
% $$$ %Sum of mean flux above mean isotherm:
% $$$ SSTA = mean(SST,3);
% $$$ isot = 21.5;
% $$$ shfluxP = nansum(shfluxA(SSTA>=isot).*area(SSTA>=isot))/1e15
% $$$ shfluxM = nansum(shfluxA(SSTA<isot).*area(SSTA<isot))/1e15
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:3:xL;
% $$$ yvec = 1:3:yL;
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ 
% $$$ months = {[1:12], ...
% $$$                             [3], ...
% $$$                             [7], ...
% $$$                             [11]};
% $$$ 
% $$$ labels = {'(a) Annual', ...
% $$$           '(b) March', ...
% $$$           '(c) July', ...
% $$$           '(d) November'};
% $$$ 
% $$$ clim = [-200 200];
% $$$ sp = 20;
% $$$ 
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ poss = [0.1300    0.415  0.7693    0.56; ...
% $$$         0.1300    0.1    0.2343    0.2680; ...
% $$$         0.3951    0.1    0.2343    0.2680; ...
% $$$         0.6681    0.1    0.2343    0.2680];
% $$$ for i=1:length(months)
% $$$     if (i == 1)
% $$$         subplot(5,3,[1 9]);
% $$$     else
% $$$         subplot(5,3,[10 13]+(i-2));
% $$$     end
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
% $$$     Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$     Z = Z(xvec,yvec);
% $$$     Z2 = Z2(xvec,yvec);
% $$$     Z(Z==0) = NaN;
% $$$     
% $$$     % Map projection:
% $$$ % $$$     map = 1;
% $$$ % $$$     ax = worldmap('World');
% $$$ % $$$     setm(ax, 'Origin', [0 260 0])
% $$$ % $$$     contourfm(Y,X,Z,[-1e10 -500:20:500 1e10],'linestyle','none');
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$     hold on;
% $$$ % $$$     [c,h] = contourm(Y,X,Z2,[-3:2:21 25:2:35],'-k');
% $$$     [c,h] = contour(X,Y,Z2,[-3:2:35],'-k');
% $$$     clabel(c,h);        
% $$$ % $$$     hand = clabelm(c,h);        
% $$$ % $$$     set(hand,'FontSize',10,'BackgroundColor','none');
% $$$     if (i == 1)
% $$$         [c,h] = contour(X,Y,Z2,[21.5 21.5],'-k','linewidth',2);
% $$$         clabel(c,h)
% $$$ % $$$         [c,h] = contourm(Y,X,Z2,[23 23],'-k','linewidth',2);
% $$$ % $$$         hand = clabelm(c,h);
% $$$ % $$$         set(hand,'FontSize',10,'BackgroundColor','none');
% $$$     end
% $$$     caxis(clim);
% $$$     if (i==1)
% $$$         cb = colorbar;
% $$$         ylabel(cb,'Wm$^{-2}$');
% $$$     end
% $$$     ylim([-75 75]);
% $$$     if (i>1)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (i<=2)
% $$$         ylabel('Latitude ($^\circ$N)');
% $$$     end
% $$$     if (i>1)
% $$$ % $$$         text(-277,53,labels{i},'BackgroundColor','w');
% $$$         text(-277,64,labels{i},'BackgroundColor','w','margin',0.01);
% $$$         set(gca,'xtick',[-240:60:60]);
% $$$     else
% $$$ % $$$         text(-279,55,labels{i},'BackgroundColor','w');
% $$$         text(-279,69,labels{i},'BackgroundColor','w','margin',0.01);
% $$$         set(gca,'xtick',[-270:30:60]);
% $$$     end        
% $$$     set(gca,'ytick',[-60:30:60]);
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'color','k');
% $$$ 
% $$$ % $$$     setm(gca, 'MlabelParallel', 'south');
% $$$ % $$$     setm(gca, 'MlabelParallel', 'south');   
% $$$ % $$$     land = shaperead('landareas.shp', 'UseGeoCoords', true);
% $$$ % $$$     geoshow(land, 'FaceColor', [0 0 0])
% $$$ end 
% $$$ colormap(cmap);
% $$$ 
% $$$ 
% $$$ %% Plot SST difference and shflux difference plots:
% $$$ 
% $$$ % $$$ % Load MOM01 data:
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/';
% $$$ % $$$ model = 'MOM01';
% $$$ % $$$ outputs = [111 222];
% $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ ndays = diff(time_snap);
% $$$ % $$$ 
% $$$ % $$$ [xL,yL] = size(lon)
% $$$ % $$$ %If MOM01, fix NaN's in grid:
% $$$ % $$$ if (strfind(model,'01'))
% $$$ % $$$     lonv = lon(:,500);
% $$$ % $$$     latv = nanmean(lat,1);
% $$$ % $$$     [tmp ind] = min(abs(latv - 60))
% $$$ % $$$     [lon,lat] = ndgrid(lonv,latv(1:ind));
% $$$ % $$$ end
% $$$ % $$$ MOM01lon = lon;
% $$$ % $$$ MOM01lat = lat;
% $$$ % $$$ 
% $$$ % $$$ load([base model sprintf('_output%03d_SurfaceVars.mat', ...
% $$$ % $$$                          outputs(1))]);
% $$$ % $$$ MOM01shflux = shflux(:,1:ind,:);
% $$$ % $$$ MOM01SST = SST(:,1:ind,:);
% $$$ 
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ % $$$ model = 'MOM025_kb3seg';
% $$$ % $$$ outputs = [75:79];
% $$$ % $$$ 
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ % $$$ model = 'MOM025_kb1em5';
% $$$ % $$$ outputs = 94;
% $$$ % $$$ 
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = [8:12];
% $$$ % $$$ outputs = [14]
% $$$ % $$$ 
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em6/mat_data/';
% $$$ % $$$ model = 'MOM025_kb1em6';
% $$$ % $$$ outputs = 30;
% $$$ 
% $$$ % $$$ % Load ACCESS-OM2:
% $$$ % $$$ model = 'ACCESS-OM2_1deg_jra55_ryf8485_kds50_may';
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/mom/mat_data/';
% $$$ % $$$ outputs = 36;
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ 
% $$$ % Average all years:
% $$$ if (length(yrs)>1)
% $$$     shflux = reshape(shflux,[length(shflux(:,1,1)) length(shflux(1,:,1)) 12 nyrs]);
% $$$     shflux = mean(shflux(:,:,:,yrs),4);
% $$$     SST = reshape(SST,[length(SST(:,1,1)) length(SST(1,:,1)) 12 nyrs]);
% $$$     SST = mean(SST(:,:,:,yrs),4);
% $$$     ndays = ndays(1:12);
% $$$ else
% $$$     shfluxa = shflux;
% $$$     SSTa = SST;
% $$$     SSTa = SSTa(:,:,13:24);
% $$$     ndays = ndays(1:12);
% $$$     for i=2:length(outputs)
% $$$         load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$         shfluxa = shfluxa+shflux;
% $$$         SSTa = SSTa+SST;
% $$$     end
% $$$     shflux = shfluxa/length(outputs);
% $$$     SST = SSTa/length(outputs);
% $$$ end
% $$$ 
% $$$ if (max(max(max(SST)))>100)
% $$$     SST = SST-273.15;
% $$$ end
% $$$ 
% $$$ % Calculate bias from first run:
% $$$ % $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ % $$$ ndays = diff(time_snap);
% $$$ if (rr >=2)
% $$$     SSTcur = SST;
% $$$     load([base RUNS{1}{1} sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$     SST = reshape(SST,[length(SST(:,1,1)) length(SST(1,:,1)) 12 nyrs]);
% $$$     SST = mean(SST(:,:,:,yrs),4);
% $$$     if (max(max(max(SST)))>100)
% $$$     SST = SST-273.15;
% $$$     end
% $$$     SSTbias = monmean(SSTcur,3,ndays) - monmean(SST,3,ndays);
% $$$ else
% $$$ % WOA13 SST:
% $$$ WOAname = '/srv/ccrc/data03/z3500785/WOA13/woa13_decav_t00_04v2.nc';
% $$$ WOASST = ncread(WOAname,'t_an',[1 1 1 1],[1440 720 1 1]);
% $$$ [WOAlon,WOAlat] = ndgrid(ncread(WOAname,'lon'),ncread(WOAname,'lat'));
% $$$ 
% $$$ %Shift longitudes:
% $$$ [tmp ind] = min(abs(WOAlon(:,1)-80));
% $$$ WOASST = cat(1,WOASST(ind+1:end,:),WOASST(1:ind,:));
% $$$ WOAlon = cat(1,WOAlon(ind+1:end,:)-360,WOAlon(1:ind,:));
% $$$ WOAlat = cat(1,WOAlat(ind+1:end,:),WOAlat(1:ind,:));
% $$$ 
% $$$ % Calculate bias from WOA:
% $$$ SSTbias = monmean(SST,3,ndays)-interp2(WOAlon',WOAlat',WOASST',lon,lat,'linear');
% $$$ % $$$ for i=1:12
% $$$ % $$$     SST(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01SST(:,:,i)',lon,lat,'linear')-SST(:,:,i);
% $$$ % $$$     shflux(:,:,i) = interp2(MOM01lon',MOM01lat',MOM01shflux(:,:,i)',lon,lat,'linear')-shflux(:,:,i);
% $$$ % $$$ end
% $$$ end
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:2:xL;
% $$$ yvec = 1:2:yL;
% $$$ 
% $$$ % $$$ figure;
% $$$ % $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ % $$$ set(gcf,'defaulttextfontsize',15);
% $$$ % $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ % $$$ poss = [0.15  0.56 0.3 0.34; ...
% $$$ % $$$         0.48  0.56 0.3 0.34; ...
% $$$ % $$$         0.15  0.1 0.3 0.34; ...
% $$$ % $$$         0.48  0.1 0.3 0.34];
% $$$ clvls = [-1e10 -1:0.05:1 1e10];
% $$$ subplot(2,3,rr);
% $$$ contourf(lon(xvec,yvec),lat(xvec,yvec),SSTbias(xvec,yvec),clvls,'linestyle', ...
% $$$          'none');
% $$$ hold on;
% $$$ % $$$ contour(lon,lat,SSTbias,[-3 -2 -1 1 2 3],'-k');
% $$$ if (rr == 1)
% $$$     contour(lon,lat,SSTbias,[-1:0.5:1],'-k');
% $$$ end
% $$$ set(gca,'color','k');
% $$$ % $$$ title('$\kappa_B=10^{-5}$ - WOA13 SST Year 2 ($^\circ$C)');
% $$$ if (rr == 1)
% $$$     title([strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''),' may','') ' - WOA SST ($^\circ$C)']);
% $$$ else
% $$$     title([strrep(strrep(strrep(strrep(RUNS{rr}{1},'_',' '),'ACCESS-OM2 ','AOM'),'deg jra55',''),' may','') ' - kds50 SST ($^\circ$C)']);
% $$$ end    
% $$$ if (rr == 1)
% $$$     caxis([-1 1]);
% $$$ else
% $$$     caxis([-0.5 0.5]);
% $$$ end
% $$$ colorbar;
% $$$ colormap(redblue(length(clvls)-1));
% $$$ set(gca,'FontSize',15);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ % $$$ set(gca,'Position',poss(rr,:));
% $$$ end
% $$$ 
% $$$ 
% $$$ %
% $$$ for i=1:length(months)
% $$$     if (i == 1)
% $$$         subplot(5,3,[1 9]);
% $$$     else
% $$$         subplot(5,3,[10 13]+(i-2));
% $$$     end
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     Z = monmean(shflux(:,:,months{i}),3,ndays(months{i}));
% $$$     Z2 = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$     Z = Z(xvec,yvec);
% $$$     Z2 = Z2(xvec,yvec);
% $$$ % $$$     contourf(X,Y,Z2,[-1e10 -5:0.25:5 1e10],'linestyle','none');
% $$$     contourf(X,Y,Z,[-1e10 -500:10:500 1e10],'linestyle','none');
% $$$     hold on;
% $$$ % $$$     quiver(lon(xvec2,yvec2),lat(xvec2,yvec2),tau_x(xvec2,yvec2),tau_y(xvec2,yvec2),3,'-k');
% $$$ % $$$     if (i==1)
% $$$ % $$$         [c,h] = contour(X,Y,Z,[-100:10:-10],'--k');
% $$$ % $$$         [c,h] = contour(X,Y,Z,[10:10:100],'-k');
% $$$ % $$$     end
% $$$ % $$$     else
% $$$ % $$$         [c,h] = contour(X,Y,Z2,[-3:4:35],'-k');
% $$$ % $$$     end
% $$$ % $$$     clabel(c,h);
% $$$     caxis([-100 100]);
% $$$     if (i==1)
% $$$         cb = colorbar;
% $$$ % $$$         ylabel(cb,'$^\circ$C');
% $$$         ylabel(cb,'Wm$^{-2}$');
% $$$     end
% $$$     ylim([-75 60]);
% $$$     if (i>1)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (i<=2)
% $$$         ylabel('Latitude ($^\circ$N)');
% $$$     end
% $$$     if (i>1)
% $$$         text(-276,53,labels{i},'BackgroundColor','w');
% $$$     else
% $$$         text(-278,55,labels{i},'BackgroundColor','w');
% $$$     end        
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'color','k');
% $$$ end 
% $$$ colormap(redblue);
% $$$ 
% $$$ %%% Outcrop area plot:
% $$$ 
% $$$ Aout = zeros(size(Te));
% $$$ for i=1:length(Te)
% $$$     i
% $$$     Aout(i) = nansum(area(mean(SST,3)>Te(i)));
% $$$ end
% $$$ figure;
% $$$ plot(Aout/1e6,Te,'-k','linewidth',2);
% $$$ xlabel('Outcrop area (km$^2$)');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ plot(-diff(Aout,[],1)/dT/1e6,T,'-k','linewidth',2);
% $$$ xlabel('Outcrop area (km$^2$/$^\circ$C)');
% $$$ ylabel('Temperature ($^\circ$C)');
% $$$ 
% $$$ 
% $$$ %Calculate fractions of surface heat flux regionally:
% $$$ 
% $$$ %Region choice:
% $$$ LATsplit = 10;
% $$$ reg = [-150 -90 -5 5];
% $$$ EEPinds = lat > reg(3) & lat < reg(4) & lon > reg(1) & lon < reg(2);
% $$$ 
% $$$ Eqinds = abs(lat)<=LATsplit;
% $$$ Ninds = lat>LATsplit;
% $$$ Sinds = lat<-LATsplit;
% $$$ 
% $$$ for i=1:12
% $$$     SHFtmp = shflux(:,:,i);
% $$$     SHFEq(i) = nansum(nansum(area(Eqinds).*SHFtmp(Eqinds),1),2);
% $$$     SHFN(i) = nansum(nansum(area(Ninds).*SHFtmp(Ninds),1),2);
% $$$     SHFS(i) = nansum(nansum(area(Sinds).*SHFtmp(Sinds),1),2);
% $$$     SHFEEP(i) = nansum(nansum(area(EEPinds).*SHFtmp(EEPinds),1),2);
% $$$ 
% $$$     Atmp = area;
% $$$     Atmp(isnan(SHFtmp)) = NaN;
% $$$     AREAtotal(i) = nansum(nansum(Atmp));
% $$$     AREAEEP(i) = nansum(nansum(Atmp(EEPinds)));
% $$$     AREAEq(i) = nansum(nansum(Atmp(Eqinds)));
% $$$     AREAN(i) = nansum(nansum(Atmp(Ninds)));
% $$$     AREAS(i) = nansum(nansum(Atmp(Sinds)));
% $$$ end
% $$$ SHFAll = SHFEq+SHFN+SHFS;
% $$$ 
% $$$ %Display output in terminal for table:
% $$$ str = {['Area fractions model ' model] ; ...
% $$$ sprintf(' NH area = %3.2f',monmean(AREAN,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' SH area = %3.2f',monmean(AREAS,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' Eq area = %3.2f',monmean(AREAEq,2,ndays)/monmean(AREAtotal,2,ndays)) ; ...
% $$$ sprintf(' EEP area = %3.2f',monmean(AREAEEP,2,ndays)/monmean(AREAtotal,2,ndays))}
% $$$ str = {['Annual totals model ' model]  ; ...
% $$$ sprintf(' Total = %3.2f',monmean(SHFAll,2,ndays)/1e15) ; ...
% $$$ sprintf(' NH = %3.2fPW (%3.0f)',monmean(SHFN,2,ndays)/1e15,monmean(SHFN,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' SH = %3.2fPW (%3.0f)',monmean(SHFS,2,ndays)/1e15,monmean(SHFS,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' Eq = %3.2fPW (%3.0f)',monmean(SHFEq,2,ndays)/1e15,monmean(SHFEq,2,ndays)/monmean(SHFAll,2,ndays)*100) ; ...
% $$$ sprintf(' EEP = %3.2fPW (%3.0f)',monmean(SHFEEP,2,ndays)/1e15,monmean(SHFEEP,2,ndays)/monmean(SHFAll,2,ndays)*100)}
% $$$ 
% $$$ % Plot regional heat fluxes as a function of season:
% $$$ figure;
% $$$ set(gcf,'Position',get(0,'ScreenSize'));
% $$$ set(gcf,'defaultlinelinewidth',2);
% $$$ 
% $$$ % $$$ subplot(3,1,1);
% $$$ plot(1:12,SHFEEP/1e15,'--r','LineWidth',4);
% $$$ hold on;
% $$$ plot(1:12,SHFEq/1e15,'--k','LineWidth',4);
% $$$ plot(1:12,SHFN/1e15,':b','LineWidth',4);
% $$$ plot(1:12,SHFS/1e15,':','color',[0 0.5 0],'LineWidth',4);
% $$$ plot(1:12,(SHFN+SHFS)/1e15,':k','LineWidth',4);
% $$$ plot(1:12,SHFAll/1e15,'-k','LineWidth',4);
% $$$ xlabel('Month');
% $$$ ylabel(['PW']);
% $$$ title(['Surface Heat Flux']);
% $$$ leg = legend('Eastern Equatorial Pacific', ...
% $$$              'Equatorial',['Northern Hemisphere $>' num2str(LATsplit) '^\circ$N'], ...
% $$$              ['Southern Hemisphere $<' num2str(LATsplit) '^\circ$S'], ...
% $$$              ['Outside Equatorial $>\pm' num2str(LATsplit) '^\circ$'], ...
% $$$              ['Total']);
% $$$ set(leg,'Position',[0.1481    0.2656    0.2314    0.2340]);
% $$$ ylim([-20 20]);
% $$$ xlim([1 12]);
% $$$ set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% $$$ set(gca,'xtick',[1:12]);
% $$$ set(gca,'Position',[0.1300    0.2451    0.5949    0.6799]);
% $$$ set(gca,'FontSize',25);
% $$$ grid on;
% $$$ 
% $$$ 
% $$$ 
% $$$ %%% Plot Seasonal cycle of wind stress and SST
% $$$ 
% $$$ % Load Base Variables:
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ model = 'MOM025_kb3seg';
% $$$ outputs = [75:79];
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ % $$$ model = 'MOM025';
% $$$ % $$$ outputs = [8:12];
% $$$ % $$$ model = 'MOM01';
% $$$ % $$$ outputs = [333];
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ SSTa = SST;
% $$$ tauxa = taux;
% $$$ tauya = tauy;
% $$$ for i=2:length(outputs)
% $$$     i
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     SSTa = SSTa+SST;
% $$$     tauxa = tauxa+taux;
% $$$     tauya = tauya+tauy;
% $$$ end
% $$$ SST = SSTa/length(outputs);
% $$$ taux = tauxa/length(outputs);
% $$$ tauy = tauya/length(outputs);
% $$$ 
% $$$ ylims = [-30 30];
% $$$ xlims = [-260 -60];
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ [tmp ln1] = min(abs(lon(:,1)-xlims(1)));
% $$$ [tmp ln2] = min(abs(lon(:,1)-xlims(2)));
% $$$ [tmp lt1] = min(abs(lat(1,:)-ylims(1)));
% $$$ [tmp lt2] = min(abs(lat(1,:)-ylims(2)));
% $$$ xvec = (ln1-2):1:(ln2-2);
% $$$ yvec = (lt1-2):1:(lt2-2);
% $$$ [xL,yL] = size(lonu);
% $$$ [tmp ln1] = min(abs(lonu(:,1)-xlims(1)));
% $$$ [tmp ln2] = min(abs(lonu(:,1)-xlims(2)));
% $$$ [tmp lt1] = min(abs(latu(1,:)-ylims(1)));
% $$$ [tmp lt2] = min(abs(latu(1,:)-ylims(2)));
% $$$ xvec2 = (ln1-2):15:(ln2-2);
% $$$ yvec2 = (lt1-2):15:(lt2-2);
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ 
% $$$ months = {[1],[3],[5],[7],[9],[11]};
% $$$ % $$$           [1:12], ...
% $$$ % $$$           [3], ...
% $$$ % $$$           [7], ...
% $$$ % $$$           [11]};
% $$$ 
% $$$ labels = {'Jan', ...
% $$$           'Mar', ...
% $$$           'May', ...
% $$$           'Jul', ...
% $$$           'Sep', ...
% $$$           'Nov'};
% $$$ % $$$ 'Annual', ...
% $$$ % $$$           'March', ...
% $$$ % $$$           'July', ...
% $$$ % $$$           'November'};
% $$$ 
% $$$ % $$$ %Colormap:
% $$$ % $$$ clim = [0 0.1];
% $$$ % $$$ sp = 0.01;
% $$$ % $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ % $$$ npts = length(cpts)
% $$$ % $$$ cmap = parula(npts-3);
% $$$ % $$$ cmap = parula(npts-3);
% $$$ % $$$ % $$$ cmap(end,:) = [0.97 0.97 0.8];
% $$$ % $$$ cmap(end,:) = [1 1 1];
% $$$ % $$$ cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ % $$$ cmap = flipud(cmap);
% $$$ 
% $$$ 
% $$$ clim = [20 30];
% $$$ sp = 1;
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = colormap(redblue(npts-3));
% $$$ cmap = colormap(flipud(lbmap(npts-3,'RedBlue')));
% $$$ 
% $$$ figure;
% $$$ %set(gcf,'Position',[3          59        1916         914]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ set(gcf,'Position',[322          58        1247         945]);
% $$$ % $$$ poss = [0.1300    0.4553    0.7693    0.4697; ...
% $$$ % $$$         0.1300    0.1389    0.2343    0.2680; ...
% $$$ % $$$         0.3951    0.1389    0.2343    0.2680; ...
% $$$ % $$$         0.6681    0.1389    0.2343    0.2680];
% $$$ poss = [0.0867    0.6910    0.4    0.26; ...
% $$$         0.5221    0.6910    0.4    0.26; ...
% $$$         0.0867    0.3926    0.4    0.26; ...
% $$$         0.5221    0.3926    0.4    0.26; ...
% $$$         0.0867    0.0899    0.4    0.26; ...
% $$$         0.5221    0.0899    0.4    0.26];
% $$$ for i=1:length(months)
% $$$ % $$$     if (i == 1)
% $$$ % $$$         subplot(5,3,[1 9]);
% $$$ % $$$     else
% $$$ % $$$         subplot(5,3,[10 13]+(i-2));
% $$$ % $$$     end
% $$$     subplot(3,2,i);
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     Z = monmean(SST(:,:,months{i}),3,ndays(months{i}));
% $$$     Z2 = monmean(taux(:,:,months{i}),3,ndays(months{i}));
% $$$     Z3 = monmean(tauy(:,:,months{i}),3,ndays(months{i}));
% $$$     Z = Z(xvec,yvec);
% $$$     Z(Z==0) = NaN;
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$ % $$$     contourf(lonu(xvec,yvec),latu(xvec,yvec),sqrt(Z2(xvec,yvec).^2+ ...
% $$$ % $$$                                                   Z3(xvec,yvec).^2), ...
% $$$ % $$$              [-1e10 0:0.01:0.5 1e10],'linestyle','none');
% $$$     hold on;
% $$$ % $$$     [c,h] = contour(X,Y,Z,[-3:2:35],'-k');
% $$$ % $$$     clabel(c,h);
% $$$     quiver(lonu(xvec2,yvec2),latu(xvec2,yvec2),Z2(xvec2,yvec2),Z3(xvec2,yvec2),0.75,'-k');
% $$$     caxis(clim);
% $$$     xlim(xlims);
% $$$     ylim(ylims);
% $$$     hold on;
% $$$     plot([-150 -90 -90 -150 -150],[-5 -5 5 5 -5],'-','color',[0 0.5 ...
% $$$                         0],'linewidth',2);
% $$$     if (i>=5)
% $$$         xlabel('Longitude ($^\circ$E)');
% $$$     else
% $$$         set(gca,'xticklabel',[]);
% $$$     end
% $$$     if (i==1 | i == 3 | i ==5 )
% $$$         ylabel('Latitude ($^\circ$N)');
% $$$     else
% $$$         set(gca,'yticklabel',[]);
% $$$     end
% $$$ % $$$     if (i>1)
% $$$         text(-257,-25,labels{i},'BackgroundColor','w');
% $$$ % $$$     else
% $$$ % $$$         text(-278,35,labels{i},'BackgroundColor','w');
% $$$ % $$$     end        
% $$$     if (i==2 | i == 4 | i ==6)
% $$$         cb = colorbar;
% $$$         ylabel(cb,'$^\circ$C');
% $$$     end
% $$$     set(gca,'Position',[poss(i,:)]);
% $$$     set(gca,'color','k');
% $$$     LabelAxes(gca,i,20,0.008,0.93);
% $$$ end 
% $$$ colormap(cmap);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%% Components of mixing plot:
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(gcf,'defaulttextfontsize',17);
% $$$ set(gcf,'defaultaxesfontsize',17);
% $$$ 
% $$$ poss = [0.15  0.56 0.3 0.34; ...
% $$$         0.48  0.56 0.3 0.34; ...
% $$$         0.15  0.15 0.3 0.34; ...
% $$$         0.48  0.15 0.3 0.34];
% $$$ 
% $$$ obj = matfile([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ LAND = obj.SST(:,:,1);
% $$$ 
% $$$ clim = [-40 0];
% $$$ sp = 2;
% $$$ doWMT = 0; % plot WMT instead of flux
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ if (doWMT)
% $$$     cmap = flipud(lbmap(npts-3,'RedBlue'));
% $$$     cmap = redblue(npts-3);
% $$$     for i=1:(npts-3)
% $$$         if (cmap(i,:) == 1.0)
% $$$             cmap(i,:) = [0.94 0.94 0.94];
% $$$         end
% $$$     end
% $$$ else
% $$$     cmap = parula(npts-3);
% $$$     cmap = parula(npts-3);
% $$$     cmap(end,:) = [0.97 0.97 0.8];
% $$$     cmap(end-1,:) = (cmap(end-1,:)+cmap(end,:))/2;
% $$$ end
% $$$ tmp = LAND;
% $$$ tmp(isnan(LAND)) = clim(1)-sp/2;
% $$$ tmp(~isnan(LAND)) = NaN;
% $$$ LAND = tmp;
% $$$ cmap(2:(end+1),:) = cmap;
% $$$ cmap(1,:) = [0 0 0];
% $$$ climn = [clim(1)-sp clim(2)];
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:1:xL;
% $$$ yvec = 1:1:yL;
% $$$ txtmonth = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
% $$$ months = [1:12];
% $$$ labels = {'(a) Background', ...
% $$$           '(b) Shear Instability', ...
% $$$           '(c) KPP Boundary Layer', ...
% $$$           '(d) Internal Wave'};
% $$$ VARS = {'FlMkppiw','FlMkppish','FlMkppbl','FlMwave'};
% $$$ 
% $$$ for ii=1:4
% $$$     subplot(2,2,ii);
% $$$     
% $$$ VAR = VARS{ii};
% $$$ TYPE = 'VertInt';
% $$$ Tl = 21.5;
% $$$ name = [base model sprintf('_output%03d',outputs(1)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$ eval(['load(name,''' VAR ''');']);
% $$$ eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$ eval([VAR 'a = ' VAR ';']);
% $$$ for i=2:length(outputs)
% $$$     name = [base model sprintf('_output%03d',outputs(i)) '_' TYPE '_T' strrep(num2str(Tl),'.','p') 'C.mat']
% $$$     eval(['load(name,''' VAR ''');']);
% $$$     eval([VAR '(isnan(' VAR ')) = 0.0;']);
% $$$     eval([VAR 'a = ' VAR 'a + ' VAR ';']);
% $$$ end
% $$$ eval([VAR ' = ' VAR 'a/length(outputs);']);
% $$$ eval([VAR '(' VAR '==0) = NaN;']);
% $$$ eval(['FlM = ' VAR ';']);
% $$$ 
% $$$     X = lon(xvec,yvec);
% $$$     Y = lat(xvec,yvec);
% $$$     if (length(months)>1)
% $$$         tmp = FlM;
% $$$         tmp(isnan(tmp)) = 0.0;
% $$$         Z = monmean(tmp(:,:,months),3,ndays(months));
% $$$         Z(Z == 0) = NaN;
% $$$     else
% $$$         Z = FlM(:,:,months);
% $$$     end
% $$$     Z = Z(xvec,yvec);
% $$$     if (doWMT)
% $$$         Z = Z*86400;
% $$$     end
% $$$     
% $$$     Z(Z<clim(1)) = clim(1);
% $$$     contourf(X,Y,Z,cpts,'linestyle','none');
% $$$     hold on;    
% $$$     contourf(X,Y,LAND(xvec,yvec),[clim(1)-sp clim(1)],'linestyle','none');
% $$$     caxis(climn);
% $$$     hold on;
% $$$ 
% $$$     if (ii>=3)
% $$$     xlabel('Longitude ($^\circ$E)');
% $$$     end
% $$$     if (ii==1 | ii==3)
% $$$     ylabel('Latitude ($^\circ$N)');
% $$$     end
% $$$     if (ii==2 | ii==4)
% $$$         cb = colorbar;
% $$$         if (~doWMT)
% $$$             ylabel(cb,'Wm$^{-2}$');
% $$$         else
% $$$             ylabel(cb,'m/day');
% $$$         end            
% $$$         ylim(cb,clim);
% $$$     end
% $$$     text(-278,-40.5,labels{ii},'BackgroundColor','w','Margin',0.5);
% $$$ % $$$     title(labels{ii});
% $$$     ylim([-45 45]);
% $$$     set(gca,'ytick',[-45:15:45]);
% $$$     set(gca,'xtick',[-240:60:60]);
% $$$     set(gca,'Position',poss(ii,:));
% $$$ end 
% $$$ colormap(cmap);
% $$$ 
% $$$ 
% $$$ %%% Plot shflux difference between runs with SST contours:
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb1em5/mat_data/';
% $$$ model = 'MOM025_kb1em5';
% $$$ outputs = 94;
% $$$ % $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag_kb3seg/mat_data/';
% $$$ % $$$ model = 'MOM025_kb3seg';
% $$$ % $$$ outputs = [75:79];
% $$$ 
% $$$ load([base model sprintf('_output%03d_BaseVars.mat',outputs(1))]);
% $$$ ndays = diff(time_snap);
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux1 = monmean(shfluxa/length(outputs),3,ndays);
% $$$ SST1 = monmean(SSTa/length(outputs),3,ndays);
% $$$ clear shfluxa SSTa;
% $$$ 
% $$$ base = '/srv/ccrc/data03/z3500785/MOM_HeatDiag/mat_data/';
% $$$ model = 'MOM025';
% $$$ outputs = [8:12];
% $$$ % $$$ outputs = [12]
% $$$ 
% $$$ % Load Variable and calculate mean:
% $$$ load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(1))]);
% $$$ shfluxa = shflux;
% $$$ SSTa = SST;
% $$$ for i=2:length(outputs)
% $$$     load([base model sprintf('_output%03d_SurfaceVars.mat',outputs(i))]);
% $$$     shfluxa = shfluxa+shflux;
% $$$     SSTa = SSTa+SST;
% $$$ end
% $$$ shflux2 = monmean(shfluxa/length(outputs),3,ndays);
% $$$ SST2 = monmean(SSTa/length(outputs),3,ndays);
% $$$ clear shfluxa SSTa;
% $$$ 
% $$$ [xL,yL] = size(lon);
% $$$ xvec = 1:1:xL;
% $$$ yvec = 291:1:706;
% $$$ 
% $$$ clim = [-25 25];
% $$$ sp = 2.5;
% $$$ 
% $$$ cpts = [-1e10 clim(1):sp:clim(2) 1e10];
% $$$ npts = length(cpts)
% $$$ cmap = redblue(npts-3);
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[3          59        1916         914]);
% $$$ set(0,'defaulttextfontsize',20);
% $$$ set(0,'defaultaxesfontsize',20);
% $$$ poss = [0.1300    0.415  0.7693    0.56; ...
% $$$         0.1300    0.1    0.2343    0.2680; ...
% $$$         0.3951    0.1    0.2343    0.2680; ...
% $$$         0.6681    0.1    0.2343    0.2680];
% $$$ subplot(5,3,[1 9]);
% $$$ X = lon(xvec,yvec);
% $$$ Y = lat(xvec,yvec);
% $$$ Z = shflux1(xvec,yvec)-shflux2(xvec,yvec);
% $$$ contourf(X,Y,Z,cpts,'linestyle','none');
% $$$ hold on;
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec)-SST2(xvec,yvec),[-1.5:0.125:-0.125],'--k');
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec)-SST2(xvec,yvec),[0.125:0.125:1.5],'-k');
% $$$ [c,h] = contour(X,Y,SST1(xvec,yvec),[21.5 21.5],'-k','linewidth',3);
% $$$ % $$$ [c,h] = contour(X,Y,SST2(xvec,yvec),[21.5 21.5],'--','color',[0 0.5 ...
% $$$ % $$$                     0],'linewidth',2);
% $$$ caxis(clim);
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Wm$^{-2}$');
% $$$ ylim([-45 45]);
% $$$ set(gca,'ytick',[-45:15:45]);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ set(gca,'xtick',[-240:60:60]);
% $$$ set(gca,'Position',[poss(1,:)]);
% $$$ set(gca,'color','k');
% $$$ colormap(cmap);
