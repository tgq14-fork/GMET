clear
%tgq: calculate the ratio (prcp) or difference (temperature) between
%station data and ensemble mean, and save those file in

EnsPath='/home/gut428/GMET/eCAI_EMDNA/output/NA_reg_ens';
StnClimoPath = '/home/gut428/GMET/StnInput_climatology'; % store climatological mean of prcp/temp
StnAnomPath = '/home/gut428/GMET/StnInput_dailyanomaly'; % store the anomaly (ratio or difference) or prcp/temp
StnDailyPath = '/home/gut428/GMET/StnInput_daily'; % original station time series
StnGridInfoPath = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo'; % grid information
StnInfoFile='/home/gut428/EMDNA/1StationInput/StnInfo.mat';

year=2018:2018;
month=1:1;

% read station list
load(StnInfoFile,'StnInfo');
stnID=StnInfo.ID;
nsta = length(stnID);

% read grid information of the study area
infile_grid = sprintf('%s/gridinfo_whole.nc',StnGridInfoPath);
lat = ncread(infile_grid,'latitude');
lon = ncread(infile_grid,'longitude');
lat1d = lat(:,1);
lon1d = lon(1,:);

% date of daily
datedaiy=datenum(1979,1,1):datenum(2018,12,31);
datedaiy=datestr(datedaiy,'yyyymm');
datedaiy=num2cell(datedaiy,2);
datedaiy=str2double(datedaiy);

for i = 1:nsta
    fprintf('Processing %d--%d\n',i,nsta);
    %tgq: for each station, find the nearest grid cells
    %extract the mean of all ensembles, and calculate the ratio (pcp_a)
    %between station prcp and ensemble mean prcp
    infile_stndaily = sprintf('%s/%s.nc',StnDailyPath,stnID{i});
    outfile_stnanom = sprintf('%s/%s.nc',StnAnomPath,stnID{i});
    if ~exist(outfile_stnanom,'file')
        copyfile(infile_stndaily,outfile_stnanom); % just copy to avoid creating new variables. data will be over-written
    end
    
    temp=ncinfo(infile_stndaily);
    invar=cell(length(temp.Variables),1);
    for j=1:length(temp.Variables)
        invar{j}=temp.Variables(j).Name;
    end
    
    % read stn lat/lon and its x/y
    slat = ncread(infile_stndaily,'latitude');
    slon = ncread(infile_stndaily,'longitude');
    [~, y] = min(abs(lat1d-slat));
    [~, x] = min(abs(lon1d-slon));
    
    for yy=1:length(year)
        EnsPathyy = [EnsPath,'/',num2str(year(yy))];
        for mm=1:length(month)
            indexdaily= find(datedaiy==year(yy)*100+mm);
            locstart=indexdaily(1);
            loclen=length(indexdaily);
            
            %tgq: climate mean of prcp tmean trange using all ensemble members
            %based on results from ensemble using TIME_MODE=climo
            infile_mean = sprintf('%s/ENS_CLIMO_%02d.nc',EnsPathyy,mm);
            pcp_clim = ncread(infile_mean,'pcp');
            tmean_clim = ncread(infile_mean,'t_mean');
            trange_clim = ncread(infile_mean,'t_range');
            tmax_clim = tmean_clim + trange_clim/2;
            tmin_clim = tmean_clim - trange_clim/2;

            % process prcp
            if ismember('prcp',invar)
                pcp_d = ncread(infile_stndaily,'prcp',locstart,loclen);
                pcp_a = pcp_d./pcp_clim(x,y);
                ncwrite(outfile_stnanom,'prcp',pcp_a,locstart,loclen);
            end
            
            %tgq: for temperature, calculate the difference between station
            %temperature and ensemble mean temperature.
            if ismember('tmin',invar)
                tmin_d = ncread(infile_stndaily,'tmin',locstart,loclen);
                tmin_a = tmin_d - tmin_clim(x,y);
                ncwrite(outfile_stnanom,'tmin',tmin_a,locstart,loclen);
            end
            
            if ismember('tmax',invar)
                tmax_d = ncread(infile_stndaily,'tmax',locstart,loclen);
                tmax_a = tmax_d - tmax_clim(x,y);
                ncwrite(outfile_stnanom,'tmax',tmax_a,locstart,loclen);
            end
        end
        
    end
end