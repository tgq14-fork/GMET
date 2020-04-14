clc;clear;
%tgq: calculate the ratio (prcp) or difference (temperature) between
%station data and ensemble mean, and save those file in

EnsPath='/home/gut428/GMET/eCAI_EMDNA/CV_output_1/NA_reg_ens';
StnAnomPath = '/home/gut428/GMET/StnInput_dailyanomaly'; % store the anomaly (ratio or difference) or prcp/temp
StnDailyPath = '/home/gut428/GMET/StnInput_daily'; % original station time series
StnGridInfoPath = '/home/gut428/GMET/eCAI_EMDNA/StnGridInfo'; % grid information
StnInfoFile='/home/gut428/EMDNA/1StationInput/StnInfo.mat';

year=2018:2018;
month=1:12; % month must be continuous

% read station list
load(StnInfoFile,'StnInfo');
stnID=StnInfo.ID;
stnLLE=StnInfo.lle;
nsta = length(stnID);

% read grid information of the study area
infile_grid = sprintf('%s/gridinfo_whole.nc',StnGridInfoPath);
lat2d = ncread(infile_grid,'latitude');
lon2d = ncread(infile_grid,'longitude');
elev = ncread(infile_grid,'elev');
lat1d = lat2d(1,:);
lon1d = lon2d(:,1);

% station index in the grid
indy=zeros(nsta,1);
indx=zeros(nsta,1);
for i=1:nsta
    [~, indy(i)] = min(abs(lat1d-stnLLE(i,1)));
    [~, indx(i)] = min(abs(lon1d-stnLLE(i,2)));
end
% some stations along boundary may have nan values, their indx indy
% need adjustment.
indxy=sub2ind(size(lat2d),indx,indy);
indtemp=find(isnan(elev(indxy)));
if ~isempty(indtemp)
    % usually adding one row or column can solve this problem
    indx(indtemp)=indx(indtemp)+1;
    indxy=sub2ind(size(lat2d),indx,indy);
end


% date of daily
datedaiy=datenum(1979,1,1):datenum(2018,12,31);
datedaiy=datestr(datedaiy,'yyyymm');
datedaiy=num2cell(datedaiy,2);
datedaiy=str2double(datedaiy);

for yy=1:length(year)
    EnsPathyy = [EnsPath,'/',num2str(year(yy))];
    % find the index of the months
    indexdaily= find(ismember( datedaiy, year(yy)*100+month));
    locstart=indexdaily(1);
    loclen=length(indexdaily);
    
    % read ensemble climo data
    pcp_clim=nan*zeros(loclen,nsta);
    tmin_clim=nan*zeros(loclen,nsta);
    tmax_clim=nan*zeros(loclen,nsta);
    flag=1;
    for mm=1:length(month)
        infile_mean = sprintf('%s/ENS_CLIMO_%02d.nc',EnsPathyy,mm);
        pcp_clim0 = ncread(infile_mean,'pcp'); pcp_clim0(pcp_clim0==-999)=nan;
        tmean_clim0 = ncread(infile_mean,'t_mean'); tmean_clim0(tmean_clim0==-999)=nan;
        trange_clim0 = ncread(infile_mean,'t_range'); trange_clim0(trange_clim0==-999)=nan;
        tmax_clim0 = tmean_clim0 + trange_clim0/2;
        tmin_clim0 = tmean_clim0 - trange_clim0/2;
        
        
        ndays=sum(datedaiy==year(yy)*100+mm);
        pcp_clim(flag:flag+ndays-1,:)=repmat(pcp_clim0(indxy)',ndays,1);
        tmin_clim(flag:flag+ndays-1,:)=repmat(tmin_clim0(indxy)',ndays,1);
        tmax_clim(flag:flag+ndays-1,:)=repmat(tmax_clim0(indxy)',ndays,1);
        flag=flag+ndays;
    end
    
    % for each station, calculate anomaly and save
    parfor i = 1:nsta
        fprintf('Processing %d--%d\n',i,nsta);
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
        
        % process prcp
        if ismember('prcp',invar)
            pcp_d = ncread(infile_stndaily,'prcp',locstart,loclen);
            pcp_a = pcp_d./pcp_clim(:,i);
            ncwrite(outfile_stnanom,'prcp',pcp_a,locstart);
        end
        
        %tgq: for temperature, calculate the difference between station
        %temperature and ensemble mean temperature.
        if ismember('tmin',invar)
            tmin_d = ncread(infile_stndaily,'tmin',locstart,loclen);
            tmin_a = tmin_d - tmin_clim(:,i);
            ncwrite(outfile_stnanom,'tmin',tmin_a,locstart);
        end
        
        if ismember('tmax',invar)
            tmax_d = ncread(infile_stndaily,'tmax',locstart,loclen);
            tmax_a = tmax_d - tmax_clim(:,i);
            ncwrite(outfile_stnanom,'tmax',tmax_a,locstart);
        end
    end
end