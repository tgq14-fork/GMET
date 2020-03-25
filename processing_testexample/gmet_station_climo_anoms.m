clear all;
%tgq: calculate the ratio (prcp) or difference (temperature) between
%station data and ensemble mean, and save those file in 

c_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/stndata/climo'; % store climatological mean of prcp/temp
a_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/stndata/anom'; % store the anomaly (ratio or difference) or prcp/temp
d_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/stndata'; % original station time series

g_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/inputs'; % grid information

%% load data

fname = sprintf('%s/list.txt',c_path);
fid = fopen(fname);
list = textscan(fid,'%s');
fclose(fid);
nsta = length(list{1});

%grid
fname = sprintf('%s/gridinfo.0625.w_subset.nc',g_path);
lat = ncread(fname,'latitude');
lon = ncread(fname,'longitude');
lat1d = lat(:,1);
lon1d = lon(1,:);

%tgq: climate mean of prcp tmean trange using all ensemble members
fname = sprintf('%s/../outputs/ens_forc_ecai_climo.mean.nc',g_path);
pcp_clim = ncread(fname,'pcp');
tmean_clim = ncread(fname,'t_mean');
trange_clim = ncread(fname,'t_range');
tmax_clim = tmean_clim + trange_clim/2;
tmin_clim = tmean_clim - trange_clim/2;


for i = 1:nsta
    %tgq: for each station, find the nearest grid cells
    %extract the mean of all ensembles, and calculate the ratio (pcp_a)
    %between station prcp and ensemble mean prcp
    try
        fname = sprintf('%s/%s',c_path,char(list{1}(i)));
        pcp = ncread(fname,'prcp');
        fname = sprintf('%s/%s',d_path,char(list{1}(i)));
        pcp_d = ncread(fname,'prcp');
        slat = ncread(fname,'latitude');
        slon = ncread(fname,'longitude');
        
        [~, y] = min(abs(lat1d-slat));
        [~, x] = min(abs(lon1d-slon));
        
        pcp_a = pcp_d./pcp_clim(x,y);
        
        fname = sprintf('%s/%s',a_path,char(list{1}(i)));
        ncwrite(fname,'prcp',pcp_a);
    catch
        fprintf(1,'No prcp: %s\n',char(list{1}(i)));
    end
    
    %tgq: for temperature, calculate the difference between station
    %temperature and ensemble mean temperature.
    try
        fname = sprintf('%s/%s',c_path,char(list{1}(i)));
        tmax = ncread(fname,'tmax');
        tmin = ncread(fname,'tmin');
   
        fname = sprintf('%s/%s',d_path,char(list{1}(i)));
        tmax_d = ncread(fname,'tmax');
        tmin_d = ncread(fname,'tmin');
    
        [~, y] = min(abs(lat1d-slat));
        [~, x] = min(abs(lon1d-slon));

        tmax_a = tmax_d - tmax_clim(x,y);
        tmin_a = tmin_d - tmin_clim(X,y);
   
        fname = sprintf('%s/%s',a_path,char(list{1}(i)));
        ncwrite(fname,'tmax',tmax_a);
        ncwrite(fname,'tmin',tmin_a);
    catch
    end
end
