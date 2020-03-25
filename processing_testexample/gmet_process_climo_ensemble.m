clear all;
%tgq: run before gmet_station_climo_anoms.m
%read all ensemble members, and sort prcp from high to low, and calculate
%the std of tmean and trange
%save in netcdf

%d_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/outputs/climo';
d_path = '/glade/p/ral/hap/anewman/GMET_test_ecai/outputs';

nens = 100;

for i = 1:nens
    fname = sprintf('%s/ens_forc_ecai_climo.%03d.nc',d_path,i);
    ens.pcp(i,:,:) = squeeze(ncread(fname,'pcp'));
    ens.tmean(i,:,:) = squeeze(ncread(fname,'t_mean'));
    ens.trange(i,:,:) = squeeze(ncread(fname,'t_range'));
end

%% compute stats
ens.pcp_sort = sort(ens.pcp,1);
ens.tmean_std = std(ens.tmean,0,1);
ens.trange_std = std(ens.trange,0,1);

%% write to output
%pcp
fname = sprintf('%s/ENS_CLIMO_SORTED_01.nc',d_path);

%ncwrite(fname,'pcp_sorted',permute(ens.pcp_sort,[2 3 1]));
ncwrite(fname,'pcp_sorted',ens.pcp_sort);

%temp
fname = sprintf('%s/ENS_CLIMO_01.nc',d_path);

ncwrite(fname,'t_mean',squeeze(mean(ens.tmean,1)));
ncwrite(fname,'t_range',squeeze(mean(ens.trange,1)));

fname = sprintf('%s/ENS_UNCERT_01.nc',d_path);
ncwrite(fname,'t_mean_stddev',squeeze(ens.tmean_std));
ncwrite(fname,'t_range_stddev',squeeze(ens.trange_std));
