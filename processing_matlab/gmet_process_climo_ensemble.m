clear
%tgq: run before gmet_station_climo_anoms.m
%read all ensemble members, and sort prcp from high to low, and calculate
%the std of tmean and trange
%save in netcdf
%Don't change output file names because they are fixed in eCAI codes!

EnsPath='/home/gut428/GMET/eCAI_EMDNA/output/NA_reg_ens';
year=2018:2018;
month=1:12;
nens = 100;

for yy=1:length(year)
    EnsPathyy = [EnsPath,'/',num2str(year(yy))];
    for mm=1:length(month)
        
        outfile_stddev = sprintf('%s/ENS_UNCERT_%02d.nc',EnsPathyy,mm);
        outfile_psort = sprintf('%s/ENS_CLIMO_SORTED_%02d.nc',EnsPathyy,mm);
        outfile_mean = sprintf('%s/ENS_CLIMO_%02d.nc',EnsPathyy,mm);
        
        if ~exist(outfile_stddev,'file') || ~exist(outfile_psort,'file') || ~exist(outfile_mean,'file')
            for i = 1:nens
                fprintf('Reading member %d--%d for %04d-%03d\n',i,nens,year(yy),mm);
                ensfile = sprintf('%s/ensemble_climo_%04d%02d.%03d.nc',EnsPathyy,year(yy),mm,i);
                ensdata.pcp(i,:,:) = squeeze(ncread(ensfile,'pcp'));
                ensdata.tmean(i,:,:) = squeeze(ncread(ensfile,'t_mean'));
                ensdata.trange(i,:,:) = squeeze(ncread(ensfile,'t_range'));
            end
            
            % compute stats
            ensdata.pcp_sort = sort(ensdata.pcp,1);
            ensdata.tmean_std = std(ensdata.tmean,0,1);
            ensdata.trange_std = std(ensdata.trange,0,1);
            
            % write to output
            [ens,x,y]=size(ensdata.pcp);
            %pcp sorted
            try
                nccreate(outfile_psort,'pcp_sorted','Datatype','single','Dimensions',{'ens',ens,'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('pcp_sorted exists in %s\n',outfile_psort);
            end
            ncwrite(outfile_psort,'pcp_sorted',ensdata.pcp_sort);
            
            %pcp and temp mean
            try
                nccreate(outfile_mean,'pcp','Datatype','single','Dimensions',{'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('pcp exists in %s\n',outfile_mean);
            end
            ncwrite(outfile_mean,'pcp',squeeze(mean(ensdata.pcp,1)));
            try
                nccreate(outfile_mean,'t_mean','Datatype','single','Dimensions',{'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('t_mean exists in %s\n',outfile_mean);
            end
            ncwrite(outfile_mean,'t_mean',squeeze(mean(ensdata.tmean,1)));
            try
                nccreate(outfile_mean,'t_range','Datatype','single','Dimensions',{'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('t_range exists in %s\n',outfile_mean);
            end
            ncwrite(outfile_mean,'t_range',squeeze(mean(ensdata.trange,1)));
            
            % temp stddev
            try
                nccreate(outfile_stddev,'t_mean_stddev','Datatype','single','Dimensions',{'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('t_mean_stddev exists in %s\n',outfile_stddev);
            end
            ncwrite(outfile_stddev,'t_mean_stddev',squeeze(ensdata.tmean_std));
            try
                nccreate(outfile_stddev,'t_range_stddev','Datatype','single','Dimensions',{'x',x,'y',y},'DeflateLevel',9);
            catch
                fprintf('t_range_stddev exists in %s\n',outfile_stddev);
            end
            ncwrite(outfile_stddev,'t_range_stddev',squeeze(ensdata.trange_std));
        end
    end
end