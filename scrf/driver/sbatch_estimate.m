% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear

% gridinfo='/home/gut428/scratch/GMET/StnGridInfo/gridinfo_whole.nc';
% outpath='/home/gut428/scratch/GMET/';
gridinfo='/home/gut428/scratch/GMET/StnGridInfo/gridinfo_whole.nc';
scrfpath='/home/gut428/scratch/GMET/EMDNA_out/SCRF';
inpath='/home/gut428/scratch/GMET/GMET_OIinput';
outpath='/home/gut428/scratch/GMET/EMDNA_out/Estimate';
outsuffix='-scale1-pop3';

Path_script='/Users/localuser/Downloads/scripts';
sens=1;
eens=100;



for i=2016:2016
    outpathi=[outpath,'/',num2str(i),outsuffix];
    for j=1:12
        file1=['est_',num2str(i*100+j),'.txt'];
        outfile1=[Path_script,'/',file1];
        infile=[inpath,'/reg_',num2str(i*100+j),'.nc'];
        outbase=[outpathi,'/EMDNA_',num2str(i*100+j)];
        fid=fopen(outfile1,'w');
        fprintf(fid,'&PARAMS\n');
        fprintf(fid,['start_time	= 1\n']);
        fprintf(fid,['ntimes		= ',num2str(eomday(i,j)),'\n']);
        fprintf(fid,['start_ens	= ',num2str(sens),'\n']);
        fprintf(fid,['stop_ens	= ',num2str(eens),'\n']);
        fprintf(fid,['clen		= 999\n']);
        fprintf(fid,['grid_name	= "',gridinfo,'"\n']);
        fprintf(fid,['out_forc_name_base	= "',outbase,'"\n']);
        fprintf(fid,['in_regr_name	= "',infile,'"\n']);
        fprintf(fid,['time_mode       = "',scrfpath,'"\n']); % different time mode
        fprintf(fid,'/\n');
        fclose(fid);
        
        outfile2=[Path_script,'/run_est_',num2str(i*100+j),'.sh'];
        fidout=fopen(outfile2,'w');
        fprintf(fidout,'#!/bin/bash\n');
        fprintf(fidout,['#SBATCH --job-name=est_',num2str(i*100+j),'\n']);
        fprintf(fidout,'#SBATCH --account=rpp-kshook\n');
        fprintf(fidout,['#SBATCH --time=0-2:00:00\n']);
        fprintf(fidout,'#SBATCH --mem=10G\n');
        
        fprintf(fidout,['mkdir -p ', outpathi,'\n']);
        exefile=['/home/gut428/scratch/GMET/EMDNA_out/generate_estimate',outsuffix,'.exe'];
        fprintf(fidout,['chmod a+x ',exefile,'\n']);
        fprintf(fidout,[exefile,' ',file1,'\n']);
        %             fprintf(fidout,'rm *.out\n');
        fclose(fidout);
    end
end