% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear

% gridinfo='/home/gut428/scratch/GMET/StnGridInfo/gridinfo_whole.nc';
% outpath='/home/gut428/scratch/GMET/';
gridinfo='/home/gut428/scratch/GMET/StnGridInfo/gridinfo_whole.nc';
outpath='/home/gut428/scratch/GMET/EMDNA_out/SCRF';

Path_script='/Users/localuser/Downloads/scripts';
sens=1;
eens=3;

for i=sens:eens
    file1=['rndnum_',num2str(i),'.txt'];
    outfile1=[Path_script,'/',file1];
    fid=fopen(outfile1,'w');
    fprintf(fid,'&PARAMS\n');
    fprintf(fid,['start_time	= 999\n']);
    fprintf(fid,['ntimes		= 999\n']);
    fprintf(fid,['start_ens	= ',num2str(i),'\n']);
    fprintf(fid,['stop_ens	= ',num2str(i),'\n']);
    fprintf(fid,['clen		= 999\n']);
    fprintf(fid,['grid_name	= "',gridinfo,'"\n']);
    fprintf(fid,['out_forc_name_base	= "',outpath,'"\n']);
    fprintf(fid,['in_regr_name	= "xxxxx.nc"\n']);
    fprintf(fid,['time_mode       = "999"\n']); % different time mode
    fprintf(fid,'/\n');
    fclose(fid);
    
    

    outfile2=[Path_script,'/run_rndnum_',num2str(i),'.txt'];
    fidout=fopen(outfile2,'w');
    fprintf(fidout,'#!/bin/bash\n');
    fprintf(fidout,['#SBATCH --job-name=rndnum_',num2str(i),'\n']);
    fprintf(fidout,'#SBATCH --account=rpp-kshook\n');
    fprintf(fidout,['#SBATCH --time=0-30:00:00\n']);
    fprintf(fidout,'#SBATCH --mem=20G\n');
    
    exefile=['/home/gut428/scratch/GMET/EMDNA_out/generate_rndnum.exe'];
    fprintf(fidout,['chmod a+x ',exefile,'\n']);
    fprintf(fidout,[exefile,' ',file1,'\n']);
    %             fprintf(fidout,'rm *.out\n');
    fclose(fidout); 
end