% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear

Path_script='/scratch/gwf/gwf_cmt/gut428/SCRF/rndnum_scripts';
exefile='/scratch/gwf/gwf_cmt/gut428/SCRF/generate_rndnum_exp2p.exe';
varall={'prcp','tmean','trange','tdew'};
for v=1:4
    var=varall{v};
    
    start_date=19500101;
    stop_date=19500102;
    ensnum=[1,1];
    
    if ismember(var,{'prcp','trange','wind'})
        cross_cc_flag=-1; % >0: use cross option
        weight_judge=1.5;
    else
        cross_cc_flag=-1;
        weight_judge=0.7;
    end
    
    exp2p_file=sprintf('/home/gut428/EMGLB/8DeterministicEstimate/Clen/exp2p_global_parameter_%s.nc',var);
    if strcmp(var,'prcp') || strcmp(var,'wind')
        if cross_cc_flag>0
            cc_file='/home/gut428/EMGLB/8DeterministicEstimate/auto_cross_corr/cross_cc_prcp_vs_trange.nc';
        else
            cc_file=sprintf('/home/gut428/EMGLB/8DeterministicEstimate/auto_cross_corr/auto_cc_Spearman_%s.nc',var);
        end
        cross_file_prefix='/scratch/gwf/gwf_cmt/gut428/SCRF/rndnum_trange/scrf_';
    else
        cc_file=sprintf('/home/gut428/EMGLB/8DeterministicEstimate/auto_cross_corr/auto_cc_Pearson_%s.nc',var);
        cross_file_prefix='xxx';
    end
    
    grid_name='/scratch/gwf/gwf_cmt/gut428/SCRF/topography_attribute_global.nc';
    path_spcorr=sprintf('/scratch/gwf/gwf_cmt/gut428/SCRF/spcc_struct_%s_wj_%s',var,num2str(weight_judge));
    out_spcorr_prefix=[path_spcorr, '/spcorr_'];
    path_rndnum=sprintf('/scratch/gwf/gwf_cmt/gut428/SCRF/rndnum_%s_wj_%s',var,num2str(weight_judge));
    out_rndnum_prefix=[path_rndnum,'/scrf_'];
    
    % make dir if necessary
    if ~exist(Path_script,'dir'); mkdir(Path_script); end
    if ~exist(path_spcorr,'dir'); mkdir(path_spcorr); end
    if ~exist(path_rndnum,'dir'); mkdir(path_rndnum); end
    
    for i=ensnum(1):ensnum(2)
        start_ens=i;
        stop_ens=i;
        
        % generate namelist file
        % file1=['namelist_',var,'_ens_',num2str(i),'.txt'];
        file1=['namelist_',var,'_ens_',num2str(i),'_wj_',num2str(weight_judge),'.txt'];
        outfile1=[Path_script,'/',file1];
        fid=fopen(outfile1,'w');
        fprintf(fid,'&PARAMS\n');
        fprintf(fid,['start_date = ',num2str(start_date),'\n']);
        fprintf(fid,['stop_date	= ',num2str(stop_date),'\n']);
        fprintf(fid,['start_ens	= ',num2str(start_ens),'\n']);
        fprintf(fid,['stop_ens	= ',num2str(stop_ens),'\n']);
        fprintf(fid,['cross_cc_flag	= ',num2str(cross_cc_flag),'\n']);
        fprintf(fid,['weight_judge	= ',num2str(weight_judge),'\n']);
        fprintf(fid,['exp2p_file	= "',exp2p_file,'"\n']);
        fprintf(fid,['cross_file_prefix	= "',cross_file_prefix,'"\n']);
        fprintf(fid,['cc_file	= "',cc_file,'"\n']);
        fprintf(fid,['grid_name	= "',grid_name,'"\n']);
        fprintf(fid,['out_spcorr_prefix	= "',out_spcorr_prefix,'"\n']);
        fprintf(fid,['out_rndnum_prefix	= "',out_rndnum_prefix,'"\n']);
        fprintf(fid,'/\n');
        fclose(fid);
        
        % generate slurm scripts
        % outfile2=[Path_script,'/run_',var,'_ens_',num2str(i),'.sh'];
        outfile2=[Path_script,'/run_',var,'_ens_',num2str(i),'_wj_',num2str(weight_judge),'.sh'];
        fidout=fopen(outfile2,'w');
        fprintf(fidout,'#!/bin/bash\n');
        fprintf(fidout,['#SBATCH --job-name=rndnum_',var,'_',num2str(i),'\n']);
        fprintf(fidout,'#SBATCH --account=hpc_c_giws_clark\n');
        fprintf(fidout,['#SBATCH --time=0-3:00:00\n']);
        fprintf(fidout,'#SBATCH --mem=20G\n');
        fprintf(fidout,['chmod a+x ',exefile,'\n']);
        fprintf(fidout,[exefile,' ',file1,'\n']);
        % fprintf(fidout,'rm *.out\n');
        fclose(fidout);
    end
end