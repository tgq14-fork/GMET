nx=1800; ny=3600;
ngrids=nx*ny;
slim=49; vv=5;
file='/Users/localuser/Research/EMGLB/RandomNumber/spcc_struct_prcp/spcorr_month_ 1-old';

fid=fopen(file,'rb');
len=fread(fid,1,'uint32');
sp_wght_var=fread(fid,ngrids*slim,'float64');
sp_wght_var=reshape(sp_wght_var,nx, ny, slim);
fclose(fid);


% sp_ipos_var=fread(fid,ngrids*slim,'int32');
% sp_ipos_var=reshape(sp_ipos_var,nx, ny,slim);
% 
% sp_jpos_var=fread(fid,ngrids*slim,'int32');
% sp_jpos_var=reshape(sp_jpos_var,nx, ny,slim);
% 
% sp_num_var=fread(fid,ngrids,'int32');
% sp_num_var=reshape(sp_num_var,nx, ny);
% 
% sp_sdev_var=fread(fid,ngrids,'float64');
% sp_sdev_var=reshape(sp_sdev_var,nx, ny);
% 
% iorder1d=fread(fid,ngrids,'int32');
% jorder1d=fread(fid,ngrids,'int32');



