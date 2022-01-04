clc;clear all;close all;
fid=fopen('mydata.grd','r');
cdum=fread(fid,4,'uint8=>char');
fscanf(fid,'\n');%读取下一行
nnx_nny=fscanf(fid,'%d',2);
fscanf(fid,'\n');
lon=fscanf(fid,'%f',2);
fscanf(fid,'\n');
lat=fscanf(fid,'%f',2);
fscanf(fid,'\n');
min_max=fscanf(fid,'%f',2);
fscanf(fid,'\n');
value_r=fscanf(fid,'%f',[im,jm]);%im为grd文件的列数，jm为grd文件的行数
fclose(fid);