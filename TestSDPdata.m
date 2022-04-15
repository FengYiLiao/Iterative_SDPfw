%TestSDPdata
clc;clear; 
addpath('Module\');
dx = 2;
Maxiter = 10;
load('SedumiData\n10m10.mat');
At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);
Objs =  InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,dx,Maxiter);