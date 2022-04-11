%TestSDPdata

clc;clear; 
dx = 1;
load('SedumiData\n10m10.mat');
At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);
Objs =  InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,dx);