%%TestBQO
clc;clear; 
addpath('Module\');
dx = 10;
Maxiter = 10;
load('SedumiData\BQO\100_1.mat');
Objs =  InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,dx,Maxiter);