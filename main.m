clear all
close all

load('Cortex.mat');
load('data.mat');
load('Gain.mat');
%% 
x0=zeros(15002,1);
z1=double(Cortex.VertConn>0);
G0=speye(15002)+(z1-diag(sum(z1,2)))/16;
Gk=G0*G0*G0;
G=sparse(double(Gk>0));


nbV = diag(sum(Cortex.VertConn,2))-0.9*Cortex.VertConn;
invnbv=inv(nbV);
[s_est]=  SI_DCB(B,Gain*invnbv,G);
s_est=invnbv*s_est;
