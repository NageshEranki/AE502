%{
        Determining a,i,e 
%}
clear;clc;

J2 = 0.00196;
R = 3390;
mu = 4.282e4;

rp=400+R;

% T=86164.20003652/3;
T=24*3600+39*60+35;

a = (T*sqrt(mu)/(2*pi))^(2/3);

ra = 2*a-rp;

e=(ra-rp)/(ra+rp);

I=1.1071;


Omegadot = -1.5*(2*pi/T)*J2*(R/a)^2*cos(I)/(1-e^2)^2

% Omegadot = Omegadot*(180/pi)*86164.20003652