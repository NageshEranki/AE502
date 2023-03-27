%{
        Determining a,i,e 
%}
clear;clc;

%   J2 zonal harmonic
% J2 = 0.00108;
J2 = 0.00196;

%   Planet radius (km)
% Re = 6370;
Re = 3390;

%   Gravitational param (km3/s2)
% mu = 3.986e5;
mu = 4.282e4;

%   Radius of perigee (km)
% rp=(600:12000)+Re;
rp=(400)+Re;

%   Time period (s)
% T=86164.20003652/3;
T = 24*3600+39*60+35;

%   Semi-major axis (km)
a = (T*sqrt(mu)/(2*pi))^(2/3);

%   Radius of apogee (km)
ra = 2*a-rp;

%   Eccentricity
e=(ra-rp)./(ra+rp);

%   Specific angular momentum
h = sqrt(rp.*mu.*(1+e));

%   Inclination (rad)
I=1.1071;
% I = 2.034;


% Omegadot = -1.5*(2*pi./T).*J2.*(Re./a).^2*cos(I)./(1-e.^2).^2;
Omegadot = -1.5*(2*pi./T).*J2.*(Re./a).^2*cos(I)./(1-e.^2).^2;

%   Nodal precession (deg/day)
% Omegadot = Omegadot*(180/pi)*86164.20003652;
Omegadot = Omegadot*(180/pi)*T;

% theta=0:0.01:2*pi;
% r=(h.^2./mu)*1./(1+e.*cos(theta));
% polarplot(theta,r)
% figure(1)
% subplot(2,1,1)
% plot(rp-Re,e)
% title("Eccentricity")
% subplot(2,1,2)
% plot(rp-Re,Omegadot)
% title("Nodal precession (deg/day)")
% xlabel('Perigee altitude (km)')


%   True anomalies that can possibly observe 
%   the former USSR (rad)
f1= 160*pi/180;
f2= 180*pi/180;

%   Convert to Eccentric anomalies (rad)
E1=2*atan(sqrt((1-e)./(1+e)).*tan(f1/2));
E2=2*pi-2*atan(sqrt((1-e)./(1+e)).*tan(f2/2));

%   Convert to Mean anomalies (rad)
M1 = E1 - e.*sin(E1);
M2 = E2 - e.*sin(E2);

%   Delta 't' between those mean anomalies
t1 = M1./(2*pi)*T;
t2 = M2./(2*pi)*T;

% figure(2)
% plot(rp-Re,2*(t2-t1)/60)
% title("Dwell time over apogee (min)")
% xlabel('Perigee altitude (km)')