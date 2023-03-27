% HW2: prob 3
clear;clc

%  Parameters
mu = 3.986e5;           % Km3/s2
R = 6370;               % Km
J2 = 0.00108;           % 

% Define orbital elements
a = 26600;              % km
i = 1.10654;            % rad
e = 0.74;               % ()
w = 5 * pi/180;         % rad
omega = 90 * pi/180;    % rad
M0 = 10 * pi/180;       % rad

%  Compute necessary elements

h = sqrt(mu*a*(1-e^2));
f0 = mean2true(M0,e);

%  Compute initial state vector

[r0,v0] = sv_from_coe([h,e,omega,i,w,f0],mu);

%  Begin Cowell's method

%   Integration time
T = 100;    %   Days
tspan = 0:200:(T*24*3600);
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y]=ode45( @(t,y) odeprob2(y,mu,J2,R) , tspan , [r0';v0'] , options );

%   Convert results in state vectors to orbital elements
coe=coe_from_sv(y(:,1:3),y(:,4:6),mu);

%   Calculate the rest of the elements
a=(coe(:,1).^2/mu).*1./(1-coe(:,4).^2);
E=2*atan(sqrt((1-coe(:,4))./(1+coe(:,4))).*tan(coe(:,6)./2));
M=E-e.*sin(E);

figure(1)
plot(t/(24*3600),a)
title('Semi-major axis (km)')
xlabel('Time (days)')

figure(2)
plot(t/(24*3600),coe(:,2)*(180/pi));
title('Inclination angle (deg)')
xlabel('Time (days)')

figure(3)
plot(t/(24*3600),coe(:,4));
title('Eccentricity')
xlabel('Time (days)')

figure(4)
plot(t/(24*3600),coe(:,5)*(180/pi));
title('Argument of perigee (deg)')
xlabel('Time (days)')

figure(5)
plot(t/(24*3600),coe(:,3)*(180/pi));
title('RAAN (deg)')
xlabel('Time (days)')
