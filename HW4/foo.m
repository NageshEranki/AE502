clear;clc;

global mu data

mu = 398600;

%   Topocentric data from Champaign
data=readmatrix('AE502_HW4.csv');

%   Lat/Lon/elev of Champaign, IL :)
%   Units in deg, km
lat = 40.1163;  lon = -88.2434; elev = .233;

%   Compute local sidereal times first (theta, in deg)
times=datetime(data(:,end),'ConvertFrom','modifiedjuliandate');
theta=siderealTime(juliandate(times)) + lon;

%   Locate the observer on the Earth's surface.
%   Model used : WGS84 ellipsoid

%   Model params
Re = 6378.1370;     % Equatorial radius (km)
f = 0.00335;        % Oblateness

%   Observer positions(km) in geocentric equatorial coordinate frame
R = zeros(size(data,1),3);

R(:,1) = (Re/sqrt(1-(2*f-f^2)*sind(lat)^2)+elev).*cosd(lat).*cosd(theta);
R(:,2) = (Re/sqrt(1-(2*f-f^2)*sind(lat)^2)+elev).*cosd(lat).*sind(theta);
R(:,3) = ((Re*(1-f)^2)/sqrt(1-(2*f-f^2)*sind(lat)^2)+elev).*sind(lat);


%   Relative position vector directions
rho = zeros(size(data,1),3);

rho(:,1) =  cosd(data(:,2)).*cosd(data(:,1));
rho(:,2) =  cosd(data(:,2)).*sind(data(:,1));
rho(:,3) =  sind(data(:,2));


%   Perform initial orbit determination
%   using Gauss' method with three pairs
%   of RA/DEC measurements


rho1=rho(1,:);rho2=rho(2,:);rho3=rho(3,:);
R1=R(1,:);R2=R(2,:);R3=R(3,:);
t=data(:,5)*86400;

t1=t(1);t2=t(2);t3=t(3);

[r20,v20]=gauss(rho1,rho2,rho3, ...
    R1,R2,R3, ...
    t1,t2,t3);

%   Least-Squares(linearized) iteration to determine
%   best fit state vector (r20,v20)

z=data(2:end,1:2); z=z(:);  %   Obs data

tol = 1e-4;

% while abs(dcost)>tol
for j=1:100
   
    %-----------------------------------------
    %   Model prediction func 'h'
    %-----------------------------------------
    h0=h(r20,v20,t,R);
    % rs=zeros(size(data,1)-3,3);
    % 
    % for i=4:size(data,1)
    %     r=curtis3_4(r20,v20,t(i)-t(2));
    %     rs(i-3,:)=r;
    % end
    % h0 = r2radec(rs); h0=h0(:);
    %-----------------------------------------

    
    
    %-----------------------------------------
    %   Compute the normal equations matrix 'H'
    %-----------------------------------------
    dr=0.0005; dv=0.0005;
    H = zeros(16,6);
    hx1=h(r20+dr*[1,0,0],v20,t,R)-h(r20+dr*[-1,0,0],v20,t,R);
    hx2=h(r20+dr*[0,1,0],v20,t,R)-h(r20+dr*[0,-1,0],v20,t,R);
    hx3=h(r20+dr*[0,0,1],v20,t,R)-h(r20+dr*[0,0,-1],v20,t,R);
    hx4=h(r20,v20+dv*[1,0,0],t,R)-h(r20,v20+dv*[-1,0,0],t,R);
    hx5=h(r20,v20+dv*[0,1,0],t,R)-h(r20,v20+dv*[0,-1,0],t,R);
    hx6=h(r20,v20+dv*[0,0,1],t,R)-h(r20,v20+dv*[0,0,-1],t,R);
    H(:,1)=(hx1)/(2*dr);
    H(:,2)=(hx2)/(2*dr);
    H(:,3)=(hx3)/(2*dr);
    H(:,4)=(hx4)/(2*dv);
    H(:,5)=(hx5)/(2*dv);
    H(:,6)=(hx6)/(2*dv);
    %-----------------------------------------

    dx = linsolve(H'*H,H'*(z-h0));

    r20=r20+dx(1:3)';
    v20=v20+dx(4:6)';

    dot(z-h0,z-h0)
    
end



function out=res(r20,v20)

    global data

    times=datetime(data(:,end),'ConvertFrom','modifiedjuliandate');
    t=juliandate(times)*86400;
    %   Predict future values by propagating
    %   r20,v20
    rs=zeros(size(data,1)-3,3);
    
    for i=4:size(data,1)
        r=curtis3_4(r20,v20,t(i)-t(2));
        rs(i-3,:)=r;
    end
    
    %   Model predictions
    h = r2radec(rs);

    %   Compute initial residuals
    z=data(4:end,1:2);
    z=z(:);h=h(:);

    out=dot(z-h,z-h);
end

% R1  =   [3489.8 , 3430.2 , 4078.5];
% R2  =   [3460.1 , 3460.1 , 4078.5];
% R3  =   [3429.9 , 3490.1 , 4078.5];
% 
% rho1 =  [0.71643,0.68074 , -0.1527];
% rho2 =  [0.56897,0.79531,-0.20917];
% rho3 =  [0.41841,0.87,-0.26059];
% 
% t1  =   0;
% t2  =   118.1;
% t3  =   237.58;
% 
% [r,v,r_old,v_old]=gauss(rho1,rho2,rho3,R1,R2,R3,t1,t2,t3);


% clear all; clc
% global mu
% deg = pi/180;
% mu = 398600;
% Re = 6378;
% f  = 1/298.26;
% % %...Data declaration for Example 5.11:
% H = 1;
% phi = 40*deg;
% t = [ 0 118.104 237.577];
% 
% ra = [ 43.5365 54.4196 64.3178]*deg;
% dec = [-8.78334 -12.0739 -15.1054]*deg;
% theta = [ 44.5065 45.000 45.4992]*deg;
% % %...
% % %...Equations 5.64, 5.76 and 5.79:
% fac1 = Re/sqrt(1-(2*f - f*f)*sin(phi)^2);
% fac2 = (Re*(1-f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*sin(phi);
% for i = 1:3
% R(i,1) = (fac1 + H)*cos(phi)*cos(theta(i));
% R(i,2) = (fac1 + H)*cos(phi)*sin(theta(i));
% R(i,3) = fac2;
% rho(i,1) = cos(dec(i))*cos(ra(i));
% rho(i,2) = cos(dec(i))*sin(ra(i));
% rho(i,3) = sin(dec(i));
% end
% %...Algorithms 5.5 and 5.6:
% [r, v, r_old, v_old] = gauss(rho(1,:), rho(2,:), rho(3,:), ...
% R(1,:), R(2,:), R(3,:), ...
% t(1), t(2), t(3));
