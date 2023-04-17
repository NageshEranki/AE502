clear;clc;

%{
        Propagating the orbit of the satellite
        around the Earth in the presence of the
        Sun.
%}

%   Perturbation param
w = 0.01;

%   Initial states (Keplerian, SIunits, radian only!)
a0 = 1; e0 = 0.5; i0 = deg2rad(45);
O0 = 0; wp0= 0;   M0 = 0;
%   Compute (initial?) mean motion
n0 = sqrt(1/a0^3);

%   Compute initial state (Delaunay elements)
%   Actions
L0 = n0*a0^2;
G0 = L0*sqrt(1-e0^2);
H0 = G0*cos(i0);
%   Angles
l0 = M0 + wp0 + O0;
g0 = wp0 + O0;
h0 = O0;

%   Define time series
t = 0:0.1:100;

%   Propagate the angle vars
%   (Easy-peasy lemon squeezy :D )
l = 1./(L0^3).*t+l0;
g = g0*ones(size(l,1),size(l,2));
h = w.*t + h0;

L = L0*ones(size(l,1),size(l,2));     %   Actions do not change
G = G0*ones(size(l,1),size(l,2));
H = H0*ones(size(l,1),size(l,2));

%   Recover Keplerian elements from
%   the propagated Delaunay elements
O = h;
wp = g-h;
M = l-g;
a = L.^2;
e = sqrt(1-G.^2./L.^2);
i = acos(H./G);

%   In addition, require angular mom and true anomaly

hang = sqrt(a.*(1-e.^2));
f = mean2true(M,e(1));

%   Calc state vector from keplerian
%   elements via '402 equations' or
%   Curtis' algorithms

r=zeros(length(M),3);
%v=zeros(length(M),3);

for j = 1:length(M)
    [rj,~] = sv_from_coe([hang(j),e(j),O(j),i(j),wp(j),f(j)],1);
    r(j,:)= rj;   %v(j,:) = vj;
end

figure(2)
plot3(r(:,1),r(:,2),r(:,3))
hold on
plot3(0,0,0,'.k','MarkerSize',10)
xlabel('x')
ylabel('y')
zlabel('z')
