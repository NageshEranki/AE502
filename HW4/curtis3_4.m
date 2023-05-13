%{
    My Two-Body propagator, using
    Curtis Algo 3-4

        Input args:
            1. r0 : Initial pos vector
            2. v0 : Initial velocity vector
            3. Dt : Time of flight
        Output:
            1. r  : Final pos vector (Km)
            2. v  : Final velocity vector (Km/s)
%}

function [r,v] = curtis3_4(r0,v0,Dt)

    %{ 
        Development notes
    
        1. All units in SI (km,km/s)
        2. Hard code initial physical quantities

    %}

    %   Heliocentric gravitational param.
    % mu = 1.32712440042e20;
    
    %   Geocentric grav. param.
    mu = 398600;
    
    %   Magnitudes of r0/v0/r/v
    r0mag = norm(r0);
    v0mag = norm(v0);

    %   Radial comp. of v0
    vr0mag = dot(r0,v0)/r0mag;

    %   Reciprocal of semi-major axis
    alpha = 2/r0mag - v0mag^2/mu;

    %   Calculate universal anomaly chi
    %   using m-function curtis3_3()
    chi = curtis3_3(Dt,r0mag,vr0mag,alpha);

    %   Compute Lagrange coefficients f/g
    f = 1 - chi^2*C(alpha*chi^2)/r0mag;
    g = Dt - 1/sqrt(mu)*chi^3*S(alpha*chi^2);

    %   Compute final pos vector
    r = f*r0 + g*v0;

    %   Final distance
    rmag = norm(r);

    %   Compute Lagrange coefficients fdot/gdot
    fdot = sqrt(mu)*(alpha*chi^3*S(alpha*chi^2)-chi)/(r0mag*rmag);
    gdot = 1-chi^2*C(alpha*chi^2)/rmag;

    %   Compute final velocity
    v = fdot*r0 + gdot*v0;

end