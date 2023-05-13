%{
    My Two-Body propagator, using
    Curtis 3.3.

    Input args: 
        1. Dt :     Time of flight (s)
        2. r0 :     Inital dist to obj (Km)
        3. Vr0:     Initial radial speed of obj (Km/s)
        4. alpha:   Reciprocal of semi-major axis

    Output:
        1. chi
%}

function chi = curtis3_3( Dt, r0, Vr0, alpha )

    %{ 
        Development notes
    
        1. Units: Km, s 
        2. Hard code initial physical quantities

    %}

    %   Iteration tolerance
    tol = 1e-4;

    %   Heliocentric gravitational param.
    % mu = 1.32712440042e20;

    %   Geocentric grav. param.
    mu = 398600;

    %   Init estimate for chi, ratio
    chi = sqrt(mu)*abs(alpha)*Dt;
    %ratio = 1;

    %   Iteration begin:
    while 1
        %   Arg for Stumpff functions
        z = alpha*chi^2;
        %   Numerator in newton step
        num = (r0*Vr0/sqrt(mu))*chi^2*C(z) + (1-alpha*r0)*chi^3*S(z) + r0*chi - sqrt(mu)*Dt;
        %   Denominator in newton step
        den = (r0*Vr0/sqrt(mu))*chi*(1-alpha*chi^2*S(z)) + (1-alpha*r0)*chi^2*C(z) + r0;
        %   Update ratio, chi
        ratio = num/den;
        %   If ratio within tol, break
        if abs(ratio) > tol
            chi = chi - ratio;
        else
            break
        end
    end

end

%{
function res = C(z)
    %{
        Stumpff function 'C(z)'
    %}
    
    if (z>0)
        
        res = (1-cosh(sqrt(z)))/z;

    elseif (z<0)
        
        res = (cosh(sqrt(-z))-1)/(-z);

    else
    
        res = 0.5;

    end
end
%}

%{
function res = S(z)
    %{
        Stumpff function 'S(z)'
    %}

    if (z>0)
        res = (sqrt(z)-sin(sqrt(z)))/sqrt(z)^3;
    elseif (z<0)
        res = (sinh(sqrt(-z))-sqrt(-z))/sqrt(-z)^3;
    else
        res = 1.0/6;
    end
end
%}