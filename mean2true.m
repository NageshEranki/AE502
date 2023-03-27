%{
    Mean to true anomaly

    In args:    e : Eccentricity
                M : Mean anomaly , rad
    out    :    theta: true anomaly, rad
%}

function theta = mean2true(M,e)

    %   First calculate the Eccentric anomaly

    if M<pi
        E = M + e/2;
    else
        E = M - e/2;
    end
    
    ratio = (E-e*sin(E)-M)/(1-e*cos(E));
    if abs(ratio)>1e-10
        E = E - ratio;
    end

    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end