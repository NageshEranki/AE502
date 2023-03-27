function coe = coe_from_sv(r,v,mu)
    %{
        Input: 
                State vector in (km)/(km/s)
                Gravitational param in (km3/s2)
        Output:
                Orbital elements:
                    h:      Angular momentum
                    i:      Inclination angle (rad)
                    Omega:  RAAN
                    e:      Eccentricity
                    w:      Arugment of periapsis
                    f:      True anomaly
    %}

    %   Magnitudes of dist/velocities
    rmag = vecnorm(r,2,2);
    vmag = vecnorm(v,2,2);

    %   radial v
    vr = dot(r,v,2)./rmag;

    %   Angular momentum vector
    h = cross(r,v,2);

    %   mag of h
    hmag = vecnorm(h,2,2);

    %   inclination
    i = acos(h(:,3)./hmag);

    %   Node line calc
    K = zeros(size(r,1),3);
    K(:,3) = 1;

    N = cross(K,h,2);

    %   Mag of N
    Nmag = vecnorm(N,2,2);

    %   RAAN
    Omega = (acos(N(:,1)./Nmag)).*(N(:,2)>=0) + (2*pi-acos(N(:,1)./Nmag)).*(N(:,2)<0);
    

    %   Eccentricity vector
    e = 1/mu*( (vmag.^2-mu./rmag).*r-(rmag.*vr).*v );

    %   Eccentricity mag
    emag = vecnorm(e,2,2);

    %   Arg of perigee
    w = (acos((dot(N,e,2))./(Nmag.*emag))).*(e(:,3)>=0) + (2*pi-acos((dot(N,e,2))./(Nmag.*emag))).*(e(:,3)<0);

    %   True anomaly
    f = acos((dot(e,r,2))./(emag.*rmag)).*(vr>=0) + (2*pi-acos((dot(e,r,2))./(emag.*rmag))).*(vr<0);

    coe(:,1) = hmag;
    coe(:,2) = i;
    coe(:,3) = Omega;
    coe(:,4) = emag;
    coe(:,5) = w;
    coe(:,6) = f;

end