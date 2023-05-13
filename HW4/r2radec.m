function radec=r2radec(r)
    %{

        Geocentric Right Ascension/Declination

        Input:
            
            r:  Distance vector (km) in non-rotating
                Geocentric coords

        Output:
            
            RA  : Right-ascension (deg)
            DEC : Declination (deg)
    %}

    rmag = vecnorm(r,2,2);

    %   Direction cosines
    l=r(:,1)./rmag; m=r(:,2)./rmag; n=r(:,3)./rmag;

    dec = asind(n);

    ra = acosd(l./cosd(dec)).*(m>0) + (360-acosd(l./cosd(dec))).*(m<=0);
  
    radec=[ra dec];
end