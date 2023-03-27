function dydt=odeprob2(y,mu,J2,R)

    dydt = zeros(6,1);
    
    dydt(1) = y(4);
    dydt(2) = y(5);
    dydt(3) = y(6);
    
    r = norm(y(1:3));

    dydt(4) = -mu*y(1)/r^3 + (3/2)*(J2*mu*R^2)/r^4*(y(1)/r)*(5*(y(3)/r)^2-1);
    dydt(5) = -mu*y(2)/r^3 + (3/2)*(J2*mu*R^2)/r^4*(y(2)/r)*(5*(y(3)/r)^2-1);
    dydt(6) = -mu*y(3)/r^3 + (3/2)*(J2*mu*R^2)/r^4*(y(3)/r)*(5*(y(3)/r)^2-3);

end