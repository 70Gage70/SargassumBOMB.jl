function [alfa,tau,R,f]=parametros(delta,a,lat0)
    f = 2 * 2*pi * sin(lat0*pi/180);        % [1/d]
    rho = 1027 * 1e9;                             % [kg/km^3]
    rhoa = 1.2 * 1e9;                             % [kg/km^3]
    nu = 1e-6 * .24*.36;                          % [km^2/d]
    nua = 1.5e-5 * .24*.36;                       % [km^2/d]
    mu = nu*rho;
    mua = nua*rhoa;
    gamma = mua/mu;
    
    i=sqrt(-1);
    phi = (i*sqrt(1-(2./delta-1).^2)+2./delta-1).^(1/3);
    Phi = i*sqrt(3)/2*(1./phi-phi)-1./(2*phi)-phi/2+1;
    Phi = real(Phi);
    R = (1 - Phi/2)/(1 - Phi/6);
    Psi = 1/pi*acos(1-Phi)-1/pi*(1-Phi).*sqrt(1-(1-Phi).^2);
    
    alfa = gamma*Psi./(gamma*Psi+(1-Psi));
    tau = ((1 - Phi/6)) / ((1 + (gamma-1)*Psi)*delta^4) * a^2/(3*mu/rho);
end