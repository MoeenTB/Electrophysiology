function[I_Ca,E_Ca,X_HH_Ca,tauCa,gCa_max,gCa] = I_Ca_func(dt,F,R,T,cCai,cCao,S,V,X_HH_Ca,Vr)
    gCa_max=(0.01/1e+12)*S; % g:Siemens , S = Surface
    E_Ca=R*T/(2*F)*log(cCao/cCai); % v
    vm=1000*(V-Vr); % all v in volts - but for the following formulas vm shall be in mv
    alpha(1)=(0.055*(-27-(vm)))/(exp((-27-(vm))/3.8) - 1);
    alpha(2)=(0.000457*exp((-13-(vm))/50));
    beta(1)=(0.94*exp((-75-(vm))/17));
    beta(2)=(0.0065/(exp((-(vm)-15)/28)+1));
    tauCa=1./(alpha+beta); % ms 
    x_inf=alpha.*tauCa;
    X_HH_Ca=(1-dt./tauCa).*X_HH_Ca+dt./tauCa.*x_inf;  % dt is in ms
    gCa=gCa_max*X_HH_Ca(1)^2*X_HH_Ca(2);
    I_Ca=gCa*(V-E_Ca); % A
end