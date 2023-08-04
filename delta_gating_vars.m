function [dn, dm, dh] = delta_gating_vars(vm, gating_vars_old)
%DELTA_GATING_VARS delta_n without delta_t
%   Detailed explanation goes here
    n_old  = gating_vars_old(1);
    m_old = gating_vars_old(2);
    h_old = gating_vars_old(3);
    
    a_n = .01*(10-vm)/(exp((10-vm)/10)-1);
    b_n = .125*exp(-vm/80);
    a_m = .1*(25-vm)/(exp(.1*(25-vm))-1);
    b_m = 4*exp(-vm/18);
    a_h = .07*exp(-vm/20);
    b_h = 1./(exp((30-vm)/10)+1);
    
    dn = a_n * (1-n_old) + b_n * n_old;
    dm = a_m * (1-m_old) + b_m * m_old;
    dh = a_h * (1-h_old) + b_h * h_old;
end

