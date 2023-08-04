function [n, m, h] = gating_vars(vm, gating_vars_old, delta_t, init)
%GATING_VARS Summary of this function goes here
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
    
    n_ss = a_n./(a_n + b_n);
    m_ss = a_m./(a_m + b_m);
    h_ss = a_h./(a_h + b_h);
    n_tau = 1./(a_n + b_n);
    m_tau = 1./(a_m + b_m);
    h_tau = 1./(a_h + b_h);
    
    if (init == 0)
        n = n_ss*delta_t/n_tau - n_old*(delta_t/n_tau - 1);
        m = m_ss*delta_t/m_tau - m_old*(delta_t/m_tau - 1);
        h = h_ss*delta_t/h_tau - h_old*(delta_t/h_tau - 1);
    else
        n = n_ss;
        m = m_ss;
        h = h_ss;
    end
end

