cclearvars;
close all;
%% Constants
% General
F = 96487; % Coulombs/mole
R = 8.314; % J/K-mole
RT_F = 25.8; % mV at 25 degrees Celsius

% Dimensions
Vol = 1e3; % um^3
S = 1e2; % um^2

% Channels
D_K = 30; % channels per um^2
D_Na = 330; % channels per um^2
lambd_K = 12; % pS
lambd_Na = 4; % pS
N_K = S * D_K; % number of channels
N_Na = S * D_Na; % number of channels
g_bar_K = 36; % mS per cm^2
g_bar_Na = 120; % mS per cm^2
g_L = 0.3; % mS per cm^2

% Electrical Properties
E_L = -49.2; % mV
Cm = 1; % uF per cm^2
Vr = -60; % mV

% Pump Currents
I_pump_NaK = 1e-18/S*1e14; % uA/cm^2
I_pump_Na = 3*I_pump_NaK; % uA/cm^2
I_pump_K = -2*I_pump_NaK; % uA/cm^2
I_PMCA = 1e-14/S*1e14; % uA/cm^2

% Simulation parameters
delta_t = 1e-2; % msec (1 delta_t = 0.01 ms)
T1 = 10e2; % msec (1e2 delta_t = 1 ms)
T2 = 50e2; % msec
T3 = 50e2; % msec
T = T1;%+T2+T3;
t = delta_t:delta_t:T;
N = length(t);

%% External Bath Concentrations
% Bath C1:
C_NaCl_C1 = 140;
C_KCl_C1 = 5.4;
C_CaCl2_C1 = 1.8;
% C_MgCl2 = 1;
% C_HEPES_C1 = 10;
% C_Glucose = 5.5;
C_Na_e_C1 = C_NaCl_C1;
C_K_e_C1 = C_KCl_C1;
C_Ca_e_C1 = C_CaCl2_C1;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

% Bath C3:
C_NaCl_C3 = 80;
C_KCl_C3 = 65.4;
C_CaCl2_C3 = 1.8;
% C_HEPES_C3 = 10;
C_Na_e_C3 = C_NaCl_C3;
C_K_e_C3 = C_KCl_C3;
C_Ca_e_C3 = C_CaCl2_C3;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

% Bath C5:
C_NaCl_C5 = 10;
C_KCl_C5 = 135.4;
C_CaCl2_C5 = 1.8;
% C_HEPES = 0;
C_Na_e_C5 = C_NaCl_C5;
C_K_e_C5 = C_KCl_C5;
C_Ca_e_C5 = C_CaCl2_C5;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

%% State Variables Initialization
% Creating the matrices
Vm = Vr*ones(1, T);
C_Na_i = 15*ones(1, T);
C_K_i = 120*ones(1, T);
C_Ca_i = 20e-6*ones(1, T); % 20 nM
n = ones(1, N);
m = ones(1, N);
h = ones(1, N);
I_K = zeros(1, N);
I_Na = zeros(1, N);
I_Ca = zeros(1, N);
I_L = zeros(1, N);
I_t = zeros(1, N);
% Init at t = 0
% Vm(1:2) = Vr;
% C_Na_i(1) = 15; % mM
% C_K_i(1) = 120; % mM
% C_Ca_i(1) = 20e-3; % mM
% C_Cl_i(1) = 10; % mM
[n(1), m(1), h(1)] = gating_vars(0, [0, 0, 0], delta_t, 1);

%% Simulation
I2C_cte = S/(Vol*F*10); %The 10 comes from dimensionality analysis
for i = 2:N-1
    if (t(i) < T1)
        C_K_e = C_K_e_C1;
        C_Na_e = C_Na_e_C1;
        C_Ca_e = C_Ca_e_C1;
    else
        if (t(i) < T1+T2)
            C_K_e = C_K_e_C3;
            C_Na_e = C_Na_e_C3;
            C_Ca_e = C_Ca_e_C3;
        else
            C_K_e = C_K_e_C5;
            C_Na_e = C_Na_e_C5;
            C_Ca_e = C_Ca_e_C5;
        end    
    end
    vm = Vm(i) - Vr;
    [n(i), m(i), h(i)] = ...
        gating_vars(vm, [n(i-1), m(i-1), h(i-1)], delta_t, 0);
    E_K = nernstVoltage(C_K_i(i), C_K_e, 1);
    E_Na = nernstVoltage(C_Na_i(i), C_Na_e, 1);
    % E_Cl = nernstVoltage(C_Cl_i, C_Cl_e, 1);
    % E_Ca = nernstVoltage(C_Ca_i, C_Ca_e, 2);
    p_K = n(i)^4;
    p_Na = m(i)^3*h(i);
    P_Ca = 1; % NEEDS UPDATING % This needs to be the same dimension.
    % Dimension of P_Ca needs to be cm/sec.
    I_K(i) = g_bar_K*p_K*(Vm(i) - E_K); % mS/cm^2*mV = uA/cm^2
    I_Na(i) = g_bar_Na*p_Na*(Vm(i) - E_Na); % mS/cm^2*mV = uA/cm^2
    x = 2*Vm*F/RT_F;
    I_Ca(i) = 2*P_Ca*x*(C_Ca_e - C_Ca_i*exp(x))/(1-exp(x));
    % OR: I_Ca(i) = I_Ca_func(...) % Pay attention to unit conversions.
    I_L(i) = g_L*(Vm(i) - E_L);
    I_t(i) = I_K(i) + I_Na(i) + I_L(i) + I_pump_NaK;
    Vm(i+1) = Vm(i) - delta_t/Cm*I_t(i);
    C_K_i(i+1) = C_K_i(i) - delta_t*(I_K(i) + I_pump_K)*I2C_cte;
    C_Na_i(i+1) = C_Na_i(i) - delta_t*(I_Na(i) + I_pump_Na)*I2C_cte;
    C_Ca_i(i+1) = Ca_Ca_i(i) - delta_t*(I_Ca(i) + I_pump_Ca)*I2C_cte/2;
    % The division by 2 is because of Ca having valence (z) = +2.
end
%% Plotting and Printing
plot(t, Vm)
xlabel("Time (ms)")
ylabel("Transmembrane Voltage (mV)")
title("V_m(t)")
figure
plot(t, C_K_i)
xlabel("Time (ms)")
ylabel("Concentrations (mM)")
title("C_i K^+ and Na^+(t)")
figure
plot(t, C_Na_i)
xlabel("Time (ms)")
ylabel("Concentrations (mM)")
title("C_i K^+ and Na^+(t)")
% hold on
% plot(t, C_Na_i)
% legend(["C_i K^+", "C_i Na^+"]);
% disp("Q28:")
% disp("I_K (A) at t1:")
% disp(I_K1*1e-15); % A
% disp("I_Na (A) at t1:")
% disp(I_Na1*1e-15); % A
% disp("I_K (A) at t2:")
% disp(I_K2*1e-15); % A
% disp("I_Na (A) at t2:")
% disp(I_Na2*1e-15); % A
