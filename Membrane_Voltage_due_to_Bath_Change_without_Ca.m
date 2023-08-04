clearvars;
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
% the N_p*lambd_p was very close to g_bar_p for both ions.

% Electrical Properties
E_L = -49.2; % mV
Cm = 1; % uF per cm^2
Vr = -60; % mV

% Pump Currents
I_pump_NaK = 1e-18/S*1e14; % uA/cm^2
I_pump_Na = 3*I_pump_NaK; % uA/cm^2
I_pump_K = -2*I_pump_NaK; % uA/cm^2
% I_PMCA = 1e-14; % uA/cm^2

% Simulation parameters
delta_t = 1e-2; % msec (1 delta_t = 0.01 ms = 10 micros)
T1 = 5e1; % msec
T2 = 5e1; % msec
T3 = 5e1; % msec
T = T1+T2+T3;
t = delta_t:delta_t:T;
N = length(t);

%% External Bath Concentrations
% Bath C1:
C_NaCl_C1 = 140;
C_KCl_C1 = 5.4;
% C_CaCl2_C1 = 1.8;
% C_MgCl2 = 1;
% C_HEPES_C1 = 10;
% C_Glucose = 5.5;
C_Na_e_C1 = C_NaCl_C1;
C_K_e_C1 = C_KCl_C1;
% C_Ca_e_C1 = C_CaCl2_C1;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

% Bath C3:
C_NaCl_C3 = 80;
C_KCl_C3 = 65.4;
% C_CaCl2_C3 = 1.8;
% C_HEPES_C3 = 10;
C_Na_e_C3 = C_NaCl_C3;
C_K_e_C3 = C_KCl_C3;
% C_Ca_e_C3 = C_CaCl2_C3;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

% Bath C5:
C_NaCl_C5 = 10;
C_KCl_C5 = 135.4;
% C_CaCl2_C5 = 1.8;
% C_HEPES = 0;
C_Na_e_C5 = C_NaCl_C5;
C_K_e_C5 = C_KCl_C5;
% C_Ca_e_C5 = C_CaCl2_C5;
% C_Cl_e = C_NaCl+C_KCl+2*C_CaCl2+2*C_MgCl2;
% Osmolarity = C_NaCl+C_KCl+C_CaCl2+2*C_MgCl2;

%% State (and other) Variables Initialization
% Creating the matrices
Vm = Vr*ones(1, T);
C_Na_i = 15*ones(1, T); % mM
C_K_i = 120*ones(1, T); % mM
n = ones(1, N);
m = ones(1, N);
h = ones(1, N);
I_K = zeros(1, N);
I_Na = zeros(1, N);
I_L = zeros(1, N);
I_t = zeros(1, N);
E_K = zeros(1, N);
E_Na = zeros(1, N);
% Init at t = 0
% Vm(1:2) = Vr;
% C_Na_i(1) = 15; % mM
% C_K_i(1) = 120; % mM
% C_Ca_i(1) = 20e-3; % mM
% C_Cl_i(1) = 10; % mM
[n(1), m(1), h(1)] = gating_vars(0, [0, 0, 0], delta_t, 1);

%% Simulation
I2C_cte = S/(Vol*F*10); %The 10 comes from dimensionality analysis
for i = 1:N-1
    if (t(i) < T1)
        C_K_e = C_K_e_C1;
        C_Na_e = C_Na_e_C1;
    else
        if (t(i) < T1+T2)
            C_K_e = C_K_e_C3;
            C_Na_e = C_Na_e_C3;
        else
            C_K_e = C_K_e_C5;
            C_Na_e = C_Na_e_C5;
        end    
    end
    vm = Vm(i) - Vr;
    E_K(i) = nernstVoltage(C_K_i(i), C_K_e, 1);
    E_Na(i) = nernstVoltage(C_Na_i(i), C_Na_e, 1);
    p_K = n(i)^4;
    p_Na = m(i)^3*h(i);
    I_K(i) = g_bar_K*p_K*(Vm(i) - E_K(i)); % mS/cm^2*mV = uA/cm^2
    I_Na(i) = g_bar_Na*p_Na*(Vm(i) - E_Na(i)); % mS/cm^2*mV = uA/cm^2
    I_L(i) = g_L*(Vm(i) - E_L);
    I_t(i) = I_K(i) + I_Na(i) + I_L(i) + I_pump_NaK;
    delta_V = - delta_t/Cm*I_t(i); % V(i+1) = V(i) + delta_V(i)
    [n(i+1), m(i+1), h(i+1)] = ...
        gating_vars(vm, [n(i), m(i), h(i)], delta_t, 0);
    Vm(i+1) = Vm(i) + delta_V;
    C_K_i(i+1) = C_K_i(i) - delta_t*(I_K(i) + I_pump_K)*I2C_cte;
    C_Na_i(i+1) = C_Na_i(i) - delta_t*(I_Na(i) + I_pump_Na)*I2C_cte;
end
%% Some Calculations
I_K_diff = abs(I_K - I_pump_K);
I_K_to_I_Na_rat = I_K./I_Na;
% V_r_GHK =  % Need relative permeabilities
g_K = g_bar_K*n.^4;
g_Na = g_bar_Na*m.^3.*h;
Vr_HH = (g_K.*E_K + g_Na.*E_Na + g_L*E_L - I_pump_NaK)/...
        (g_K + g_Na + g_L);
%% Plotting and Printing
% Vm
figure
plot(t, Vm)
hold on
plot(T1*ones(1, 20), linspace(min(Vm), max(Vm), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(Vm), max(Vm), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("V_m (mV)")
title("Transmembrane Voltage vs. time")
% C_K_i
figure
plot(t, C_K_i)
hold on
plot(T1*ones(1, 20), linspace(min(C_K_i), max(C_K_i), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(C_K_i), max(C_K_i), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("[K^+]_i(mM)")
title("Internal Potassium Concentration vs. time")
% C_Na_i
figure
plot(t, C_Na_i)
hold on
plot(T1*ones(1, 20), linspace(min(C_Na_i), max(C_Na_i), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(C_Na_i), max(C_Na_i), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("[Na^+]_i(mM)")
title("Internal Sodium Concentration vs. time")
% figure
% plot(t, C_Na_i)
% xlabel("Time (ms)")
% ylabel("Concentrations (mM)")
% title("C_i K^+ and Na^+(t)")
% hold on
% plot(t, C_Na_i)
% legend(["C_i K^+", "C_i Na^+"]);
% Vm, Vr difference
figure
plot(t, abs(Vm-Vr_HH))
hold on
plot(T1*ones(1, 20), linspace(min(abs(Vm-Vr_HH)), max(abs(Vm-Vr_HH)), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(abs(Vm-Vr_HH)), max(abs(Vm-Vr_HH)), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("|V_m - V_{rest, HH}| (mV)")
title("Difference of V_m and the V_{rest} from HH model over time")
% I_K to I_Na ratio
figure
plot(t, I_K_to_I_Na_rat)
hold on
plot(T1*ones(1, 20), linspace(min(I_K_to_I_Na_rat), max(I_K_to_I_Na_rat), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(I_K_to_I_Na_rat), max(I_K_to_I_Na_rat), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("I_K / I_Na (mV)")
title("Ratio of I_K to I_{Na} over time")
% I_K
figure
plot(t, I_K)
hold on
plot(T1*ones(1, 20), linspace(min(I_K), max(I_K), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(I_K), max(I_K), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("I_K (\muA/cm^2)")
title("Transmembrane Potassium Current vs. time")
disp("I_pump_K (uA/cm^2) = ");
disp(I_pump_K);
% I_Na
figure
plot(t, I_Na)
hold on
plot(T1*ones(1, 20), linspace(min(I_Na), max(I_Na), 20), 'r-', 'LineWidth', 0.1)
plot((T1+T2)*ones(1, 20), linspace(min(I_Na), max(I_Na), 20), 'r-', 'LineWidth', 0.1)
xlabel("Time (ms)")
ylabel("I_Na (\muA/cm^2)")
title("Transmembrane Sodium Current vs. time")
disp("I_pump_Na (uA/cm^2) = ");
disp(I_pump_Na);