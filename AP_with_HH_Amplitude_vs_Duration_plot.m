clearvars;
close all;
%% Constants
Cm = 1; % uF/cm^2
Vr = -60; % mV
delta_t = 1e-3; % ms (delta_t = 1 us)
total_T = 100; % ms
t = 0:delta_t:total_T;
% N = total_T / delta_t;
N = length(t);
% Conductances
g_bar_K = 36; % mS/cm^2
g_bar_Na = 120; % mS/cm^2
g_L = 0.3; % mS/cm^2
% Nernst Voltages
E_K = -72.1; % mV
E_Na = 52.4; % mV
E_L = -49.2; % mV
% State Variables Initialization (from table 13.3 of the book)
% Vm0 = -11.5; % mV
% n0 = 0.378;
% m0 = 0.417;
% h0 = 0.477;

%% Vector Initialization
Vm = Vr*ones(1, N);
n = ones(1, N);
m = ones(1, N);
h = ones(1, N);
[n(1), m(1), h(1)] = gating_vars(0, [0, 0, 0], delta_t, 1);
I_K = zeros(1, N);
I_Na = zeros(1, N);
I_L = zeros(1, N);
I_t = zeros(1, N);

%% Stimulation
% Step Current
step_delay = 0e1; % ms
A1 = 13.25*5; % uA/cm^2
D1 = 32e-1; % ms
T = 2e0; % ms % The time between the pulses
A2 = A1;
D2 = D1;
% Some values might result in errors due to numerical errors of storing
% numbers.
I_pulse_single1 = [zeros(1, step_delay/delta_t), ...
                  A1*ones(1, D1/delta_t),...
                  zeros(1, T/delta_t)];
I_pulse_single2 = [A2*ones(1, D2/delta_t),...
                  zeros(1, N-(step_delay+D1+D2+T)/delta_t)];
I_stim = [I_pulse_single1, I_pulse_single2];

%% Simulation
for i = 1:N-1
    vm = Vm(i) - Vr;
    p_K = n(i)^4;
    p_Na = m(i)^3*h(i);
    I_K(i) = g_bar_K*p_K*(Vm(i) - E_K); % mS/cm^2*mV = uA/cm^2
    I_Na(i) = g_bar_Na*p_Na*(Vm(i) - E_Na); % mS/cm^2*mV = uA/cm^2
    I_L(i) = g_L*(Vm(i) - E_L);
    I_t(i) = I_K(i) + I_Na(i) + I_L(i) - I_stim(i);
    delta_V = - delta_t/Cm*I_t(i);
    [n(i+1), m(i+1), h(i+1)] = ...
        gating_vars(vm, [n(i), m(i), h(i)], delta_t, 0);
    Vm(i+1) = Vm(i) + delta_V;
end

%% Plotting and Printing
% Trying to replicate the book's 5.16 figure
% g_K = g_bar_K*n.^4;
% g_Na = g_bar_Na*m.^3.*h;
% figure
% subplot(4, 1, 1);
% plot(t, I_stim)
% xlabel("Time (ms)")
% ylabel("Stimulation Current (\muA/cm^2)")
% title("I_{stim}(t)"+ " (A = "+ num2str(A1)+...
%         ", D(pulse duration) = "+ num2str(D1)+...
%         "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")
% subplot(4, 1, 2);
% plot(t, Vm);
% xlabel("Time (ms)")
% ylabel("Transmembrane Voltage (mV)")
% title("V_m(t)"+ " (A = "+ num2str(A1)+...
%         ", D(pulse duration) = "+ num2str(D1)+...
%         "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")
% subplot(4, 1, 3);
% hold on;
% plot(t, g_Na);
% plot(t, g_K, '--');
% xlabel("Time (ms)")
% ylabel("g_K and g_{Na} (mS/cm^2)")
% title("Conductances over time"+ " (A = "+ num2str(A1)+...
%         ", D(pulse duration) = "+ num2str(D1)+...
%         "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")
% legend(["g_{Na}", "g_K"]);
% subplot(4, 1, 4);
% hold on;
% plot(t, I_Na);
% plot(t, I_K, '--');
% xlabel("Time (ms)")
% ylabel("I_K and I_{Na} (\muA/cm^2)")
% title("Currents over time"+ " (A = "+ num2str(A1)+...
%         ", D(pulse duration) = "+ num2str(D1)+...
%         "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")
% legend(["I_{Na}", "I_K"]);

% The I_T Plot
AA = 13.25;
A_arr1 = [AA/4, AA/3, AA/2, AA, AA*2, AA*3, AA*4, AA*5];
T_arr1 = [13.7, 11.2, 8.9, 5.9, 3.8, 2.8, 2.2, 1.8];
DD1 = 3.2;
A_arr2 = [53, 100, 150, 200, 300, 500];
T_arr2 = [13.43, 9.82, 7.98, 6.87, 5.6, 4.45];
DD2 = 0.2;
figure
plot(T_arr1, A_arr1, '*-')
xlabel("Time between the two pulses (ms)")
ylabel("Amplitude of pulses (A) (\muA/cm^2)")
title("The relationship between pulse strength and the time required between stimulations"+...
        ", D(pulse duration) = "+ num2str(DD1)+"(ms)")
figure
plot(T_arr2, A_arr2, '*-')
xlabel("Time between the two pulses (ms)")
ylabel("Amplitude of pulses (A) (\muA/cm^2)")
title("The relationship between pulse strength and the time required between stimulations"+...
        ", D(pulse duration) = "+ num2str(DD2)+"(ms)")
figure
plot([T_arr1, T_arr2, 13.35], [DD1*A_arr1, DD2*A_arr2, 0.4*26.5], '*');
xlabel("Time between the two pulses (ms)")
ylabel("Strength of pulses (A*D) (nC/cm^2)")
title("The relationship between pulse strength and the time required between stimulations")
% legend(["I_{Na}", "I_K"]);