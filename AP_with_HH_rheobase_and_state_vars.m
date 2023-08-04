clearvars;
close all;
%% Constants
Cm = 1; % uF/cm^2
Vr = -60; % mV
delta_t = 1e-3; % ms (delta_t = 1 us)
total_T = 40; % ms
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
step_delay = 1000;
amplitude = 53; % uA/cm^2
pulse_duration = 200; % us
% pulse_interval = 10; % ms
% reps = floor((N-step_delay)/(pulse_duration+pulse_interval));
% I_step = [zeros(1, step_delay), ...
%           amplitude*ones(1, N-step_delay)];
I_pulse_single = [zeros(1, step_delay), ...
                  amplitude*ones(1, pulse_duration),...
                  zeros(1, N-step_delay-pulse_duration)];
% I_pulse_rep = [amplitude*ones(1, pulse_duration),...
%                zeros(1, pulse_interval)];
% I_pulse_rep = repmat(I_pulse_rep, 0, reps);
% I_pulse_rep = [zeros(1, step_delay),...
%                I_pulse_rep,...
%                zeros(1, mod(N-step_delay, pulse_duration+pulse_interval))];
I_stim = I_pulse_single;

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
%     [dn, dm, dh] = delta_gating_vars(vm, [n(i-1), m(i-1), h(i-1)]);
%     n(i)  = n(i-1) + delta_t * dn;
%     m(i) = m(i-1) + delta_t * dm;
%     h(i) = h(i-1) + delta_t * dh;
    Vm(i+1) = Vm(i) + delta_V;
end

%% Plotting and Printing
% Part a
[peak, t_peak] = max(Vm);
subplot(2, 1, 1);
plot(t, Vm)
xlabel("Time (ms)")
ylabel("Transmembrane Voltage (mV)")
title("V_m(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(2, 1, 2);
plot(t, I_stim)
xlabel("Time (ms)")
ylabel("Stimulation Current (mV)")
title("I_{stim}(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
    
disp("Peak = ")
disp(peak)
disp("T_peak = ")
disp(t_peak*delta_t)

% Part b
figure
subplot(4, 1, 1);
plot(t, Vm)
% xlim([0, 5]);
xlabel("Time (ms)")
ylabel("Transmembrane Voltage (mV)")
title("V_m(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(4, 1, 2);
plot(t, n)
% xlim([0, 5]);
xlabel("Time (ms)")
ylabel("n(-)")
title("n(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(4, 1, 3);
plot(t, m)
% xlim([0, 5]);
xlabel("Time (ms)")
ylabel("m(-)")
title("m(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(4, 1, 4);
plot(t, h)
% xlim([0, 5]);
xlabel("Time (ms)")
ylabel("h(-)")
title("h(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")

% Trying to replicate the book's 5.15 figure
figure
subplot(3, 1, 1);
plot(t, I_stim)
xlabel("Time (ms)")
ylabel("Stimulation Current (mV)")
title("I_{stim}(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(3, 1, 2);
plot(t, Vm);
xlabel("Time (ms)")
ylabel("Transmembrane Voltage (mV)")
title("V_m(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
subplot(3, 1, 3);
hold on;
plot(t, n);
plot(t, m);
plot(t, h);
xlabel("Time (ms)")
ylabel("gating variable (-)")
title("gating variables(t)"+ " (A = "+ num2str(amplitude)+...
        ", T = "+ num2str(total_T)+...
        "(ms), delay = "+num2str(step_delay)+ "(\mus))")
legend(["n", "m", "h"]);