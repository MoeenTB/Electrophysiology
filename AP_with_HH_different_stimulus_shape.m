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
d = 0e3; % us % step_delay
A = 5; % uA/cm^2  % A = amplitude
D = 4e3; % us % D = pulse duration
tri_duration = 2*D; % A*D = pC/cm^2
T = 15e3; % us  % time Between pulses
% pulse_interval = 10e3; % 10 ms
% reps = floor((N-step_delay)/(pulse_duration+pulse_interval));
% I_step = [zeros(1, step_delay), ...
%           amplitude*ones(1, N-step_delay)];
I_pulse_single1 = [zeros(1, d), ...
                   A*ones(1, D),...
                   zeros(1, N-d-D)];
I_pulse_single2 = [zeros(1, d), ...
                   A/2*ones(1, D*2),...
                   zeros(1, N-d-D*2)];
I_pulse_single3 = [zeros(1, d), ...
                   A*2*ones(1, D/2),...
                   zeros(1, N-d-D/2)];
I_impulse = [zeros(1, d), ...
             A*D*ones(1, 1),...
             zeros(1, N-d-1)];
I_impulse2 = [zeros(1, d), ...
             A*D/2*ones(1, 2),...
             zeros(1, N-d-2)];
I_triangle = [zeros(1, d), ...
              A*triang(tri_duration)', ...
              zeros(1, N-d-tri_duration)];
% x = N - step_delay;
% The time needs changing for sawtooth and sine
I_sawtooth = [zeros(1, d), ...
              A*(1+sawtooth(t(d+1:end))')/2]; 
          % Original sawtooth is in [-1, 1] so I added 1 and divided by 2 to make it in [0, 1].
I_sin = [zeros(1, d), ...
         A*D/1000/pi*(1+sin(t(d+1:end)))/2];
I_stim = I_sin;

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
figure
subplot(2, 1, 1);
plot(t, I_stim)
xlabel("Time (ms)")
ylabel("Stimulation Current (\muA/cm^2)")
title("I_{stim}(t)"+ " (A = "+ num2str(A)+...
        ", D(pulse duration) = "+ num2str(D)+...
        "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")
subplot(2, 1, 2);
plot(t, Vm);
xlabel("Time (ms)")
ylabel("Transmembrane Voltage (mV)")
title("V_m(t)"+ " (A = "+ num2str(A)+...
        ", D(pulse duration) = "+ num2str(D)+...
        "(\mus), T(inter-pulse interval) = "+num2str(T)+ "(\mus))")