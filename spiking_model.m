clear

number_points = 1000;
time_step = 1E-3;                           % 1 msec timestep
total_time = number_points*time_step;       % total time in seconds
time = 0:time_step: total_time-time_step;   


E_excit = 0E-03;                            % reversal potential for an excitatory synaptic input
E_inhibit = -80E-03;                        % reversal potential for an inhibitory synaptic input
Pmax = 1;

tau_synaptic = 3E-3;                        % time course of neurotranmsmitter decay in synaptic cleft (3 ms)
t = exp(tau_synaptic);                      % initial condition for state of neurotranmsmitter decay
Cortical_input_1 = zeros(number_points,1);  % initialize a vector to contain cortical inputs
stimulus = 0.25/time_step;                  % time of stimulus for cortical input
rg_1 = 2;                                   % constant from Dayan and Abbot's Synaptic Integrate and Fire Neuron.  
                                            % I will modulate this to change synaptic weight


% This loop generates a cortical "input".  It's essentially just a vector
% representing 1 electrical stimulus administered at the stimulus time.
% I just turned this into a postsynaptic synaptic probability 
% (based on an approximated transmitter diffusion).  
% I chose this model of representing synaptic transmission from Dayan and Abbott.

for p = 1:number_points
    
    if p == stimulus
        t = 0;
        Cortical_input_1(p) = (Pmax*t*exp(1-(t/tau_synaptic)))/tau_synaptic;
    else
        t = t+time_step;
        Cortical_input_1(p) = (Pmax*t*exp(1-(t/tau_synaptic)))/tau_synaptic;
    end
end


% Here I am running a number of "trials".  Each trial generates a
% random current for each neuron (to cause some spontaneous
% firing.  In particular, the pallidal neuron (cell B) fires really rapidly, at
% about 60Hz in vivo spontaneously.  This input, in addition to the cortical
% input is fed into a function (Integrate_and_Fire) that will store changes
% in membrane potential in the inhibitory interneuron (cell A), as well as
% information about its spike times.  In turn, the pallidal neuron (cell B)
% receives a combined internal input, cortical input and postsynaptic probability
% that comes from action potential firing in cell A.  Finally, the output
% of this computation can be passed to the thalamic neuron in DLM (cell C).
% I have tried to include some stochasticity in the inherent inputs.  Also,
% the synaptic strength at each input can be modulated.  Since the function
% can take a matrix for each synaptic input, I can make specifiy the
% strength of many inputs in each case, perhaps also in a randomized
% fashion to inject some variability.


trials = 50;                                    % number of trials
spikes_A = zeros(number_points, trials);        % initializing spike train and firing rate storing vectors
spikes_B = zeros(number_points, trials);
spikes_C = zeros(number_points, trials);
FR_A = zeros(trials,1);
FR_B = zeros(trials,1);
FR_C = zeros(trials,1);

C_fired = zeros(trials,1);

for trial = 1:trials

    I_internal_A = normrnd(3E-9,5E-10,number_points,1);      % inject a random current into each cell to approximate a normal 
    I_internal_B = normrnd(1E-9,0.5E-9,number_points,1);            % spontaneous firing rate in each cell
    I_internal_C = normrnd(2.3E-9,0.1E-9,number_points,1);

    [v_B, V_threshold_B, spike_B, spike_times_B, Prob_Syn_output_B] = Integrate_and_Fire(I_internal_B, Cortical_input_1, E_excit, 2, zeros(number_points,1), E_excit, rg_1, number_points, time_step);
    [v_A, V_threshold_A, spike_A, spike_times_A, Prob_Syn_output_A] = Integrate_and_Fire(I_internal_A, Cortical_input_1, E_excit, 2, Prob_Syn_output_B, E_inhibit, 4, number_points, time_step);
    [v_C, V_threshold_C, spike_C, spike_times_C, Prob_Syn_output_C] = Integrate_and_Fire(I_internal_C, zeros(number_points,1), E_excit, 1, Prob_Syn_output_A, E_inhibit, 2,  number_points, time_step);
    spikes_A(:,trial) = spike_times_A;
    spikes_B(:,trial) = spike_times_B;
    spikes_C(:,trial) = spike_times_C;
    FR_A(trial) = (sum(spike_A))/(time_step*number_points);
    FR_B(trial) = (sum(spike_B))/(time_step*number_points);
    FR_C(trial) = (sum(spike_C))/(time_step*number_points);
    
    C_fired(trial) = sum(spike_C(round(0.25/time_step):round(0.35/time_step))) >0;
end

Avg_FR_A = mean(FR_A);
Var_FR_A = var(FR_A);
Avg_FR_B = mean(FR_B);
Var_FR_B = var(FR_B);
Avg_FR_C = mean(FR_C);
Var_FR_C = var(FR_C);

Pfire = mean(C_fired);

figure(3)
clf
subplot(3,1,1)
plot (time, v_B,'b');
title(['Inhibitory Interneuron; Avg Firing Rate (Hz) = ' (num2str(Avg_FR_B))]);
xlabel('Time (sec)');
ylabel('Voltage (mV)');
axis([0 total_time -.1 .1]);
subplot(3,1,2)
plot (time, v_A,'b');
title(['Pallidal Neuron; Avg Firing Rate (Hz) = ' (num2str(Avg_FR_A))])';
xlabel('Time (sec)');
ylabel('Voltage (mV)');
axis([0 total_time -.1 .1]);
subplot(3,1,3)
plot (time, v_C,'b');
title(['DLM thalamic neuron; Avg Firing Rate (Hz) = ' (num2str(Avg_FR_C))])';
xlabel('Time (sec)');
ylabel('Voltage (mV)');
axis([0 total_time -.1 .1]);


figure(2)
clf

subplot(3,1,1)
hold on
for n=1:length(spikes_B)
    plot(spikes_B(n,:),1:trials,'*','MarkerSize', 3);
end
plot([0.25,0.25],[0,trials], 'm');
axis([time_step total_time 1 trials]);
title('Inhibitory Interneuron Raster');
xlabel('Time (sec)');
ylabel('Trial number');
hold off

subplot(3,1,2)
hold on
for n=1:length(spikes_A)
    plot(spikes_A(n,:),1:trials,'*','MarkerSize', 3);
end
plot([0.25,0.25],[0,trials], 'm');
axis([time_step total_time 1 trials]);
title('Pallidal Neuron Raster');
xlabel('Time (sec)');
ylabel('Trial number');
hold off

subplot(3,1,3)
hold on
for n=1:length(spikes_C)
    plot(spikes_C(n,:),1:trials,'*','MarkerSize', 3);
end
plot([0.25,0.25],[0,trials], 'm');
axis([time_step total_time 1 trials]);
title(['DLM Neuron Raster (P fire) = ' (num2str(Pfire))])';
xlabel('Time (sec)');
ylabel('Trial number');