function [V, V_threshold, spike, spike_times, Prob_Syn_output] = Integrate_and_Fire(I, Ps_1, Es_1, rg_1, Ps_2, Es_2, rg_2, number_points, time_step)

%Integrate_and_Fire takes an internal input, and up to 2 other general
%classes of inputs (although there could me many inputs within these
%classes) including information about their polarity and strength.  It will
%determine changes in the membrane potential for a unit receiving these
%inputs, over a time period determined by number_points and time_step
%according to a classical integrate and fire neuron model.  For each unit
%it will output:
    % vector V containing the resulting membrane potential of
        % that cell
    % vector V_threshold showing the threshold for firing of that cell (a
        % little trick I used to impose a somewhat arbitrary refractory period)
    % vector spike containing a binary 0 or 1 for an action potential at 
        % each time point
    % vector spike_times, listing the time (in seconds) at which one of
        % these spikes occured
    % vector Prob_Syn_output.  This was another computation that I thought 
        % was useful while looking through the Dayan and Abbott text.  It's 
        % essentially a way of modeling presence of neurotransmitter in the
        % cleft after an AP, which will clearly affect the probability of 
        % postsynaptic events, something in which I am interested in since
        % I'm passing the output of this function as the input into a
        % another cell.
        
% I'm not sure how physiological some of the constants I've chosen to use
% are, but I believe that mostly these would affect the strength of the
% inputs that need to be injected.

V = zeros(number_points,1);             % Initialize what will be my output variables
V_threshold = zeros(number_points,1);
spike = zeros(number_points,1);
spike_times = zeros(number_points,1);

V_rest = -65E-3;                        % Resting membrane potential = -65 mV
V_reset = -80E-3;                       % After an action potential, membrane potential is reset to - 80mV
AP = 100E-3;                            % I made a spike go up to 100 mV
tau = 10E-3;                            % Membrane time constant = 10ms
R_mem = 10E6;                           % Membrane resistance = 10E6
min_V_threshold = -45E-3;               % V_threshold for spike generation = -45 mV.  It can never be less than this.
max_V_threshold = 200E-3;               % After an AP is fired, I reset the threshold voltage for firing another AP to 200 mV
                                            % to simulate a refractory
                                            % period.

abs_refractory_period = 2E-3;           % The length of my refractory period, or the amount of time before V_threshold returns to
                                            % its minimum value (-45mV)
refractory_delta = ((max_V_threshold - min_V_threshold)/(abs_refractory_period))*time_step;         % How much V_threshold changes 
                                                                                                        % with each time step after
                                                                                                        % an action potential
                                          

V(1) = V_rest;                          % Initial membrane potential
V_threshold(1) = min_V_threshold;       % V_threshold starts at its normal (minimal) value

    
Pmax = 1;                               % Maximum probability of postsynaptic conductance
tau_synaptic = 3E-3;                   % Time course of neurotransmitter decay, or probability of postsynaptic conductance
Prob_Syn_ouput = zeros(number_points,1);% Initialize postsynaptic probability    
Prob_Syn_ouput(1) = 0;
t=exp(tau_synaptic);                    % Initialize t

syn_input = 1;%length();                 % Number of synaptic inputs of a particular class
PSP = zeros(syn_input, 1);              % A vector to store the membrane potential changes resulting from each input

k=1;
for n = 2:number_points
    V_threshold(n) = max((V_threshold(n-1) - refractory_delta), min_V_threshold);
    if V(n-1) == AP;
        V(n) = V_reset;
        spike(n)=0;
        t=0;
        Prob_Syn_output(n) = (Pmax*t*exp(1-(t/tau_synaptic)))/tau_synaptic;
    elseif V(n-1) >= V_threshold(n-1)
        V(n) = AP;
        V_threshold(n) = max_V_threshold;
        spike(n) =1;
        spike_times(k) = n*time_step;
        k = k+1;
        t = t+time_step;
        Prob_Syn_output(n) = (Pmax*t*exp(1-(t/tau_synaptic)))/tau_synaptic;
    else
        for i = 1:syn_input
            PSP(i) = min(V(n-1) + time_step*((V_rest-V(n-1)-((rg_1(i))*(Ps_1(n-1))*(V(n-1)-Es_1(i)))-((rg_2(i))*(Ps_2(n-1))*(V(n-1)-Es_2(i)))+(R_mem*(I(n-1))))/tau),V_threshold(n));
        end
        V(n) = sum(PSP);
        spike(n)=0;
        t = t+time_step;
        Prob_Syn_output(n) = (Pmax*t*exp(1-(t/tau_synaptic)))/tau_synaptic;
    end
end;

end

