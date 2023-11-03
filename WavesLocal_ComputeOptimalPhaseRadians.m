function [absphase] = WavesLocal_ComputeOptimalPhaseRadians(index_phasebin)

% Define cycle duration 
cycle_duration = 100;
% Create the corresponding sine wave with phase in radians
time = linspace(0,1,cycle_duration);
x = 2 * pi * 1 * time + -1.5; % the phase offest is selected based on the oscillatory disk 
% Compute phase in radians 
e = exp(1i*(x));
phase_cycle = angle(e);
% Transform negative angles in radians to positive
phase_cycle_bis = NaN(1,length(phase_cycle));
for i = 1:length(phase_cycle)
    if phase_cycle(i) < 0 % negative value 
        phase_cycle_bis(i) = phase_cycle(i) + (2*pi);
    else
        phase_cycle_bis(i) = phase_cycle(i);
    end
end
% Compute phase in radians according to the index of maximal phase bin
absphase = phase_cycle_bis(index_phasebin);

end





