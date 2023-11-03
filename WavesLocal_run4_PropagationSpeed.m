%% Waves Local - Propagation Speed
% by Camille Fakche 25/03/21
 
dirScript = 'my_path_scripts'; % where this script is
cd(dirScript);
dirData = 'my_path_data'; % where the data are
ListF = {'8','10'}; % List of frequencies 
ListPair = {'pos1pos2','pos1pos3','pos2pos3'}; % Pair of target positions

%% Step 1: Compute empirical phase difference (in degrees)

for F = 1:length(ListF)
    % Initialize phase diff matrix
    optimal_phase_diff_matrix = nan(1,3);
    disp([ListF{F} 'Hz']);
    % Load optimal phase computed on the data averaged across subjects
    disp('optimal phase computed on the data averaged across subjects');
    load([dirData  '\subjall\subjall_optimal_phase_' ListF{F} 'Hz']);
    disp(['optimal phase (in degrees): ' num2str(optimal_phase_matrix(:,1)')]);
    % Compute phase diff
    disp('optimal phase diff (in degrees)');
    for pairpos = 1:length(ListPair)
        if strcmp(ListPair{pairpos},'pos1pos2')
            optimal_phase_diff_matrix(1,pairpos) = abs(optimal_phase_matrix(3,1) - optimal_phase_matrix(2,1));
            disp([ListPair{pairpos} ': ' num2str(optimal_phase_diff_matrix(1,pairpos))]);
        elseif strcmp(ListPair{pairpos},'pos1pos3')
            optimal_phase_diff_matrix(1,pairpos) = abs(optimal_phase_matrix(3,1) - optimal_phase_matrix(1,1));
            disp([ListPair{pairpos} ': ' num2str(optimal_phase_diff_matrix(1,pairpos))]);
        elseif strcmp(ListPair{pairpos},'pos2pos3')
            optimal_phase_diff_matrix(1,pairpos) = abs(optimal_phase_matrix(2,1) - optimal_phase_matrix(1,1));
            disp([ListPair{pairpos} ': ' num2str(optimal_phase_diff_matrix(1,pairpos))]);
        end
    end
    save([dirData '\subjall\subjall_optimal_phase_diff_' ListF{F} 'Hz']);
    clear optimal_phase_matrix
end

%% Step 2: Estimate propagation speed according to the empirical phase shift (in degrees)

% Size and Distance between targets in the cortex
target_size_cortex = 0.8; % mm - diameter 
target_interval_cortex = 0.8; % mm
% Distance between targets center in the cortex
target_center_distance_cortex = (target_size_cortex/2)+target_interval_cortex+(target_size_cortex/2); % mm

for F = 1:length(ListF)
    disp([ListF{F} 'Hz']);
    load([dirData '\subjall\subjall_optimal_phase_diff_' ListF{F} 'Hz']);
    % Compute one cycle duration in ms
    time_one_cycle = 1/str2double(ListF{F})*1000;
    one_cycle = 1;
    one_cycle_degree = 360; %degrees
    % Divide by two the phase shift between position 1 & 3 because there is
    % a distance of two targets
    phase_shift_degree = optimal_phase_diff_matrix;
    phase_shift_degree(2) = phase_shift_degree(2)./2;
    average_phase_shift_degree = mean(phase_shift_degree);
    % Convert empirical phase shift in degree to cycle distance shift
    cycle_distance_shift = average_phase_shift_degree*one_cycle/one_cycle_degree;
    % Convert cycle distance shift to time lag
    time_lag = time_one_cycle*cycle_distance_shift/one_cycle; % ms
    % Compute propagation speed
    propagation_speed_mm_target_center = target_center_distance_cortex/time_lag;
    disp(['Propagation speed target center: ' num2str(propagation_speed_mm_target_center) ' m/s']);
end




