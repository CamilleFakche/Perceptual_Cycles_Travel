%% Waves-local - Optimal Phase Shift 
% by Camille Fakche 26/03/21

dirScript = 'my_path_scripts'; % where this script is
cd(dirScript);
dirData = 'D:\Travail\2018-2022_PhD\4-Waves\Article\DataTesting' ;% 'my_path_data'; % where the data are
ListSuj = {'4wmsoci' 'egipb68' 'rrm6ne2' 'oiode78' 'ucim3ab' 'gx0xqtk' '0gygnw2' 'c71rpql' 'ss32xn3' ...
    'to23poo' 'xuimlor' 'tm75twi' 'hkqi4hi' 'gipl7pg' '69pejov' 'jvzhekl' 'psdj2b7'};
ListF = {'4','6','8','10'}; % List of frequencies 
nphasebin = 7; % number of phase bin
ListPair = {'pos1pos2','pos1pos3','pos2pos3'}; % Pair of target positions
% addpath to Circular Toolbox (Free)
% Can be downland here: https://github.com/circstat/circstat-matlab
addpath('my_path/CircStat2012a');

%% Step 1: Compute optimal phase in radians/degrees on real data averaged across participants

% Sine function: One cycle + Two cycle
fit = @(b,x)  b(1).*(sin(2*pi*x./7 + b(2))) + b(3).*(sin(2*pi*x./3.5 + b(4))) + b(5);

for F = 1:length(ListF)
    disp(['Frequency ' ListF{F} ' Hz']);
    optimal_phase_matrix = nan(3,2);
    for pos = 1:3
        % Load fitted data
        load([dirData '\subjall\subjall_data_' ListF{F} 'Hz_pos' num2str(pos) '_fitted']);
        % Use fitted parameters to create a sine wave
        x = 1:nphasebin;
        xp = linspace(min(x),max(x));
        fitted_data= fit(s,xp);
        % Identify the optimal behavioral phase
        max_performance = max(fitted_data); % lenght fit : 100
        index_max_performance = find((fit(s,xp) == max_performance));
        if length(index_max_performance) > 1
            index_max_performance = index_max_performance(1);
        end
        % Estimation of the optimal behavioral phase in radians 
        [absphase] = WavesLocal_ComputeOptimalPhaseRadians(index_max_performance);
        if pos == 1; disp('Position 3'); elseif pos == 2; disp('Position 2'); elseif pos == 3; disp('Position 1'); end
        disp(['Opimal phase (in radians): ' num2str(absphase)]);
        absphase_deg = rad2deg(absphase); % degrees
        optimal_phase_matrix(pos,2) = absphase; % absphase_deg;
        optimal_phase_matrix(pos,1) = absphase_deg;
        disp(['Opimal phase (in degrees): ' num2str(absphase_deg)]);
    end
    save([dirData  '\subjall\subjall_optimal_phase_' ListF{F} 'Hz'],'optimal_phase_matrix');
end


%% Step 2: Compute optimal phase in radians/degrees on real data for each subject

% Sine function: One cycle + Two cycle
fit = @(b,x)  b(1).*(sin(2*pi*x./7 + b(2))) + b(3).*(sin(2*pi*x./3.5 + b(4))) + b(5);

for suj = 1:length(ListSuj)
    for F = 1:length(ListF)
        optimal_phase_matrix = nan(3,2);
        for pos = 1:3
            % Load fitted data
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_data_' ListF{F} 'Hz_pos' num2str(pos) '_fitted']);
            % Use fitted parameters to create a sine wave
            x = 1:nphasebin;
            xp = linspace(min(x),max(x));
            fitted_data= fit(s,xp);
            % Identify the optimal behavioral phase
            max_performance = max(fitted_data); % lenght fit : 100
            index_max_performance = find((fit(s,xp) == max_performance));
            if length(index_max_performance) > 1
                index_max_performance = index_max_performance(1);
            end
            % Estimation of the optimal behavioral phase in radians 
            [absphase] = WavesLocal_ComputeOptimalPhaseRadians(index_max_performance);
            absphase_deg = rad2deg(absphase); % degrees
            optimal_phase_matrix(pos,2) = absphase; % absphase_deg;
            optimal_phase_matrix(pos,1) = absphase_deg;
        end
        save([dirData  '\' ListSuj{suj} '\' ListSuj{suj} '_optimal_phase_' ListF{F} 'Hz'],'optimal_phase_matrix');
    end
end

%% Step 3: Rose plot 

for F = 1:length(ListF)
    disp([ListF{F} 'Hz']); 
    % Prepare data : one vector with optimal phase for all
    % subjects, position
    % Initialize matrix
    optimal_phase_onevalue = NaN(length(ListSuj),3); % one column per pos, one line per subject
    for suj = 1:length(ListSuj)
        % Load optimal phase
        load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_optimal_phase_' ListF{F} 'Hz']);
        for pos = [1 2 3]
            % Indent matrix with absolute optimal phase in
            % radians 
            optimal_phase_onevalue(suj,pos) = optimal_phase_matrix(pos,2);
        end
    end
    figure;
    for pos = 1:3
        subplot(1,3,pos)
        circ_plot(optimal_phase_onevalue(:,pos),'hist',[],20,false,true,'linewidth',2,'color','r');
    end
    title([ListF{F} 'Hz']);
end

%% Step 4: Harrison-Kanji 

% Initialize matrix
nsuj = length(ListSuj);
optimal_phase_onevalue = []; 
freq_mat = [ones(3*nsuj,1);2*ones(3*nsuj,1);3*ones(3*nsuj,1);4*ones(3*nsuj,1)];
pos_mat = [ones(nsuj,1);2*ones(nsuj,1);3*ones(nsuj,1)];
pos_mat = repmat(pos_mat,4,1);
for F = 1:length(ListF)
    for suj = 1:nsuj
        % Load optimal phase
         load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_optimal_phase_' ListF{F} 'Hz']);
        for pos = [1 2 3]
            % Indent matrix with absolute optimal phase in radians 
            optimal_phase_onevalue = [optimal_phase_onevalue;optimal_phase_matrix(pos,2)];
        end
    end
end
optimal_phase_onevalue = [optimal_phase_onevalue,freq_mat,pos_mat];
alpha = optimal_phase_onevalue(:,1); 
idp = freq_mat; 
idq = pos_mat; 
inter = 1; fn = cell(1,2); fn{1,1} = 'Frequency'; fn{1,2} = 'Position';
[pval, stats] = circ_hktest(alpha, idp, idq, inter, fn);
disp(stats);

%% Step 5: Compute phase difference in radians for each participants, freq, pair of position

for F = 1:length(ListF)
    real_optimal_phase_difference_allsuj = nan(length(ListSuj),length(ListPair));
    for pospair = 1:length(ListPair) 
        for suj = 1:length(ListSuj)
            % Load real optimal phase
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_optimal_phase_' ListF{F} 'Hz']);
            if strcmp(ListPair{pospair},'pos1pos2')
                % Choose phase in radians (1) and compute the difference
                % for real data
                real_optimal_phase_difference = circ_dist(optimal_phase_matrix(1,2), optimal_phase_matrix(2,2));
                real_pos = 'pos3pos2';
            elseif strcmp(ListPair{pospair},'pos1pos3')
                % Choose phase in radians (1) and compute the difference
                % for real data
                real_optimal_phase_difference = circ_dist(optimal_phase_matrix(1,2), optimal_phase_matrix(3,2));
                real_pos = 'pos3pos1';
            elseif strcmp(ListPair{pospair},'pos2pos3')
                % Choose phase in radians (1) and compute the difference
                % for real data
                real_optimal_phase_difference = circ_dist(optimal_phase_matrix(2,2), optimal_phase_matrix(3,2));
                real_pos = 'pos2pos1';
            end
            % Take absolute values of difference
            real_optimal_phase_difference_abs = abs(real_optimal_phase_difference);
            % Indent matrix
            real_optimal_phase_difference_allsuj(suj,pospair) = real_optimal_phase_difference_abs;
        end
    end
    % Save
    save([dirData 'subjall_optimal_phase_difference' ListF{F} 'Hz'],'real_optimal_phase_difference_allsuj');
end

%% Step 6: Watson-William on phase difference

for F = 1:length(ListF)
    disp(ListF(F));
    % Load data
    load([dirData 'subjall_optimal_phase_difference' ListF{F} 'Hz']);
    for pospair = 1:length(ListPair)
        if strcmp(ListPair{pospair},'pos1pos2')
            real_pos = 'pos3pos2';
        elseif strcmp(ListPair{pospair},'pos1pos3')
            real_pos = 'pos3pos1';
        elseif strcmp(ListPair{pospair},'pos2pos3')
            real_pos = 'pos2pos1';
        end
        data2 = 0.00001*ones(length(ListSuj),1);
        disp(real_pos);
        % Watson-Williams test - analogue ANOVA
        % H0: Population have equal means.
        % H1: Population have unequal means.
        p1 = circ_wwtest(real_optimal_phase_difference_allsuj(:,pospair),data2);
        disp(['Watson-Williams Test : p = ' num2str(p1)]);
    end
end
