%% Waves local - Fit complex sine function to data
% by Camille Fakche 26/03/21

dirScript = 'my_path_scripts'; % where this script is
cd(dirScript);
dirData = 'my_path_data'; % where the data are
ListSuj = {'4wmsoci' 'egipb68' 'rrm6ne2' 'oiode78' 'ucim3ab' 'gx0xqtk' '0gygnw2' 'c71rpql' 'ss32xn3' ...
    'to23poo' 'xuimlor' 'tm75twi' 'hkqi4hi' 'gipl7pg' '69pejov' 'jvzhekl' 'psdj2b7'};
nsuj = length(ListSuj); % number of subjects
ListF = [4, 6, 8 ,10]; % List of frequencies 
phasebin = 1:7; % List of phase bins
nphasebin = 7; % number of phase bin
ntargetpos = 3; % number of target positions 

%% Step 1: Compute hit rate for each phase bin and each subject 

for suj = 1:length(ListSuj)
    for F = 1:length(ListF)
        % Load data
        load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_result_matrix_' num2str(ListF(F)) 'Hz']);
        % result_matrix_all [nTrials * 3]
        % Column 1: Target position 1, 2, 3 (note that 1 corresponds to
        % Position 3, and 3 to Position 1)
        % Column 2: Subject's reponse 0, not perceived, 1, perceived 
        % Column 3: Phase bin of the oscillating disk at the moment of
        % target display
        for pos = 1:ntargetpos
            % Select trials with target position of interest
            index_pos = find(result_matrix_all(:,1)==pos);
            result_matrix_pos = result_matrix_all(index_pos,:);
            % Compute hit rate for each phase bin 
            data_matrix = nan(1,nphasebin);
            nTarget = nan(1,nphasebin);
            for phasebin = 1:nphasebin
                % Select row with phasebin of interest
                index_phasebin = find(result_matrix_pos(:,3)==phasebin);
                result_matrix_phasebin = result_matrix_pos(index_phasebin,:);
                ntarget = length(result_matrix_phasebin);
                indexcorr = find(result_matrix_phasebin(:,2)==1);
                indexincorr = find(result_matrix_phasebin(:,2)==0);
                ncorr = length(indexcorr); nincorr = length(indexincorr);
                hit = ncorr/ntarget;
                data_matrix(:,phasebin) = hit;
                nTarget(:,phasebin) = ntarget;
                clear hit ncorr nincorr indexcorr indexincorr ntarget result_matrix_phasebin
            end
            % Save 
            save([dirData '\' ListSuj{suj} '\'  ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'data_matrix','nTarget');
        end
    end
end

%% Step 2: Average hit rate across subjects + Compute SEM/IC95

cd(dirData);
mkdir('subjall');

for F = 1:length(ListF)
    for pos = 1:ntargetpos
        data_matrix_allsubj = NaN(nsuj,nphasebin);
        for suj = 1:length(ListSuj)
            % Load data_matrix
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
            % Indent matrix allsubj
            data_matrix_allsubj(suj,:) = data_matrix;
            clear data_matrix
        end
        % Mean across subject
        data_matrix = mean(data_matrix_allsubj,1);
        % Compute standard deviation
        STD = std(data_matrix_allsubj);
        % Compute standard error of the mean
        SEM = STD/sqrt(nsuj);
        % Compute 95 confidence interval
        IC95 = NaN(2,nphasebin);
        IC95(1,:) = data_matrix - (1.96*SEM); % Inf
        IC95(2,:) = data_matrix + (1.96*SEM); % Supp
        % Save
        save([dirData '\subjall\subjall_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'data_matrix');
        save([dirData '\subjall\subjall_SEM_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'SEM');
        save([dirData '\subjall\subjall_IC95_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'IC95');
    end
end


%% Step 3: Fit to complex sine function 

% Subj by subj or AllSubj
% ListSuj = {'subjall'};

for suj = 1:length(ListSuj)
    for F = 1:length(ListF)
        for pos = 1:ntargetpos
            % Load data
            load([dirData '\' ListSuj{suj} '\'  ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
            y = data_matrix; x = 1:7;
            yu = max(y);
            yl = min(y);
            yr = (yu-yl); % Range of ‘y’
            ym = mean(y); % Baseline of 'y'
            % Sine function: One cycle + Two cycle
            fit = @(b,x)  b(1).*(sin(2*pi*x./7 + b(2))) + b(3).*(sin(2*pi*x./3.5 + b(4))) + b(5);
            % b(1): Amplitude of the induced frequency
            % b(2): Phase offset of the induced frequency
            % b(3): Amplitude of the first harmonic
            % b(4): Phase offset of the first harmonic
            % b(5): Baseline level
            % Least-Squares cost function
            fcn = @(b) sum((fit(b,x) - y).^2);
            % Minimise Least-Squares
            s = fminsearch(fcn, [yr; -1; yr; -1; ym]);
            % Save fit parameters
            save([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos) '_fitted' ],'s');
        end
    end
end

%% Step 4: Plot real and fitted data

% Subj by subj or AllSubj
% ListSuj = {'subjall'};

fit = @(b,x)  b(1).*(sin(2*pi*x./7 + b(2))) + b(3).*(sin(2*pi*x./3.5 + b(4))) + b(5);

for suj = 1:length(ListSuj)
    disp(['Participant ' ListSuj{suj}]);
    for F = 1:length(ListF)
        disp(['Frequency ' num2str(ListF(F)) ' Hz']);
        figure;
        for pos = 1:3
            % Load fitted data
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos) '_fitted']);
            % Load real data
            load([dirData '\' ListSuj{suj} '\'  ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
            % Subplot
            subplot(1,3,pos);
            if strcmp(ListSuj{suj}, 'subjall')
                % Load IC95
                load([dirData '\subjall\subjall_IC95_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
                % Plot IC95
                if pos == 1
                    plot(x,IC95(1,:), '-b' ,'LineWidth',0.5);
                    plot(x,IC95(2,:), '-b' ,'LineWidth',0.5);
                    x2 = [x, fliplr(x)];
                    inBetween = [IC95(1,:), fliplr(IC95(2,:))];
                    fill(x2, inBetween, 'b'); 
                elseif pos == 2
                    plot(x,IC95(1,:), '-c' ,'LineWidth',0.5);
                    plot(x,IC95(2,:), '-c' ,'LineWidth',0.5);
                    x2 = [x, fliplr(x)];
                    inBetween = [IC95(1,:), fliplr(IC95(2,:))];
                    fill(x2, inBetween, 'c'); 
                elseif pos == 3
                    plot(x,IC95(1,:), '-g' ,'LineWidth',0.5);
                    plot(x,IC95(2,:), '-g' ,'LineWidth',0.5);
                    x2 = [x, fliplr(x)];
                    inBetween = [IC95(1,:), fliplr(IC95(2,:))];
                    fill(x2, inBetween, 'g'); 
                end
                alpha(0.10);
                hold on
            end
            % Plot the fit and the real data
            xp = linspace(min(x),max(x));
            if pos == 1
                scatter(1:nphasebin,data_matrix,400,'.','b'); % real data
                hold on
                plot(xp,fit(s,xp), 'b','LineWidth',2.8); % fit
                ylabel('Hit','FontSize',15);
                title('Position 3','FontSize',15);
            elseif pos == 2
                scatter(1:nphasebin,data_matrix,400,'.','c'); % real data
                hold on
                plot(xp,fit(s,xp), 'c','LineWidth',2.8); % fit
                title('Position 2','FontSize',15);
            elseif pos == 3
                scatter(1:nphasebin,data_matrix,400,'.','g'); % real data
                hold on
                plot(xp,fit(s,xp), 'g','LineWidth',2.8); % fit
                title('Position 1','FontSize',15); 
            end
            % Limits
            xlim([0 8]);
            ylim([0 1]); 
            set(gca,'XTick',1:nphasebin,'FontSize',10);
            set(gca,'XTickLabel',{'270°', ' ',' ','90°',' ',' ','270°'});
            xlabel('Phase Bin','FontSize',15);
        end
        pause
        close all
    end
end

