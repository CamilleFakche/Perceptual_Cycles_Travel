%% Waves local - Monte Carlo
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
niterations = 100;%50000; % number of iterations for Monte Carlo

%% Step 1. Create 50 000 datasets by randomly assigning a performance value

for suj = 1:length(ListSuj)
    for F = 1:length(ListF) 
        for pos = 1:ntargetpos
            % Create matrix with the data of the 50 000 virtual datasets
            virtual_data_matrix_hit = NaN(niterations,nphasebin); 
            for i = 1:niterations % 50 000 datasets
                % Load result matrix [pos * corr/incorr * phasebin]
                load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_result_matrix_' num2str(ListF(F)) 'Hz']); 
                % Select trials with target position of interest
                index = find(result_matrix_all(:,1)==pos);
                result_matrix_all = result_matrix_all(index,:); 
                % Identify average performance 
                hit_pos = sum(result_matrix_all(:,2))/length(result_matrix_all)*100;
                % Clear the real results
                result_matrix_all(:,2) = NaN; 
                % Assign randomly a O (incorr) or a 1 (corr) according to
                % the average performance 
                perf = NaN(length(result_matrix_all),1);
                for n = 1:length(result_matrix_all)
                    variable = randi([0, 100]);
                    if variable > hit_pos
                        perf(n,1) = 0;
                    elseif variable < hit_pos
                        perf(n,1) = 1;
                    end
                end
                result_matrix_all(:,2) = perf;
                % Compute hit rate for each phase bin 
                Hit = nan(1,nphasebin);
                for phasebin = 1:nphasebin
                    % Select trials with phasebin of interest
                    index = find(result_matrix_all(:,3)==phasebin);
                    result_matrix_phasebin = result_matrix_all(index,:);
                    ntarget = length(result_matrix_phasebin);
                    indexcorr = find(result_matrix_phasebin(:,2)==1);
                    indexincorr = find(result_matrix_phasebin(:,2)==0);
                    ncorr = length(indexcorr); nincorr = length(indexincorr);
                    hit = ncorr/ntarget;
                    Hit(:,phasebin) = hit;
                    clear hit ncorr nincorr indexcorr indexincorr ntarget result_matrix_phasebin
                end
               % Indent matrix 
               virtual_data_matrix_hit(i,:) = Hit;
            end 
            % Save 
            virtual_data_matrix = virtual_data_matrix_hit; 
            save([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'virtual_data_matrix');
            clear virtual_data_matrix
        end
    end
end

%% Step 2. Average 50 000 virtual datasets across subjects

for F = 1:length(ListF)
    for pos = 1:ntargetpos
        % Initialize matrix
        virtual_data_matrix_allsubj = NaN(nsuj,niterations,nphasebin);
        for suj = 1:length(ListSuj)
            % Load subj virtual datasets
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
            % Indent matrix
            virtual_data_matrix_allsubj(suj,:,:) = virtual_data_matrix;
        end
        clear virtual_data_matrix
        % Mean across subjects
        virtual_data_matrix = squeeze(mean(virtual_data_matrix_allsubj(:,:,:),1));
        % Save
        save([dirData '\subjall\subjall_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos)],'virtual_data_matrix');
        clear virtual_data_matrix virtual_data_matrix_allsubj
    end
end


%% Step 3. Fit complex sine function to virtual dataset 

% Subj by subj or AllSubj
% ListSuj = {'subjall'};

for suj = 1:length(ListSuj)
    for F = 1:length(ListF)
        for pos = 1:ntargetpos
            % Load virtual data
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos)]);
            virtual_s = NaN(niterations,5);
            for iter = 1:niterations
                y = virtual_data_matrix(iter,:); x = 1:7;
                yu = max(y);
                yl = min(y);
                yr = (yu-yl); % Range of ‘y’
                ym = mean(y); % Baseline of 'y'
                % Sine function: One cycle + Two cycle
                fit = @(b,x)  b(1).*(sin(2*pi*x./7 + b(2))) + b(3).*(sin(2*pi*x./3.5 + b(4))) + b(5);
                % Least-Squares cost function
                fcn = @(b) sum((fit(b,x) - y).^2);
                % Minimise Least-Squares
                s = fminsearch(fcn, [yr; -1; yr; -1; ym]);
                virtual_s(iter,:) = s;
            end
            % Save
            save([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos) '_fitted' ],'virtual_s');
        end
    end
end

%% Step 4. Compute p-value on the amplitude parameter

% Subj by subj or AllSubj
% ListSuj = {'subjall'};

for suj = 1:length(ListSuj)
    for F = 1:length(ListF)
        pvalue_f1 = NaN(1,ntargetpos);
        pvalue_f2 = NaN(1,ntargetpos);
        pvalue_f1f2 = NaN(1,ntargetpos);
        for pos = 1:ntargetpos
            % Load real fitted data
            load([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_data_' num2str(ListF(F)) 'Hz_pos' num2str(pos) '_fitted']);
            amp_real_fitted_data_f1 = s(1);
            amp_real_fitted_data_f2 = s(3);
            % Load virtual fitted data
            load([dirData '\' ListSuj{suj} '\'  ListSuj{suj} '_virtualdata_' num2str(ListF(F)) 'Hz_pos' num2str(pos) '_fitted' ]);
            amp_virtual_fitted_data_f1 = virtual_s(:,1);
            amp_virtual_fitted_data_f2 = virtual_s(:,3);
            % Take the absolute value of amplitudes
            amp_real_fitted_data_f1_abs = abs(amp_real_fitted_data_f1);
            amp_real_fitted_data_f2_abs = abs(amp_real_fitted_data_f2);
            amp_virtual_fitted_data_f1_abs = abs(amp_virtual_fitted_data_f1);
            amp_virtual_fitted_data_f2_abs = abs(amp_virtual_fitted_data_f2);
            % Sum of both amplitude
            amp_real_fitted_data_f1f2_abs = amp_real_fitted_data_f1_abs + amp_real_fitted_data_f2_abs;
            amp_virtual_fitted_data_f1f2_abs = amp_virtual_fitted_data_f1_abs + amp_virtual_fitted_data_f2_abs;
            % Compute the proportion of surrogate fitted amplitudes = or >
            % to the fitted amplitude of the real data
            n_amp_fitted_data_supp_f1 = 0;
            for amp = 1:length(amp_virtual_fitted_data_f1_abs)
                if amp_virtual_fitted_data_f1_abs(amp) >= amp_real_fitted_data_f1_abs
                    n_amp_fitted_data_supp_f1 = n_amp_fitted_data_supp_f1 + 1;
                end
            end
            pvalue_f1(pos) = n_amp_fitted_data_supp_f1/length(amp_virtual_fitted_data_f1_abs);
            n_amp_fitted_data_supp_f2 = 0;
            for amp = 1:length(amp_virtual_fitted_data_f2_abs)
                if amp_virtual_fitted_data_f2_abs(amp) >= amp_real_fitted_data_f2_abs
                    n_amp_fitted_data_supp_f2 = n_amp_fitted_data_supp_f2 + 1;
                end
            end
            pvalue_f2(pos) = n_amp_fitted_data_supp_f2/length(amp_virtual_fitted_data_f2_abs);
            n_amp_fitted_data_supp_f1f2 = 0;
            for amp = 1:length(amp_virtual_fitted_data_f1f2_abs)
                if amp_virtual_fitted_data_f1f2_abs(amp) >= amp_real_fitted_data_f1f2_abs
                    n_amp_fitted_data_supp_f1f2 = n_amp_fitted_data_supp_f1f2 + 1;
                end
            end
            pvalue_f1f2(pos) = n_amp_fitted_data_supp_f1f2/length(amp_virtual_fitted_data_f1f2_abs);
            
        end
        disp(ListSuj{suj});
        disp([num2str(ListF(F)) 'Hz'])
        disp('F1');
        disp('Pvalue');
        disp(pvalue_f1);
        disp('F2');
        disp('Pvalue');
        disp(pvalue_f2);
        disp('F1F2');
        disp('Pvalue');
        disp(pvalue_f1f2);
        pause
        clc
        % Save
        save([dirData '\' ListSuj{suj} '\' ListSuj{suj} '_pvalue_amp_' num2str(ListF(F)) 'Hz'],'pvalue_f1','pvalue_f2','pvalue_f1f2');
    end
    
end


