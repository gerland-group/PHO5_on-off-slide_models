function [time_rescaling, num_interations, Flag_probs, cond_Flag_probs, log2_F_M_ratio_error, log2_F_N1_N2_ratio_error, log2_F_M_norm_ratio_error] = ...
    MP_calc_FM_error_and_time_scale_deterministic(W_gains_2, error_mode, time_scale_input, Dion_param_index)
% calculates the fit errors of measured Flag/Myc histone dynamics

    [ts, Flag_assembly_prob_fun, ts_data, log2_F_M_ratios_N1_data, err_plus, err_minus, lambda, nucl, log2_F_M_ratio_sd] = get_Flag_Myc_data_Dion(Dion_param_index);
    
    [t_data, N1_N2_F_ratio_mean, error_t_data, N1_N2_F_ratio_sd, log2_N1_N2_F_ratio_mean, log2_N1_N2_F_ratio_sd] = get_Flag_Myc_data_Rufiange();
    
    N = size(W_gains_2,3);
    
    if N > 1
        error('implementation for N > 1 not done!')
    end
    
    max_interations = 15;  
    
    FM_states = [1 1 1; 2 1 1; 1 2 1; 2 2 1; 1 1 2; 2 1 2; 1 2 2; 2 2 2; ...  % List of configurations with Flag (2) and Myc (1) nucleosomes
                 1 1 0; 2 1 0; 1 2 0; 2 2 0; ...
                 1 0 1; 2 0 1; 1 0 2; 2 0 2; ...
                 0 1 1; 0 2 1; 0 1 2; 0 2 2; ...
                 0 0 1; 0 0 2; 0 1 0; 0 2 0; 1 0 0; 2 0 0; 0 0 0;];

    M_helper = zeros(27);  % used to calculate the extended Flag/Myc rate matrix from the normal nucleosome configuration rate matrix
    for i = 1:27  % from state j to state i
        for j= 1:27
            si = FM_states(i,:);
            sj = FM_states(j,:);
            uneq = si ~= sj;
            i_uneq = find(uneq);
            if sum(uneq) == 1
                if si(uneq) == 0  % disassembly
                    M_helper(i,j) = 1;
                elseif sj(uneq) == 0  % assembly
                    if si(uneq) == 1
                        M_helper(i,j) = 3;  % new myc
                    else
                        M_helper(i,j) = 4;  % new flag
                    end
                end
            elseif (all(uneq == [1 1 0]) || all(uneq == [0 1 1])) && all(si(i_uneq([1 2])) == sj(i_uneq([2 1]))) && any([si(i_uneq(1)) sj(i_uneq(1))] == 0)  % sliding
                M_helper(i,j) = 2;
            end
        end
    end
    myc_indeces = M_helper == 3;
    flag_indeces = M_helper == 4;

    time_rescaling = zeros(N,1);
    num_interations = zeros(N,1);

    Flag_probs = zeros(length(ts)-1,3,N);  % Probability of Flag nucleosome
    cond_Flag_probs = zeros(length(ts)-1,3,N);  % Probability of Flag nucleosome conditioned on a nucleosome being there

    num_t_matches = 0;
    
    for n = 1:N

        show_progress(n,N)
        W_gain = W_gains_2(:,:,n);
        
        M_without_probs = zeros(27);  % Flag/Myc rate matrix with assembly reactions not yet multiplied with Flag/Myc assembly prob.
        for i = 1:27
            for j = 1:27
                if M_helper(i,j) > 0
                    M_without_probs(i,j) = W_gain(FM_state_to_N_state(i),FM_state_to_N_state(j));
                end
            end
        end
        
        init_FM_state = zeros(27,1);
        init_FM_state([1 9 13 17 21 23 25 27]) = calc_stat_distr(W_gain);
        
        if isfinite(time_scale_input)
            scaling = time_scale_input(n);
        else
            scaling = 1;
        end
        
        [Flag_probs_temp, cond_Flag_probs_temp, FM_dyn] = calc_Flag_Myc_probs_and_ratios_deterministic(scaling(1)*M_without_probs, ts, Flag_assembly_prob_fun, ...
        myc_indeces, flag_indeces, init_FM_state, FM_states);  % initial solutions without any rescaling
        
        if ~isfinite(time_scale_input)

            for j=2:max_interations  % rescaling iterations
                [~, i_match_value] = min(abs(cond_Flag_prob_N3_data - cond_Flag_probs_temp(:,3)));
                if num_t_matches == 4  % stop after this number of iterations after the last iteration with t_match ~= t_checkpoint
                    break
                elseif cond_Flag_probs_temp(end,3) < cond_Flag_prob_N3_data  % checkpoint value not yet reached
                    [~, i_match] = min(abs(cond_Flag_prob_N3_ode_sol - cond_Flag_probs_temp(end,3)));
                    t_match = ts(i_match);
                    if t_match == 0
                        t_match = ts(2)/2;
                    end
                    scaling(j) = scaling(j-1) * ts(end) / t_match;  % we use ts(end) because we took cond_Flag_probs_temp(end,3) to compare with the ode solution
                else  % the usual case after the checkpoint value is reached somewhere
                    t_match = ts(i_match_value);
                    if t_match == t_data
                        num_t_matches = num_t_matches + 1;  % count iterations with t_match == t_checkpoint (the last iterations to increase accuracy)
                    end
                    diff_value = cond_Flag_probs_temp(i_match_value,3) - cond_Flag_prob_N3_data;
                    if diff_value > 0 && i_match_value>=2
                        slope = (cond_Flag_probs_temp(i_match_value,3) - cond_Flag_probs_temp(i_match_value-1,3)) / ts(2);
                    else
                        slope = (cond_Flag_probs_temp(i_match_value+1,3) - cond_Flag_probs_temp(i_match_value,3)) / ts(2);
                    end
                    t_match_2 = t_match - diff_value / slope;  % correct also for difference in y values using the slope (needed especially when t_match == checkpoint_time)                  
                    scaling(j) = scaling(j-1) * t_match_2 / t_data;
                end

                [Flag_probs_temp, cond_Flag_probs_temp, FM_dyn] = calc_Flag_Myc_probs_and_ratios_deterministic(scaling(j)*M_without_probs, ts, Flag_assembly_prob_fun, ...
                myc_indeces, flag_indeces, init_FM_state, FM_states);

                num_interations = j;
            end
        end
        time_rescaling(n) = scaling(end);
        Flag_probs(:,:,n) = Flag_probs_temp(2:end,:);
        cond_Flag_probs(:,:,n) = cond_Flag_probs_temp(2:end,:);

        % TEST
%         Flag_probs_det = [sum(FM_dyn(:,FM_states(:,1) == 2),2) sum(FM_dyn(:,FM_states(:,2) == 2),2) sum(FM_dyn(:,FM_states(:,3) == 2),2)];
%         Flag_probs_given_nucl_det = Flag_probs_det ./ [sum(FM_dyn(:,FM_states(:,1) > 0),2) sum(FM_dyn(:,FM_states(:,2) > 0),2) sum(FM_dyn(:,FM_states(:,3) > 0),2)];
%     
%         T = zeros(8,27);  % to transform Flag/Myc configuration distributions into 'any nucleosome' configuration distributions
%         T(1,1:8) = 1;
%         T(2,9:12) = 1;
%         T(3,13:16) = 1;
%         T(4,17:20) = 1;
%         T(5,21:22) = 1;
%         T(6,23:24) = 1;
%         T(7,25:26) = 1;
%         T(8,27) = 1;
%         
%         figure
%         plot(ts, FM_dyn)
%         ylim([0 1])
% 
%         figure
%         subplot(4,1,1)
%         plot(ts(2:end), Flag_assembly_prob_fun(ts(2:end)))
%         ylim([0 1])
% 
%         subplot(4,1,2)
%         plot(ts, T * FM_dyn')
%         ylim([0 1])
% 
%         subplot(4,1,3)
%         plot(ts, Flag_probs_det)
%         ylim([0 1])
% 
%         subplot(4,1,4)
%         plot(ts, Flag_probs_given_nucl_det)
%         ylim([0 1])

    end
    
    log2_F_M_ratios = squeeze(log2( cond_Flag_probs(:,nucl,:) ./ (1 - cond_Flag_probs(:,nucl,:)) ));
    
    log2_F_N1_N2_ratios = squeeze(log2( Flag_probs(:,1,:) ./ Flag_probs(:,2,:) ));
    
    i_data = find(ismember(ts, ts_data));
    if length(i_data) ~= length(ts_data)
        error('not all ts_data found in ts!')
    end
    
    % in the following we assume N == 1:
    
    % not correcting for sloppy offset:
    errors_temp = log2_F_M_ratios(i_data) - log2_F_M_ratios_N1_data';
    errors_temp_plus = errors_temp(errors_temp>=0).^2 ./ err_plus(errors_temp>=0)'.^2;
    errors_temp_minus = errors_temp(errors_temp<0).^2 ./ err_minus(errors_temp<0)'.^2;
    log2_F_M_ratios_N1_squared_error_sum = sum([errors_temp_plus; errors_temp_minus]);
    
    % methods to correct the sloppy offset:
    
    % diff with mean:
    %errors_temp = (log2_F_M_ratios(i_data) - mean(log2_F_M_ratios(i_data))) - (log2_F_M_ratios_N1_data' - mean(log2_F_M_ratios_N1_data));
    %log2_F_M_ratio_diffs_N1_squared_error_mean = mean((errors_temp/sigma_F_M_ratio).^2);

    
    % diff with second data point:
    %i_norm = 2;
    %is_rest = [1 3 4];
    %errors_temp = (log2_F_M_ratios(i_data(is_rest)) - log2_F_M_ratios(i_data(i_norm))) - (log2_F_M_ratios_N1_data(is_rest)' - log2_F_M_ratios_N1_data(i_norm));
    %log2_F_M_ratio_diffs_N1_squared_error_mean = mean((errors_temp/sigma_F_M_ratio).^2);

    % log likelihood for degenerate multivariate Gaussian (diff with mean, mathematically correct):
    errors_temp = (log2_F_M_ratios(i_data) - mean(log2_F_M_ratios(i_data))) - (log2_F_M_ratios_N1_data' - mean(log2_F_M_ratios_N1_data));
    D = log2_F_M_ratio_sd*eye(4);  % initial covariance matrix
    B = eye(4) - 0.25; % linear transformation of substracting the average of the four data points
    C =  B*D*B'; % covariance matrix after the linear transformation
    log2_F_M_ratio_diffs_N1_squared_error_mean = errors_temp' * pinv(C) * errors_temp;  % using the pseudoinverse
      
    Flag_prob_ratios = squeeze(Flag_probs(:,1,:) ./ Flag_probs(:,2,:));
    
    if size(ts,1) == 1
        ts_huge = ones(size(Flag_probs,3),1)*ts;
    end
    index_data = ts == t_data;
    if ~isfinite(time_scale_input)  % calc Flag Myc fit error with second F/M data set only after finding the best time scale using the first F/M data set
 
        t_min_check = t_data-error_t_data;
        t_max_check = t_data+error_t_data;
        
        if strcmp(error_mode,'rectangle')
            t_check_start = find(ts_huge(1,2:end)>=t_min_check,1);  % model has to match the Flag ratio within this time frame (use ts of the first model, they are equal anyway)
            t_check_end = find(ts_huge(1,2:end)<=t_max_check,1,'last');
            t_check = t_check_start:t_check_end;
            combined_F_M_error = min(abs( (Flag_prob_ratios(t_check,:)-N1_N2_F_ratio_mean)/N1_N2_F_ratio_sd ),[],1);  % absolut error of ratio only
        elseif strcmp(error_mode,'ellipse')
            combined_F_M_error = min(sqrt( ((Flag_prob_ratios-N1_N2_F_ratio_mean)/N1_N2_F_ratio_sd).^2 + ...
            ((ts_huge(:,2:end)'-t_data)/error_t_data).^2 ),[],1);
        end
    else  % use given time scale and calc combined Flag Myc fit error with first and second F/M data set        
        % only Dion data (sqrt mean of squared errors of different time points)
        log2_F_M_ratio_error = sqrt(log2_F_M_ratios_N1_squared_error_sum/4);
        log2_F_M_norm_ratio_error = sqrt(log2_F_M_ratio_diffs_N1_squared_error_mean);
        
        % only Rufiange data
        %F_ratio_error = (Flag_prob_ratios(index_data,:)-N1_N2_F_ratio_mean)/N1_N2_F_ratio_sd;
        log2_F_N1_N2_ratio_error = (log2(Flag_prob_ratios(index_data,:))-log2_N1_N2_F_ratio_mean)/log2_N1_N2_F_ratio_sd;
        
    end
end