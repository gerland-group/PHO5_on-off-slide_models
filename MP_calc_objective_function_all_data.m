function [F, F_terms, log10_LHs, ps, W_gains, squ_sticky_N3_errors, occs_sticky_N3, time_rescaling, num_interations, Flag_probs, cond_Flag_probs] = ...
    MP_calc_objective_function_all_data(stage, fixed_param_pos, num_con, num_reg, value_indeces_in_W, log_values_all, fixed_param_log_value, configuration_count_data, min_S_rate, log10_LH_part2, ...
    num_prefactors, time_scale_input)
    
    %%% calc log likelihood value
    log_values = log_values_all(1:num_con+size(configuration_count_data,2)*num_reg-1);  % rate log_values without fixed parameter
    W_gains = cell(size(configuration_count_data,2),1);
    if fixed_param_pos <= num_con  % fixed param is a constitutive parameter
        % params before fixed param
        vi_1_1 = 1:fixed_param_pos - 1;
        vi_2_1 = 1:fixed_param_pos - 1;
        vi_3_1 = 1:fixed_param_pos - 1;
        % params after fixed param
        vi_1_2 = fixed_param_pos:num_con-1+num_reg;
        vi_2_2 = [fixed_param_pos:num_con-1 num_con-1+num_reg+(1:num_reg)];
        vi_3_2 = [fixed_param_pos:num_con-1 num_con-1+2*num_reg+(1:num_reg)];

        [log10_LHs(1), ps(1,:), W_gains{1}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, [log_values(vi_1_1) fixed_param_log_value log_values(vi_1_2)], configuration_count_data(:,1), min_S_rate);
        [log10_LHs(2), ps(2,:), W_gains{2}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, [log_values(vi_2_1) fixed_param_log_value log_values(vi_2_2)], configuration_count_data(:,2), min_S_rate);
        [log10_LHs(3), ps(3,:), W_gains{3}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, [log_values(vi_3_1) fixed_param_log_value log_values(vi_3_2)], configuration_count_data(:,3), min_S_rate);

    else  % fixed param is an regulated parameter (its value for the first dataset will be one)
        % params before fixed param
        vi_1_1 = 1:fixed_param_pos - 1;
        % params after fixed param
        vi_1_2 = fixed_param_pos:num_con+num_reg-1;
        vi_2 = [1:num_con num_con+num_reg-1+(1:num_reg)];
        vi_3 = [1:num_con num_con+2*num_reg-1+(1:num_reg)];

        [log10_LHs(1), ps(1,:), W_gains{1}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, [log_values(vi_1_1) fixed_param_log_value log_values(vi_1_2)], configuration_count_data(:,1), min_S_rate);
        [log10_LHs(2), ps(2,:), W_gains{2}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, log_values(vi_2), configuration_count_data(:,2), min_S_rate);
        [log10_LHs(3), ps(3,:), W_gains{3}] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, log_values(vi_3), configuration_count_data(:,3), min_S_rate);
    end
    occ_1 = [sum(ps(1,[1 3 4 5])); sum(ps(1,[1 2 4 6])); sum(ps(1,[1 2 3 7]))];
    occ_2 = [sum(ps(2,[1 3 4 5])); sum(ps(2,[1 2 4 6])); sum(ps(2,[1 2 3 7]))];
    ps = ps';
    log10_LH = sum(log10_LHs);
    
    if any(diag(W_gains{1}) ~= 0) || any(diag(W_gains{2}) ~= 0) || any(diag(W_gains{3}) ~= 0)
        W_gains{1}
        W_gains{2}
        W_gains{3}
        value_indeces_in_W
        log_values
    end
    
    %%% calc sticky N-3 fit error
    squ_sticky_N3_errors = NaN * ones(1,2);
    occs_sticky_N3 = NaN * ones(2,3);
    if stage == 2  || stage == 4  
        rate_bounds = [-inf inf];
        base_fold_change_on_Boeger_data = 0;  % default 0, then the fold change is based on the model steady state occupancy
        rescale_fix_err_with_measurement_err = 1;  % default 1

        occ_1_2 = [occ_1 occ_2];
        for mutant_id = 1:2
            if base_fold_change_on_Boeger_data
                [fix_occ_N2_N3, sd_fix_occ_N2_N3] = get_fixation_occupancy_data(occupancy_prob_data([2, 3], 1), mutant_id);  % input needs to be activated N-2 and N-3 occupancy
            else
                [fix_occ_N2_N3, sd_fix_occ_N2_N3] = get_fixation_occupancy_data(occ_1_2([2 3])', mutant_id);  % input needs to be activated N-2 and N-3 occupancy
            end

            if ~rescale_fix_err_with_measurement_err
                sd_fix_occ_N2_N3 = ones(size(sd_fix_occ_N2_N3));
            end
            
            log_prefactors = log_values_all(num_con+size(configuration_count_data,2)*num_reg-1 + ((mutant_id-1)*num_prefactors+1:mutant_id*num_prefactors));

            occs_sticky_N3(mutant_id,:) = calc_occs_after_rate_manipulation(W_gains{1}, log_prefactors, rate_bounds);

            squ_sticky_N3_errors(mutant_id) = mean( (occs_sticky_N3(mutant_id, [2 3])' - fix_occ_N2_N3).^2 ./ (sd_fix_occ_N2_N3.^2) );
        end
        sticky_N3_errors =  sqrt(squ_sticky_N3_errors);
        sticky_N3_rms_error = sqrt(mean(squ_sticky_N3_errors));
    else
        sticky_N3_errors = [0 0];
        sticky_N3_rms_error = 0;
    end
    
    %%% calc Flag/Myc ratios fit error
    if stage == 3 || stage == 4
        [time_rescaling, num_interations, Flag_probs, cond_Flag_probs, log2_F_M_ratio_error, log2_F_N1_N2_ratio_error, log2_F_M_norm_ratio_error] = ...
            MP_calc_FM_error_and_time_scale_deterministic(W_gains{2}, "ellipse", time_scale_input, 0);
    else
        time_rescaling = NaN;
        num_interations = NaN;
        ts = get_Flag_Myc_data_Dion();
        Flag_probs = NaN * ones(length(ts)-1,3);
        cond_Flag_probs = Flag_probs;
        log2_F_N1_N2_ratio_error = 0;
        log2_F_M_norm_ratio_error = 0;
    end
       
    F_terms(1) = log10_LH_part2 - log10_LH;
    F_terms(2) = log10(exp(1)) * 0.5 * 4 * sticky_N3_rms_error^2;
    F_terms(3) = log10(exp(1)) * 0.5 * log2_F_M_norm_ratio_error^2;
    F_terms(4) = log10(exp(1)) * 0.5 * log2_F_N1_N2_ratio_error^2;
    F = sum(F_terms);
end
