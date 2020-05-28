function [Fs_fmincon, Fs, F_terms, log10_LHs, ps, W_gains, squ_sticky_N3_errors, occs_sticky_N3, values_con, values_reg, prefactors, ...
    exitflags, outputs, lambdas, time_rescaling, num_interations, Flag_probs, cond_Flag_probs, fmincon_gradients, fmincon_hessians, derivest_gradients, derivest_hessians] = ...
    MP_combined_fit_all_data(stage, LHs_data, configuration_count_data, C, params_con, params_reg, fixed_param, rate_bounds, min_S_rate, ...
    init_values_con, init_values_reg, init_prefactors, fit_time_scale, time_scale_bounds, F_max, prefactor_bounds_log, num_prefactors)

% maximum likelihood fit function used in the stages after stage 1
% stage = 2 for the "stage 2 fit" and stage = 4 for the "stage 3 fit" (stage = 3 is not used anymore in this version)

    n_data = size(configuration_count_data,2);
    if n_data ~=3
        error('only implemented for n_data == 3!')
    end
    
    if stage ~= 4
        fit_time_scale = 0;
    end
    
    N = size(params_reg,1);
    max_num_con = size(params_con,2);
    max_num_reg = size(params_reg,2);
    
    values_con = NaN * ones(N, max_num_con);
    values_reg = NaN * ones(N, max_num_reg, n_data);  % regulated values for each dataset
    prefactors = NaN * ones(N, 2*num_prefactors);
    
    Fs_fmincon = NaN * ones(N, 1);  % due to the stochastic part of optimization this can be different than Fs
    Fs = NaN * ones(N, 1);
    F_terms = NaN * ones(N, 4);
    exitflags = NaN * ones(N, 1);
    log10_LHs = zeros(N, n_data);
    ps = NaN * ones(8,n_data,N);
       
    squ_sticky_N3_errors = NaN * ones(2,1,N);
    time_rescaling = NaN * ones(N,1);
    num_interations = NaN * ones(N,1);
    
    ts = get_Flag_Myc_data_Dion();
    Flag_probs = NaN * ones(length(ts)-1,3,N);
    cond_Flag_probs = Flag_probs;
     
    M = zeros(n_data,1);
    for i = 1:n_data
        M(i) = multinomial(configuration_count_data(:,i));
    end
    log10_M = log10(M);
    log10_LH_part2 = -sum(log10_M) + sum(log10(LHs_data));

    fprintf('Optimizing parameters and prefactors for %d combinations.\n',N);

    bounds_log = log(rate_bounds);  % optimization will take place in log space
    
    if stage == 1 || stage == 3
        max_complexity = 7;  % works for maximum of 8 fit parameters (not counting prefactors here)
    else
        max_complexity = 7 + 2 * num_prefactors;
    end
    
    num_derivest_params = max_complexity+1;
    derivest_gradients = NaN * ones(N,num_derivest_params);
    derivest_hessians = NaN * ones(num_derivest_params, num_derivest_params, N);
    
    if fit_time_scale 
        num_fmincon_params = max_complexity+1;
    else
        num_fmincon_params = max_complexity;
    end
    fmincon_gradients = NaN * ones(N,num_fmincon_params);
    fmincon_hessians = NaN * ones(num_fmincon_params, num_fmincon_params, N);
    
    %fmincon_options = optimoptions('fmincon');  % uses default algorithm 'interior-point', takes at least 1.5 to 2 times longer in stage 1 without being more robust (in default setting even less accurate), however sometimes finds better LH-values
    fmincon_options = optimoptions('fmincon','Algorithm','sqp');
    fmincon_options.Display = 'none';  % Display set to iter gives progress output
    fmincon_options.ObjectiveLimit = 0;
    
    if stage == 1 || stage == 2  % options from the MP_combined_fit_cluster_setup.m script
        iteration_factor = 1e2;  % 10 is usually enough to avoid exitflag 0 when tolerace_factor = 1, 1e3 increases the needed time too much when tolerance_factor = 1e-4
        tolerance_factor = 1e-4;  % lower means lower tolerances, i.e. higher accuracy
    elseif stage == 3 || stage == 4  % slightly reduced requirements
        iteration_factor = 5;  % 2 takes up to 80min for a few models on the cluster, mean duration ~ 12min with tol. factor 1e-2
        tolerance_factor = 1e-4;  % lower means lower tolerances, i.e. higher accuracy
    else
        error("stage can only be 1, 2, 3 or 4")
    end
    
    fmincon_options.MaxIterations = iteration_factor * fmincon_options.MaxIterations;
    if strcmp(fmincon_options.Algorithm,'sqp')
        fmincon_options.MaxFunctionEvaluations = iteration_factor * 100*max_complexity;
    else
        fmincon_options.MaxFunctionEvaluations = iteration_factor * fmincon_options.MaxFunctionEvaluations;
    end
    fmincon_options.ConstraintTolerance = tolerance_factor * fmincon_options.ConstraintTolerance;
    fmincon_options.OptimalityTolerance = tolerance_factor * fmincon_options.OptimalityTolerance;
    fmincon_options.StepTolerance = tolerance_factor * fmincon_options.StepTolerance;

    %for i=1:N  % for debugging
    parfor(i=1:N, getParforArg)
        show_progress(i, N);
        
        params_con_temp = params_con(i,params_con(i,:)>0);  % cut off zeros at the end
        params_reg_temp = params_reg(i,params_reg(i,:)>0);
        
        init_values_con_temp = init_values_con(i,~isnan(init_values_con(i,:)));
        init_values_reg_temp = init_values_reg(i,~isnan(init_values_reg(i,:)));
        init_prefactors_temp = init_prefactors(i,:);
        
        % The number of constitutive and regulated parameters varies, as does the total number of parameters.
        num_con = length(params_con_temp);
        num_reg = length(params_reg_temp);
        
        if num_con ~= length(init_values_con_temp) || num_reg*n_data ~= length(init_values_reg_temp)
            error("MP_combined_fit_all_data: number of params and values do not match!")
        end
        
        init_log_values_all = log([init_values_con_temp, init_values_reg_temp]);
        
        if stage == 2
            init_log_values_all = [init_log_values_all, zeros(1,2*num_prefactors)];
        elseif stage == 4
            init_log_values_all = [init_log_values_all, log(init_prefactors_temp)];
        end
        if fit_time_scale
            init_log_values_all = [init_log_values_all, 0];
        end
        
        if stage == 1 || stage == 3
            b1 = bounds_log(1)*ones(1,num_con+n_data*num_reg-1);
            b2 = bounds_log(2)*ones(1,num_con+n_data*num_reg-1);
        else
            if stage == 4
                extend_bounds = 0;
            else
                extend_bounds = 0;
            end
            b1 = [(bounds_log(1)-extend_bounds)*ones(1,num_con+n_data*num_reg-1) prefactor_bounds_log(1)*ones(1,2*num_prefactors)];
            b2 = [(bounds_log(2)+extend_bounds)*ones(1,num_con+n_data*num_reg-1) prefactor_bounds_log(2)*ones(1,2*num_prefactors)];
        end
        if fit_time_scale
            b1 = [b1 log(time_scale_bounds(1))];
            b2 = [b2 log(time_scale_bounds(2))];
        end
        
        value_indeces_in_W = calc_value_indeces_in_W([params_con_temp params_reg_temp], C, size(configuration_count_data,1), 0);
        
        fixed_param_log_value = 0;
        fixed_param_pos = find([params_con_temp params_reg_temp params_reg_temp params_reg_temp]==fixed_param, 1);
        if isempty(fixed_param_pos)
            error('fixed_param not found in param combination!')
        end
        init_log_values_all = init_log_values_all([1:fixed_param_pos-1 fixed_param_pos+1:length(init_log_values_all)]);  % exclude fixed param
               
        if stage == 4 || stage == 3
            tic
        end
        
        if ~fit_time_scale
            f = @(log_values_all) MP_calc_objective_function_all_data(stage, fixed_param_pos, num_con, num_reg, value_indeces_in_W, log_values_all, ...
                fixed_param_log_value, configuration_count_data, min_S_rate, log10_LH_part2, num_prefactors, NaN);
        else
            f = @(log_values_all) MP_calc_objective_function_all_data(stage, fixed_param_pos, num_con, num_reg, value_indeces_in_W, log_values_all(1:end-1), ...
                fixed_param_log_value, configuration_count_data, min_S_rate, log10_LH_part2, num_prefactors, exp(log_values_all(end)));
        end
        
        [log_values_all_temp, Fs_fmincon(i), exitflags(i), outputs(i), lambdas(i), grad_temp, hessian_temp] = ...
            fmincon(f, init_log_values_all, [], [], [], [], b1, b2, [], fmincon_options );

        [Fs(i), F_terms(i,:), log10_LHs(i,:), ps(:,:,i), W_gains_temp, squ_sticky_N3_errors(:,:,i), occs_sticky_N3(:,:,i), ...
            time_rescaling(i), num_interations(i), Flag_probs(:,:,i), cond_Flag_probs(:,:,i)] = f(log_values_all_temp);
        
        % also calculate gradient and hessian with derivest matlab package (this takes time, 3 to 5 times the optimization before)
        if stage == 4 && fit_time_scale && Fs_fmincon(i) <= F_max
            
            [grad,err,finaldelta] = gradest(f, log_values_all_temp);
            [hess,err] = hessian(f, log_values_all_temp);
        
            derivest_gradients(i,:) = [grad NaN * ones(1,num_derivest_params-length(grad))];
            derivest_hessians(:,:,i) = [hess NaN * ones(size(hess,1),num_derivest_params-size(hess,2)); ...
                NaN * ones(num_derivest_params-size(hess,1), num_derivest_params);];
        end
        
        if fit_time_scale
           log_values_all_temp = log_values_all_temp(1:end-1);  % remove time scale from vector, it is saved also in time_rescaling
        end
        
        log10_LHs(i,:) = log10_LHs(i,:) + log10_M';
        
        log_values_temp = log_values_all_temp(1:num_con+n_data*num_reg-1);        
        values_con_temp = NaN*ones(1,max_num_con);
        values_reg_temp = NaN*ones(n_data,max_num_reg);
        if fixed_param_pos <= num_con  % fixed param is a constitutive parameter
            values_con_temp(1,:) = [log_values_temp(1:fixed_param_pos-1) fixed_param_log_value log_values_temp(fixed_param_pos:num_con-1) NaN*ones(1,max_num_con-num_con) ];
            values_reg_temp(:,:) = [ [log_values_temp(num_con+(0:num_reg-1)); log_values_temp(num_con+(num_reg:2*num_reg-1)); log_values_temp(num_con+2*num_reg:end)] ...
                                     NaN.*ones(n_data,max_num_reg-num_reg) ];  % first half corresponds to the act dataset, the second half to the repressed
        else  % fixed param is an regulated parameter (its value for the first dataset will be one)
            values_con_temp(1,:) = [log_values_temp(1:num_con) NaN*ones(1,max_num_con-num_con) ];
            values_reg_temp(:,:) = [ [log_values_temp(num_con+1:fixed_param_pos-1) fixed_param_log_value log_values_temp(fixed_param_pos:num_con+num_reg-1); ...
                                      log_values_temp(num_con+num_reg-1 + (1:num_reg)); log_values_temp(num_con+2*num_reg-1 + (1:num_reg))] NaN*ones(n_data,max_num_reg-num_reg) ];
        end
        values_con(i,:) = exp(values_con_temp);
        values_reg(i,:,:) = exp(values_reg_temp)';  % transpose neccessary to get indeces right
        
        if stage == 2 || stage == 4
            prefactors(i,:) = exp(log_values_all_temp(num_con+n_data*num_reg:end));
        else
            prefactors(i,:) = init_prefactors_temp;
        end
        
        W_gains_1(:,:,i) = W_gains_temp{1};
        W_gains_2(:,:,i) = W_gains_temp{2};
        W_gains_3(:,:,i) = W_gains_temp{3};
        
        fmincon_gradients(i,:) = [grad_temp' NaN * ones(1,num_fmincon_params-length(grad_temp))];
        fmincon_hessians(:,:,i) = [hessian_temp NaN * ones(size(hessian_temp,1),num_fmincon_params-size(hessian_temp,2)); ...
            NaN * ones(num_fmincon_params-size(hessian_temp,1),num_fmincon_params);];
        
        if stage == 4 || stage == 3
            toc
        end
    end
    
    W_gains{1} = W_gains_1;
    W_gains{2} = W_gains_2;
    W_gains{3} = W_gains_3;
         
    if any(exitflags~=1)
        warning('Check exitflags! Some are unequal to 1!')
    end
            
end
