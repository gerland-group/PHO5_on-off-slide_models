function [ LHs, values_con, values_reg, exitflags ] = MP_combined_fit( data, C, params_con, params_reg, fixed_param, rate_bounds, min_S_rate, opt_function_options)
% First stage of the fit procedure. Calculates the maximum likelihood and the parameter values (process rates). 
% Inital rates are random thus several runs are recommended.

    n_data = size(data,2);    
    N = size(params_reg,1);
    max_num_con = size(params_con,2);
    max_num_reg = size(params_reg,2);
    values_con = NaN * ones(N, max_num_con);
    values_reg = NaN * ones(N, max_num_reg, n_data);  % regulated values for each dataset
    LHs = zeros(N,1);
    exitflags = zeros(N,1);
   
    M = zeros(n_data,1);
    for i = 1:n_data
        M(i) = multinomial(data(:,i));
    end

    fprintf('Optimizing parameters for %d combinations.\n',N);

    bounds_log = log(rate_bounds);  % optimization will take place in log space

    %for i=1:N  % for debugging
    parfor(i=1:N, getParforArg)  % runs a normal for loop if no parallel pool is activated

        show_progress(i, N);
        
        params_con_temp = params_con(i,params_con(i,:)>0);  % cut off zeros at the end
        params_reg_temp = params_reg(i,params_reg(i,:)>0);
        
        % The number of constitutive and regulated parameters varies, as does the total number of parameters.
        num_con = length(params_con_temp);
        num_reg = length(params_reg_temp);
        
        init_values = bounds_log(1)+rand(1,num_con+n_data*num_reg-1)*(bounds_log(2)-bounds_log(1));  % initial values equally distributed in log space between rate_bounds
        %init_values = [-4.1200    2.1908   -2.1265   -0.7107    0.4409    4.0778];  % random initial values for reproducable testing
        
        b1 = bounds_log(1)*ones(1,num_con+n_data*num_reg-1);
        b2 = bounds_log(2)*ones(1,num_con+n_data*num_reg-1);
        
        value_indeces_in_W = calc_value_indeces_in_W([params_con_temp params_reg_temp], C, size(data,1), 0);
        
        % The value for the fixed param will be set to 1 -> one parameter less.
        fixed_param_value = log(1);
        fixed_param_pos = find([params_con_temp params_reg_temp]==fixed_param, 1);
        if isempty(fixed_param_pos)
            error('fixed_param not found in param combination!')
        end
             
        values_con_temp = NaN*ones(1,max_num_con);
        values_reg_temp = NaN*ones(n_data,max_num_reg);
        
        switch n_data
            
            case 2
                if fixed_param_pos <= num_con  % fixed param is a constitutive parameter
                    % params before fixed param
                    vi_1_1 = 1:fixed_param_pos - 1;
                    vi_2_1 = 1:fixed_param_pos - 1;
                    % params after fixed param
                    vi_1_2 = fixed_param_pos:num_con-1+num_reg;
                    vi_2_2 = [fixed_param_pos:num_con-1 num_con-1+num_reg+(1:num_reg)];

                    [values_temp, LHs(i), exitflags(i), output, lambda, grad, hessian] = fmincon(@(values) MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_1_1) fixed_param_value values(vi_1_2)], data(:,1), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_2_1) fixed_param_value values(vi_2_2)], data(:,2), min_S_rate), init_values, [], [], [], [], b1, b2, [], opt_function_options );

                    values_con_temp(1,:) = [values_temp(1:fixed_param_pos-1) fixed_param_value values_temp(fixed_param_pos:num_con-1) NaN*ones(1,max_num_con-num_con) ];
                    values_reg_temp(:,:) = [ [values_temp(num_con+(0:num_reg-1)); values_temp(num_con+num_reg:end)] NaN.*ones(n_data,max_num_reg-num_reg) ];  % first half corresponds to the act dataset, the second half to the repressed

                else  % fixed param is an regulated parameter (its value for the first dataset will be one)
                    % params before fixed param
                    vi_1_1 = 1:fixed_param_pos - 1;
                    % params after fixed param
                    vi_1_2 = fixed_param_pos:num_con+num_reg-1;
                    vi_2 = [1:num_con num_con+num_reg-1+(1:num_reg)];

                    init_values(fixed_param_pos+num_reg-1) = 1;  % set the initial value of the previously fixed parameter to 1 to avoid bad numerics

                    [values_temp, LHs(i), exitflags(i), output, lambda, grad, hessian] = fmincon(@(values) MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_1_1) fixed_param_value values(vi_1_2)], data(:,1), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, values(vi_2), data(:,2), min_S_rate), init_values, [], [], [], [], b1, b2, [], opt_function_options );

                    values_con_temp(1,:) = [values_temp(1:num_con) NaN*ones(1,max_num_con-num_con) ];
                    values_reg_temp(:,:) = [ [values_temp(num_con+1:fixed_param_pos-1) fixed_param_value values_temp(fixed_param_pos:num_con+num_reg-1); values_temp(num_con+num_reg:end)] NaN*ones(n_data,max_num_reg-num_reg) ];  % first half corresponds to the act dataset, the second half to the repressed
                end

                values_con(i,:) = exp(values_con_temp);
                values_reg(i,:,:) = exp(values_reg_temp)';  % transpose neccessary to get indeces right
                
            case 3
                 if fixed_param_pos <= num_con  % fixed param is a constitutive parameter
                    % params before fixed param
                    vi_1_1 = 1:fixed_param_pos - 1;
                    vi_2_1 = 1:fixed_param_pos - 1;
                    vi_3_1 = 1:fixed_param_pos - 1;
                    % params after fixed param
                    vi_1_2 = fixed_param_pos:num_con-1+num_reg;
                    vi_2_2 = [fixed_param_pos:num_con-1 num_con-1+num_reg+(1:num_reg)];
                    vi_3_2 = [fixed_param_pos:num_con-1 num_con-1+2*num_reg+(1:num_reg)];

                    [values_temp, LHs(i), exitflags(i), output, lambda, grad, hessian] = fmincon(@(values) MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_1_1) fixed_param_value values(vi_1_2)], data(:,1), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_2_1) fixed_param_value values(vi_2_2)], data(:,2), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_3_1) fixed_param_value values(vi_3_2)], data(:,3), min_S_rate), init_values, [], [], [], [], b1, b2, [], opt_function_options );

                    values_con_temp(1,:) = [values_temp(1:fixed_param_pos-1) fixed_param_value values_temp(fixed_param_pos:num_con-1) NaN*ones(1,max_num_con-num_con) ];
                    values_reg_temp(:,:) = [ [values_temp(num_con+(0:num_reg-1)); values_temp(num_con+(num_reg:2*num_reg-1)); values_temp(num_con+2*num_reg:end)] NaN.*ones(n_data,max_num_reg-num_reg) ];  % first half corresponds to the act dataset, the second half to the repressed

                else  % fixed param is an regulated parameter (its value for the first dataset will be one)
                    % params before fixed param
                    vi_1_1 = 1:fixed_param_pos - 1;
                    % params after fixed param
                    vi_1_2 = fixed_param_pos:num_con+num_reg-1;
                    vi_2 = [1:num_con num_con+num_reg-1+(1:num_reg)];
                    vi_3 = [1:num_con num_con+2*num_reg-1+(1:num_reg)];

                    
                    init_values(fixed_param_pos+num_reg-1) = 1;  % set the initial value of the previously fixed parameter to 1 to avoid bad numerics
                    init_values(fixed_param_pos+2*num_reg-1) = 1;  % set the initial value of the previously fixed parameter to 1 to avoid bad numerics
                    
                    [values_temp, LHs(i), exitflags(i), output, lambda, grad, hessian] = fmincon(@(values) MP_calc_LH_ps_W_gain(value_indeces_in_W, [values(vi_1_1) fixed_param_value values(vi_1_2)], data(:,1), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, values(vi_2), data(:,2), min_S_rate) + MP_calc_LH_ps_W_gain(value_indeces_in_W, values(vi_3), data(:,3), min_S_rate), init_values, [], [], [], [], b1, b2, [], opt_function_options );

                    values_con_temp(1,:) = [values_temp(1:num_con) NaN*ones(1,max_num_con-num_con) ];
                    values_reg_temp(:,:) = [ [values_temp(num_con+1:fixed_param_pos-1) fixed_param_value values_temp(fixed_param_pos:num_con+num_reg-1); values_temp(num_con+num_reg-1 + (1:num_reg)); values_temp(num_con+2*num_reg-1 + (1:num_reg))] NaN*ones(n_data,max_num_reg-num_reg) ];
                end

                values_con(i,:) = exp(values_con_temp);
                values_reg(i,:,:) = exp(values_reg_temp)';  % transpose neccessary to get indeces right
                
            case 4
                error('Fitting 4 data sets not yet implemented!')
                
        end
        
        % for testing:
%         params_con_temp
%         params_reg_temp
%         values_temp
%         grad
%         hessian
%         [V, D] = eig(hessian)
%         min(D(D~=0))
        
    end
    
    LHs = -LHs + sum(log10(M));
    
    if any(exitflags~=1)
        warning('Check exitflags! Some are unequal to 1!')
    end
            
end
