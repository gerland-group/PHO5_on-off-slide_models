function [ LHs_single, ps, W_gains, fluxes, occs, net_fluxes] = MP_combined_fit_calc_obs( data, C, params_con, params_reg, values_con, values_reg, max_index, min_S_rate )
% calculates observables like fluxes from fitted parameter values

    n_data = size(data,2);  
    N = min(size(params_reg,1), max_index);
    LHs_single = zeros(N,n_data);
    ps = zeros(N,8,n_data);
    occs = zeros(N,3,n_data);
    W_gains = cell(n_data,1);
    fluxes = cell(n_data,1);
    net_fluxes = cell(n_data,1);
   
    M = zeros(n_data,1);
    for i = 1:n_data
        M(i) = multinomial(data(:,i));
    end
      
    for i=1:N
        
        params_con_temp = params_con(i,params_con(i,:)>0);  % cut off zeros at the end
        params_reg_temp = params_reg(i,params_reg(i,:)>0);
        
        value_indeces_in_W = calc_value_indeces_in_W([params_con_temp params_reg_temp], C, size(data,1), 0);
        
        num_con = length(params_con_temp);
        num_reg = length(params_reg_temp);
        
        for j=1:n_data
            values_temp = [values_con(i,1:num_con) values_reg(i,1:num_reg,j)];
            [LHs_single(i,j), ps(i,:,j), W_gains{j}(:,:,i)] = MP_calc_LH_ps_W_gain(value_indeces_in_W, log(values_temp), data(:,j), min_S_rate);
            LHs_single(i,j) = log10(M(j)) - LHs_single(i,j);
            fluxes{j}(:,:,i) = W_gains{j}(:,:,i).*(ones(1,8)'*ps(i,:,j));
            net_fluxes{j}(:,:,i) = calc_net_fluxes(fluxes{j}(:,:,i));
            ps_temp = ps(i,:,j);
            occs(i,:,j) = [sum(ps_temp([1 3 4 5])); sum(ps_temp([1 2 4 6])); sum(ps_temp([1 2 3 7]))];
        end
    end         
end
