function [ts, Flag_assembly_prob_fun, ts_data, log2_F_M_ratios_N1_data, err_plus, err_minus, lambda, nucl, log2_F_M_ratio_sd] = get_Flag_Myc_data_Dion(param_index)
    
    ts_data = [1 1.5 2 2.5] - 0.25;
    
    if ~exist('param_index','var')
        param_index = 0;
    end
    
    nucl = 1;  % probe at N-1 position
    lambda = 0.0230754*60;  % inferred turnover rate in 1/h at N-1 (Chr2:431049-431108, from Suppl. Table S2 Dion2007)
    log2_F_M_ratios_N1_data = [-0.416667165399388, 1.23527284184153, 1.86857649796473, 2.59602638116156];  % after redoing the fit to find the normalization constants
    err_plus = 0.85 * ones(1,4);
    err_minus = 1.3 * ones(1,4);
    
    log2_F_M_ratio_sd = 0.4;  % estimated standard deviation of measurement in log2, corresponding std in Rufiange et al. data is 0.2043
    
    ts = 0:0.002:3;  % used to calculate Flag Myc results
    p_max = 0.9375;  % limit of the probability that a new nucleosome is Flag tagged (from personal correspondence, Dion2007)
    beta_F = 0.01*60;  % Flag degradation rate in 1/h
    
    if param_index == 0  % default
        Flag_assembly_prob_fun = @(t) 1./(1 + (1/p_max-1)./(1-exp(-beta_F*t)));  % probability that a new nucleosome is going to be Flag tagged
    else  % used for testing other nucleosome pool parameters (not needed anymore)
        p_factors = [1];
        beta_F_factors = [10 1/10];

        k = 1;
        for i=1:length(p_factors)
            for j=1:length(beta_F_factors)
                q_setup(k,:) = [p_factors(i) beta_F_factors(j)];
                k = k+1;
            end
        end
        p_setup = q_setup .* (ones(size(q_setup,1),1) * [p_max beta_F]);
        Flag_assembly_prob_fun = @(t) 1./(1 + (1/p_setup(param_index,1)-1)./(1-exp(-p_setup(param_index,2)*t)));
    end
    
end
