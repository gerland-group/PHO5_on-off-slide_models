function [log_LH_reduced, ps, W_gain] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, log_values, data, min_S_rate)
% calculates the likelihood, the steady state and the rate matrix with zero diagonal (W_gain)
% version without boundary treatment
  
    values = exp(log_values);  % transform to non-log values
    W_gain = zeros(size(data,1));
    values_with_zero = [0 values];
    W_gain(:) = values_with_zero(value_indeces_in_W(:)+1);  % adding 0 to values and shifting indeces by one writes zeros where no parameter sits

    % The following enforces a minimal rate > 0 for sliding reactions:
    S_glo_ind = [0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0];
    W_gain(S_glo_ind>0) = max([W_gain(S_glo_ind>0); min_S_rate*ones(1,8)],[],1);

    ps = calc_stat_distr(W_gain);
    if ~all(ps>=0)
        error('Probabilities not well-defined!')
    end
    log_LH_reduced = sum(data.*log10(ps));  % no additional minus sign here
end
