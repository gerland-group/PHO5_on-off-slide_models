function [log_LH_reduced, ps, W_gain] = MP_calc_LH_ps_W_gain(value_indeces_in_W, log_values, data, min_S_rate, rate_bounds)
% calculates the likelihood, the steady state and the rate matrix with zero diagonal (W_gain)
% version with boundary treatment and additional minus sign
    if ~exist('rate_bounds','var') || ( all(log_values>rate_bounds(1)) && all(log_values<rate_bounds(2)) )  % for combined fits, the rate bounds are enforces not here, but by fmincon
        [log_LH_reduced, ps, W_gain] = MP_calc_LH_ps_W_gain_no_bounds(value_indeces_in_W, log_values, data, min_S_rate);
        log_LH_reduced = -log_LH_reduced;  % additional minus sign here
    else  % this case should only occur in individual fitting
        log_LH_reduced = 1000;
        ps = NaN * ones(length(data),1);
        W_gain = NaN * ones(length(data));
    end

end
