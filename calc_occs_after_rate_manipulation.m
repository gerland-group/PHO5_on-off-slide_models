function occupancy = calc_occs_after_rate_manipulation(W_gain, log_rate_factors, rate_bounds, nucls)
% calculates occupancy values resulting from rate manipulation at reactions involving the N-3 site   
% if a sliding rate is zero at the beginning, it stays zero
% new rates are cut to be within the provided rate_bounds

    W_gain_new = W_gain;
    rate_factors = exp(log_rate_factors);
    
    if length(rate_factors) == 12
        W_gain_new(4,1) = W_gain(4,1) * rate_factors(7);
        W_gain_new(1,4) = W_gain(1,4) * rate_factors(1);

        W_gain_new(6,2) = W_gain(6,2) * rate_factors(8);
        W_gain_new(2,6) = W_gain(2,6) * rate_factors(2);

        W_gain_new(5,3) = W_gain(5,3) * rate_factors(9);
        W_gain_new(3,5) = W_gain(3,5) * rate_factors(3);

        W_gain_new(8,7) = W_gain(8,7) * rate_factors(10);
        W_gain_new(7,8) = W_gain(7,8) * rate_factors(4);

        W_gain_new(4,3) = W_gain(4,3) * rate_factors(11);  % sliding
        W_gain_new(3,4) = W_gain(3,4) * rate_factors(5);  % sliding

        W_gain_new(6,7) = W_gain(6,7) * rate_factors(12);  % sliding
        W_gain_new(7,6) = W_gain(7,6) * rate_factors(6);  % sliding
    
    elseif length(rate_factors) == 4

        W_gain_new(4,1) = W_gain(4,1) * rate_factors(3);  % D3
        W_gain_new(1,4) = W_gain(1,4) * rate_factors(1);  % A3

        W_gain_new(6,2) = W_gain(6,2) * rate_factors(3);  % D3
        W_gain_new(2,6) = W_gain(2,6) * rate_factors(1);  % A3

        W_gain_new(5,3) = W_gain(5,3) * rate_factors(3);  % D3
        W_gain_new(3,5) = W_gain(3,5) * rate_factors(1);  % A3

        W_gain_new(8,7) = W_gain(8,7) * rate_factors(3);  % D3
        W_gain_new(7,8) = W_gain(7,8) * rate_factors(1);  % A3

        W_gain_new(4,3) = W_gain(4,3) * rate_factors(4);  % S32
        W_gain_new(3,4) = W_gain(3,4) * rate_factors(2);  % S23

        W_gain_new(6,7) = W_gain(6,7) * rate_factors(4);  % S32
        W_gain_new(7,6) = W_gain(7,6) * rate_factors(2);  % S23

    else
        error("number of rate_factors not 12 or 4!")
    end
    
    if any(isfinite(rate_bounds))
        W_gain_new(W_gain_new~=0) = min(max(W_gain_new(W_gain_new~=0), rate_bounds(1)), rate_bounds(2));
    end
    
    ps_temp = calc_stat_distr(W_gain_new);
    occupancy = [sum(ps_temp([1 3 4 5])); sum(ps_temp([1 2 4 6])); sum(ps_temp([1 2 3 7]))];
    
    if nargin > 3
        occupancy = occupancy(nucls);
    end
end
