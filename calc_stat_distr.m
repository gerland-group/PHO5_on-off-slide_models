function ps = calc_stat_distr(W_gain)  
% wrapper for the state reduction algorithm which gains a lot of speed when compiled

    ps = calc_stat_distr_state_reduction(W_gain);
    %ps = calc_stat_distr_state_reduction_mex(W_gain);  % use MATLAB coder to build the mex function
    
end

