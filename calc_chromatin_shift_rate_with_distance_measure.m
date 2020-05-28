function [shift_rate, shift_rate_with_distance_measure, ts_shift, shift_probs] = calc_chromatin_shift_rate_with_distance_measure(W_gain, ps_start, plot_flag)
% alternative estimate of the effective chromatin opening/closing rate using distance measures between configuration distributions


    my_ode_options = odeset('NonNegative',1,'MaxStep',0.01);
    calc_distance = @(x,y) norm(sqrt(x) - sqrt(y)) / sqrt(2);  % Hellinger distance
    %calc_distance = @(x,y) norm(x - y);  % Euclidean distance

    ps_end = calc_stat_distr(W_gain)';
    
    [shift_rate, W] = calc_chromatin_shift_rate(W_gain);
    
    ts_shift = (0:0.01:5).*(1/shift_rate);
    [~, shift_probs] = ode23(@(t,y) W*y, ts_shift, ps_start, my_ode_options);
    distances = NaN * ones(size(ts_shift));
    for l=1:length(ts_shift)
        distances(l) = calc_distance(shift_probs(l,:), ps_end);
    end
    myfit = fit(ts_shift', distances', fittype('a*exp(-b*x)'), 'StartPoint', [1 1]);
    
    shift_rate_with_distance_measure = myfit.b;
    
    if plot_flag
        plot(ts_shift, shift_probs) 
        hold on
        plot(ts_shift,distances,'k--')
        plot(ts_shift, myfit(ts_shift), 'k:')
        xlim([min(ts_shift) max(ts_shift)])
        xlabel('time in h')
        ylabel('conf. prob. and dist.')
    end
    
end

