function [] = show_progress(i, i_max, num_steps)

    if nargin < 3
        num_steps = 10;
    end

    if mod(i,round(i_max/num_steps))==0
        fprintf('%d%% done, loop index = %d\n',int64(round(100*i/i_max)),int64(i));
    end
end

