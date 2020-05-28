function [shift_rate, W] = calc_chromatin_shift_rate(W_gain)
% the effective chromatin opening/closing rate is the negative eigenvalue closest to zero

    W = W_gain - diag(sum(W_gain));
    evs = real(eig(W));
    if abs(max(evs)) > 1e-10
        error("Highest eigenvalue not close to zero!")
    end
    shift_rate = -max(setdiff(evs, max(evs)));
    if shift_rate <= 0
        error("Shift rate negative!")
    end
end

