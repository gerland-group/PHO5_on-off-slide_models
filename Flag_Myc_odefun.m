function [result] = Flag_Myc_odefun(t, y, M_without_probs, Flag_assembly_prob_fun, myc_indeces, flag_indeces)
    M = M_without_probs;
    p = Flag_assembly_prob_fun(t);
    M(flag_indeces) = M(flag_indeces) * p;
    M(myc_indeces) = M(myc_indeces) * (1-p);
    
    M = M - diag(sum(M));
    
    result = M * y;
end

