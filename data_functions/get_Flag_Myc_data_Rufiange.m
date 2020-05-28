function [t_data, N1_N2_F_ratio_mean, error_t_data, N1_N2_F_ratio_sd, log2_N1_N2_F_ratio_mean, log2_N1_N2_F_ratio_sd] = get_Flag_Myc_data_Rufiange()

    t_data = 1.75;  % in h
    error_t_data = 0.25*t_data;  % estimated error

    % Data of the two replicates (sent by email, mean values in Fig. 3F of Rufiange et al.)
    Flag_N1 = [9.26	13.8];
    Flag_N2 = [15.68 28.55];
    F_ratios = Flag_N1./Flag_N2;
    
    N1_N2_F_ratio_mean = mean(F_ratios);
    N1_N2_F_ratio_sd = std(F_ratios);  % sd of ratios

    log2_N1_N2_F_ratio_mean = mean(log2(F_ratios));
    log2_N1_N2_F_ratio_sd = std(log2(F_ratios));
    
end
