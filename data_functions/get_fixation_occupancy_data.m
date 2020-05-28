function [fix_occ_N2_N3, sd_fix_occ_N2_N3] = get_fixation_occupancy_data(occ_N2_N3, mutant_id)

    % Philipp Korber's N-3 sticky N-3 mutant  restriction enyzme accessibility data (-Pi)
    % sticky N-3 mutants 1 and 2 fold changes and their standard deviation
  
    % [N-2 N-3] (induced)

    p1_fc = [0.6796875 0.5213675]';
    p2_fc = [0.7421875 0.6239316]';
    
    % assuming independent biological replicates
    sd_p1_fc = [0.03638965 0.23043547]';
    sd_p2_fc = [0.07905954 0.13487148]';    

    if mutant_id == 1
        fix_occ_N2_N3 = 1 - p1_fc .* (1-occ_N2_N3);
        sd_fix_occ_N2_N3 = sd_p1_fc .* (1-occ_N2_N3);
    elseif mutant_id == 2
        fix_occ_N2_N3 = 1 - p2_fc .* (1-occ_N2_N3);
        sd_fix_occ_N2_N3 = sd_p2_fc .* (1-occ_N2_N3);
    else
        error("mutant_id not valid.")
    end

end