function [N_state_num] = FM_state_to_N_state(FM_state_num)
% reduces the extendend configurations with Flag and Myc nucleosomes back to the normal nucleosome configurations

    if FM_state_num >= 1 && FM_state_num <= 27
        if FM_state_num <= 8
            N_state_num = 1;
        elseif FM_state_num <= 12
            N_state_num = 4;
        elseif FM_state_num <= 16
            N_state_num = 3;
        elseif FM_state_num <= 20
            N_state_num = 2;
        elseif FM_state_num <= 22
            N_state_num = 7;
        elseif FM_state_num <= 24
            N_state_num = 6;
        elseif FM_state_num <= 26
            N_state_num = 5;
        elseif FM_state_num == 27
            N_state_num = 8;
        end
    else
        error("FM_state_num needs to be 1,2,...,27")
    end
end
