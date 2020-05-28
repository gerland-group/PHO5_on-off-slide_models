function [ occs ] = calc_occs_from_configuration_probs( configuration_probs )

    occs = [sum(configuration_probs([1 3 4 5],:),1); sum(configuration_probs([1 2 4 6],:),1); sum(configuration_probs([1 2 3 7],:),1)];

end
