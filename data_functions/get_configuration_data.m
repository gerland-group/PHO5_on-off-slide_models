function [ configuration_count_data, configuration_prob_data, occupancy_prob_data, LHs_data, data_text ] = get_configuration_data( strain_indeces )

    % Data from Brown, Mao, Falkovskaia, Jurica and Boeger (2013)
    data_all = [141 10 29 15 6 1 7 1; ...        % Table S1
                61 8 57 16 28 2 19 12; ...
                14 8 39 11 38 16 37 50; ...
                125 25 28 9 7 2 9 5; ...
                15 4 36 13 43 5 37 50; ...
                109 20 41 12 13 5 7 2]';
    strain_infos = {'rep: pho4D pho80D TATAm', 'half-act: pho4[85-99] pho80D TATAm', 'act: pho80D TATAm', 'rep: wt', 'act: pho80D', 'rep: pho2D'};
    
    configuration_count_data = data_all(:,strain_indeces);
    
    for i=1:length(strain_indeces)
        data_text{i} = strain_infos{strain_indeces(i)};
    end

    configuration_prob_data = configuration_count_data ./ (ones(8,1)*sum(configuration_count_data));
    for i = 1:size(configuration_count_data,2)
        LHs_data(i) = multinomial(configuration_count_data(:,i))*prod(configuration_prob_data(:,i).^configuration_count_data(:,i));
    end

    occupancy_prob_data = calc_occs_from_configuration_probs(configuration_prob_data);

end
