function [ value_indeces_in_W, error_param_missing, error_param_not_often_enough ] = calc_value_indeces_in_W( params_temp, C, size_of_W, check_mode )
% calculates W matrix without rate values, but parameter indeces in params_temp (not the parameter numbers from setup_params) instead, e.g.:
% 0 1 1 1
% 2 0 3 3
% 3 4 0 4
% 5 5 5 0
% this is useful to determine if a parameter value is overwritten completely or can be replaced by another parameter

        error_param_missing = NaN;
        error_param_not_often_enough = NaN;
        [params_temp_sorted, I] = sort(params_temp);
        C_model = C(params_temp_sorted);
        value_indeces_in_W = zeros(size_of_W);
        for k = 1:length(C_model)
            for j = 1:size(C_model{k},1)
                % I(k) is the index with respect to the unsorted params, positioned where the corresponding param value will have to be:
                value_indeces_in_W(C_model{k}(j,1),C_model{k}(j,2)) = I(k);
            end
        end
        
        if check_mode
            error_param_missing = [];
            for i=1:length(params_temp)
                found_i = find(value_indeces_in_W(:) - i == 0,1);
                if isempty(found_i)
                    error_param_missing = [error_param_missing params_temp(i)];
                end
            end
            
            % Check if params occur often enough to give a new model:
            % i.e. check if there is an equivalent model with a configuration specific parameter because the less specific parameter in this model is overwritten too often
            % For site dependent parameters this check is enough. For global parameters (1,2,3) this check does not find equivalent models using site specific parameters.
            error_param_not_often_enough = [];
            for i=1:length(params_temp)
                found_i = find(value_indeces_in_W(:) - i == 0);
                if length(found_i) < (size(C{params_temp(i)},1)>1)+1
                    error_param_not_often_enough = [error_param_not_often_enough params_temp(i)];
                end
            end
            
            % Check if global params can be substituted by one site specific parameter:
            % Idea: do all substitutions and check if the global param is still there.
        end
        
end
