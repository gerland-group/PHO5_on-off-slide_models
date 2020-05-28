function indeces = calc_models_with_params_from_param_set(params, param_set)
    indeces = [];
    for i=1:size(params,1)
        if all(myismember(params(i,:),[param_set 0]))
            %params(i,:)
            indeces = [indeces i];
        end
    end
end

