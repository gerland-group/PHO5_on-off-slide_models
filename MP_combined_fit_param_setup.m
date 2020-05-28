function [ C, C_text, C_text_mat, params_con, params_reg ] = MP_combined_fit_param_setup( preset_params, params_never_reg, max_complexity, max_num_params_con, max_num_params_reg, n_data )
% Prepares all the parameter (i.e. process) combinatinations of constitutive and regulated parameters and saves them in params_con and params_reg.
% Also sorts out redundant combinations.
% C_text, C_text_mat contain the process names, C the affected reactions.

% preset_params are part of every combination. They can be constitutive or regulated.
% params_never_reg are not part of every combination, but if they are, they have to be constitutive.
% max_complexity is the number of free parameters (total parameters number - 1).
% n_data decides the number of data sets to fit by changing the regulated parameters only.

    if ~isempty(max_complexity) && max_complexity < n_data + 1
        error('max_complexity has to be at least n_data + 1 (i.e. 2 constitutive and 1 regulated parameter -> 2 + n_data values, one of which is set to 1) or []')
    end
    
    if ~isempty(max_num_params_reg) && max_num_params_reg < 1
        error('max_num_reg_params has to be >= 1 or = []')
    end
    
    if ~isempty(max_num_params_con) && max_num_params_con < 2
        error('max_num_con_params has to be >= 2 or = []')
    end
    
    if isempty(max_complexity) && ( isempty(max_num_params_reg) || isempty(max_num_params_con) )
        error('If max_complexity is empty, max_num_params_reg and max_num_params_con have to be well-defined.')
    end
    
    if ~isempty(max_complexity) && ~isempty(max_num_params_reg) && ~isempty(max_num_params_con)
        warning('max_complexity, max_num_con_params and max_num_reg_params are given. Check if that is really what you want!')
    end
    
    if max_complexity < length(preset_params)
        error('max_complexity has to be at least length(preset_params) (in this case one of the preset_params will become regulated, no additional params will be added)')
    elseif max_complexity == length(preset_params)
        warning('max_complexity == length(preset_params), no additional parameters will be added to the preset parameters')
    end
       
    if max_complexity >= 9
        warning('max_complexity >= 9 for n_data=2 leads to redundant and possible senseless combinations that are not filtered out. If n_data>2, check if it helps. Redundant parameter checks need to be adjusted!')
    end

    % Setup site specific parameters
    A = cell(1,3);
    D = cell(1,3);
    S = cell(1,4);
    A{1} = [1,2; 3,7; 4,6; 5,8;];  % Pos 1 assembly
    A{2} = [1,3; 2,7; 4,5; 6,8;];  % Pos 2 assembly
    A{3} = [1,4; 2,6; 3,5; 7,8;];  % Pos 3 assembly
    D{1} = [2,1; 7,3; 6,4; 8,5;];  % Pos 1 disassembly
    D{2} = [3,1; 7,2; 5,4; 8,6;];  % Pos 2 disassembly
    D{3} = [4,1; 6,2; 5,3; 8,7;];  % Pos 3 disassembly
    S{1} = [3,2; 5,6;];  % Middle -> top sliding
    S{2} = [3,4; 7,6;];  % Middle -> bottom sliding
    S{3} = [2,3; 6,5;];  % Top -> middle sliding
    S{4} = [4,3; 6,7;];  % Bottom -> middle sliding
    A_site_text = {'A1  ', 'A2  ', 'A3  '};
    D_site_text = {'D1  ', 'D2  ', 'D3  '};
    S_site_text = {'S21 ', 'S23 ', 'S12 ', 'S32 '};
    
    % Global assembly, disassembly and sliding
    C{1} = [A{1}; A{2}; A{3};];  C_text = {'A   '}; A_global_param_num = length(C);
    
    % Params in C need to be ordered from more global to more specific:
    % i. e. from the 4th param on, any param has to be a subset of the earlier param or has to have no intersection with the earlier param (for all earlier params)
    % Only then, more specific params overwrite less specific ones!

    A_site_params = [];
    D_site_params = [];
    AD_site_params = [];
    S_combined_params = [];
    
    % Assembly
    for i=1:length(A)
        C{end+1} = A{i};
        A_site_params = [A_site_params length(C)];
        AD_site_params = [AD_site_params length(C)];
        C_text{end+1} = A_site_text{i};
        for j=1:size(A{i},1)
            C{end+1} = A{i}(j,:);
            C_text{end+1} = sprintf('A%d-%d', A{i}(j,2), A{i}(j,1));
        end
        
    end
    
    C{end+1} = [D{1}; D{2}; D{3};];  C_text{end+1} = 'D   '; D_global_param_num = length(C);
    
    % Disassembly
    for i=1:length(D)
        C{end+1} = D{i};
        D_site_params = [D_site_params length(C)];
        AD_site_params = [AD_site_params length(C)];
        C_text{end+1} = D_site_text{i};
        for j=1:size(D{i},1)
            C{end+1} = D{i}(j,:);
            C_text{end+1} = sprintf('D%d-%d', D{i}(j,2), D{i}(j,1));
        end
    end
    
    C{end+1} = [S{1}; S{2}; S{3}; S{4};];  C_text{end+1} = 'S   '; S_global_param_num = length(C);
    C{end+1} = [S{1}; S{2}]; C_text{end+1} = 'S2* '; S_from_2_param_num = length(C);

    % Sliding
    for i=1:length(S)
        C{end+1} = S{i};
        S_combined_params = [S_combined_params length(C)];
        C_text{end+1} = S_site_text{i};
        for j=1:size(S{i},1)
            C{end+1} = S{i}(j,:);
            C_text{end+1} = sprintf('S%d-%d', S{i}(j,2), S{i}(j,1));
        end
    end
    
    % C_text_mat
    for i=1:length(C_text)
        C_text_mat(i,:) = C_text{i};
    end
    C_text_mat = ['-   '; C_text_mat];
    
    if all(preset_params-A_global_param_num) || all(preset_params-D_global_param_num)
        warning('Global assembly and disassembly are not in each combination. Check for unique steady state still needs to be implemented! Redundant parameter checks need to be adjusted!')
    end
        
    % Create parameter combinations
    possible_add_params = setdiff(1:length(C),preset_params);
    max_add_params = min([max_complexity+(2-n_data)-length(preset_params) max_num_params_con+max_num_params_reg-length(preset_params)]);  % total values = c_max+1; minus n_data-1 for one regulated parameter
    add_params = zeros(1, max_add_params);  % first row contains no additional parameters
    for n=1:max_add_params
        add_params_combs = nchoosek(possible_add_params,n);
        
        % Check for redundant combinations:
        keep_param_combs = [];
        for comb=1:size(add_params_combs,1)
            
            if n==max_add_params  % show progress for the last case (should be by far the largest)
                show_progress(comb, size(add_params_combs,1));
            end
            
            params_temp = sort([preset_params add_params_combs(comb,:)]);  % sorting is important here
            delete_comb = 0;
            
            % Delete combinations where single sliding params overwrite all or all but one reactions of a combined sliding param
            if ~delete_comb
                S_combined_in_params = myismember(params_temp,S_combined_params);
                for s=find(S_combined_in_params)
                    if s<length(params_temp) && params_temp(s+1) <= params_temp(s) + 2  % the two single sliding params come right after the combined sliding params
                        delete_comb = 1;
                        break
                    end
                end
            end
            
            % Delete combinations where single assembly/disassembly params overwrite all or all but one reactions of a combined assembly/disassembly param
            if ~delete_comb
                AD_site_in_params = myismember(params_temp,AD_site_params);
                for s=find(AD_site_in_params)
                    if s<length(params_temp)-2 && params_temp(s+3) <= params_temp(s) + 4  % the four single AD params come right afther the combined AD params
                        delete_comb = 1;
                        break
                    end
                end
            end
            
            % Delete combinations with A/D glo and all site A/D params 
            % (keep the ones with A/D glo and all but one combined A/D param as representatives 
            % of all site A/D params without A/D glo when running with preset_params [1 2])
            if ~delete_comb && ( sum(myismember([A_global_param_num A_site_params],params_temp))==4 || sum(myismember([D_global_param_num D_site_params],params_temp))==4 )
                delete_comb = 2;
            end
            
            % If S glo is preset param, delete combination with S glo and ( all combined S params or all but one combined S param and two single reaction sliding params )
            if ~delete_comb &&  myismember(S_global_param_num, preset_params) && myismember(S_global_param_num, params_temp) && ( sum(myismember(S_combined_params,params_temp)) == length(S_combined_params) || sum(myismember(S_combined_params,params_temp)) == length(S_combined_params)-1 && sum(myismember([S_combined_params(1)+[1 2] S_combined_params(2)+[1 2] S_combined_params(3)+[1 2] S_combined_params(4)+[1 2]],params_temp))>=2 )
                delete_comb = 2;
            end
            
            % If S glo is not preset param, delete combination with S glo and all or all but one combined S param or all but two combined S params and more than two reactions specific sliding params
            if ~delete_comb && ~myismember(S_global_param_num, preset_params) && myismember(S_global_param_num, params_temp) && ( sum(myismember(S_combined_params,params_temp)) >= length(S_combined_params)-1 || sum(myismember(S_combined_params,params_temp)) == length(S_combined_params)-2 && sum(myismember([S_combined_params(1)+[1 2] S_combined_params(2)+[1 2] S_combined_params(3)+[1 2] S_combined_params(4)+[1 2]],params_temp))>2  )
                delete_comb = 2;
            end
            
            % Delete combinations of similar to Aglo A1 A2 and 4 out of A3
            if ~delete_comb && myismember(A_global_param_num,params_temp) && sum(myismember(A_site_params,params_temp)) == 2
                if sum(myismember(A_site_params(~myismember(A_site_params,params_temp))+(1:4),params_temp)) == 4
                    delete_comb = 3;
                end
            end
            if ~delete_comb && myismember(D_global_param_num,params_temp) && sum(myismember(D_site_params,params_temp)) == 2
                if sum(myismember(D_site_params(~myismember(D_site_params,params_temp))+(1:4),params_temp)) == 4
                    delete_comb = 3;
                end
            end
            
            % Delete combination with "S from 2" because it could be replaced with S #2 > #3 or S #2 > #1 or a single reaction sliding param
            if ~delete_comb && myismember(S_from_2_param_num, params_temp) && ( sum(myismember(S_combined_params([1 2]), params_temp))>0 || sum(myismember(S_combined_params(1)+[1 2], params_temp))>1 || sum(myismember(S_combined_params(2)+[1 2], params_temp))>1 )
                delete_comb = 4;
            end
            
            % If S glo is not preset param, delete combination with S glo and S from 2 and more than ( zero additional combined S param or 2 single reaction sliding params )
            if ~delete_comb && ~myismember(S_global_param_num, preset_params) && myismember(S_global_param_num, params_temp) && myismember(S_from_2_param_num, params_temp) && ( sum(myismember(S_combined_params([3 4]),params_temp))>0 || sum(myismember([S_combined_params(3)+[1 2] S_combined_params(4)+[1 2]],params_temp))>2 )
                delete_comb = 5;
            end
            
            % If S glo is preset param, delete combination with S glo and S from 2 and with ( two additional combined S param or 4 single reaction sliding params or one additional combined S param and 2 single reaction sliding params )
            if ~delete_comb &&  myismember(S_global_param_num, preset_params) && myismember(S_global_param_num, params_temp) && myismember(S_from_2_param_num, params_temp) && ( sum(myismember(S_combined_params([3 4]),params_temp))>=2 || sum(myismember([S_combined_params(3)+[1 2] S_combined_params(4)+[1 2]],params_temp))>=4 || myismember(S_combined_params(3),params_temp) && sum(myismember(S_combined_params(4)+[1 2],params_temp))>=2 || myismember(S_combined_params(4),params_temp) && sum(myismember(S_combined_params(3)+[1 2],params_temp))>=2 )
                delete_comb = 5;
            end
                       
            if ~delete_comb
                % calc_value_indeces_in_W can also be used to sort out models with missing params or params that occur not often enough
                % missing params should throw an error, because these models should be sorted out before
                % some models with params that occur not often enough are explicitly allowed above
                [ value_indeces_in_W, error_param_missing, error_param_not_often_enough ] = calc_value_indeces_in_W( params_temp, C, 8, 1);
                if ~isempty(error_param_missing)
                    params_temp
                    C_text(params_temp)
                    value_indeces_in_W
                    error_param_missing
                    warning('This parameter combination results in a missing param:')
                elseif ~isempty(error_param_not_often_enough) && ~all(myismember(error_param_not_often_enough, preset_params))
                    params_temp
                    C_text(params_temp)
                    value_indeces_in_W
                    error_param_not_often_enough
                    warning('This parameter combination contains a non-preset param which can be replaced by a single reaction param:')
                end
            end
            
            % If all checks are ok
            if ~delete_comb
                keep_param_combs = [keep_param_combs comb];
            end
        end
        add_params_combs = add_params_combs(keep_param_combs,:);
        add_params = [add_params; add_params_combs zeros(size(add_params_combs,1),max_add_params-n)];
    end
    
    % Determine constitutive and regulated paramaters
    params_con = zeros(1000000,min([max_complexity+1-n_data, max_num_params_con]));  % from a case with one regulated parameter
    params_reg = zeros(1000000,min([floor(max_complexity/n_data), max_num_params_reg]));  % each model has at least one constitive parameter
    
    block_start = 1;
    block_end = 0;
    
    % Add additional params
    for i=1:size(add_params,1)
        
        show_progress(i, size(add_params,1));
        
        add_params_temp = add_params(i,add_params(i,:)>0);
        if ~isempty(params_never_reg)
            params_possibly_reg = setdiff([preset_params add_params_temp],params_never_reg);  % could be optimized to be faster
        else
            params_possibly_reg = [preset_params add_params_temp];
        end
        L = length(add_params_temp) + length(preset_params);  % total length of specific param combination (not counting multiple regulated values)
        
        for nip=1:min([L-1, floor((max_complexity+1-L)/(n_data-1)), length(params_possibly_reg), max_num_params_reg])  % nip: number of regulated params
            params_reg_combs = nchoosek(params_possibly_reg,nip);  % all combinations of nip regulated params
            
            % find left over constitutive parameters (this method is fast than using setdiff on each combination)
            params_mat = (ones(size(params_reg_combs,1),1)*[preset_params add_params_temp])';
            for ip = 1:nip    
                ind_params_mat = (params_reg_combs(:,ip)*ones(1,size(params_mat,1)))';
                params_mat = reshape(params_mat(params_mat-ind_params_mat~=0),[L-ip,size(params_reg_combs,1)]);
            end
            params_con_combs = params_mat';
            
            block_end = block_start + size(params_con_combs,1)-1;
            if block_end > size(params_con,1)  % double the memory
                if 2 * size(params_con,1) > 1000000000
                    error('Tried to allocate memory for more than 1,000,000,000 combinations')
                end
                if 2 * size(params_con,1) > 100000000
                    warning('Allocating memory for more than 100,000,000 combinations')
                end
                params_con = [params_con; zeros(size(params_con))];
                params_reg = [params_reg; zeros(size(params_reg))];
            end
            params_con(block_start:block_end,:) = [params_con_combs zeros(size(params_con_combs,1),size(params_con,2)-size(params_con_combs,2))];
            params_reg(block_start:block_end,:) = [params_reg_combs zeros(size(params_reg_combs,1),size(params_reg,2)-nip)];
            block_start = block_end + 1;

        end
    end
    
    params_con = params_con(1:block_end,:);
    params_reg = params_reg(1:block_end,:);

    fprintf('Found %d combinations of constitutive and regulated parameters.\n',size(params_con,1))
        
end
