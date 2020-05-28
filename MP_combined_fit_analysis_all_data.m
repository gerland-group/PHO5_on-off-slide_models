% This script runs stages 2 and 3 of the fit procedure and analyses the results

%% Add paths to analysis functions (needed for cluster data analysis)

set(0,'defaultfigureposition',[1700,200,800,1200]')
set(0,'defaultfigurecolor','white')
set(0,'defaultaxesfontsize',10)
warning('WARNING TEST OK, DEFAULTS SET')

addpath external_functions/DERIVESTsuite
addpath external_functions/xticklabel_rotate
addpath external_functions/format_ticks
addpath helper_functions
addpath plot_functions
addpath data_functions

load("temp_results.mat")

if ~exist('img/models', 'dir')
    mkdir('img/models')
end

%% Sort results

n_data = size(configuration_count_data,2);
[LHs,I_best_LH] = max(LHs_mat,[],2);
AICs = -2* LHs*log(10) + (sum(params_con_setup>0,2) + n_data*sum(params_reg_setup>0,2) - 1) * 2;
BICs = -2* LHs*log(10) + (sum(params_con_setup>0,2) + n_data*sum(params_reg_setup>0,2) - 1) * log(mean(sum(configuration_count_data)));  % mean(sum(configuration_count_data)) is the (approx.) number of observations of one gene molecule in both, activated and repressed, conditions
% The BIC calculation could be more rigorous by randomly throwing away excess datapoints such that both conditions have the same number of analysed molecules and then do the likelihood fits.
AIC_lower_bounds = -2*log(prod(LHs_data)) + (1:20)' * 2;
BIC_lower_bounds = -2*log(prod(LHs_data)) + (1:20)' * log(mean(sum(configuration_count_data)));  % given complexity (i.e. number of real free parameters) this is how low the BIC can get

values_con = NaN*ones(size(values_con_mat(:,:,1)));
values_reg = NaN*ones(size(values_reg_mat(:,:,:,1)));
for i=1:size(LHs,1)
    values_con(i,:) = values_con_mat(i,:,I_best_LH(i));
    values_reg(i,:,:) = values_reg_mat(i,:,:,I_best_LH(i));
end

switch 1
    case 1  % sort with respect to LHs
        [LHs, I] = sort(LHs,'descend');
        BICs = BICs(I);
        AICs = AICs(I);
    case 2  % sort with respect to BICs
        [BICs, I] = sort(BICs,'ascend');
        LHs = LHs(I);
        AICs = AICs(I);
    case 3  % sort with respect to AICs
        [AICs, I] = sort(AICs,'ascend');
        LHs = LHs(I);
        BICs = BICs(I);
end

params_con = params_con_setup(I,:);
params_reg = params_reg_setup(I,:);
complexity = sum(params_con>0,2) + n_data*sum(params_reg>0,2) - 1;
sliding = any(params_con>32,2) | any(params_reg>32,2);
values_con = values_con(I,:);
values_reg = values_reg(I,:,:);
values_con_mat = values_con_mat(I,:,:);
values_reg_mat = values_reg_mat(I,:,:,:);
LHs_mat = LHs_mat(I,:);
exitflags_sorted = exitflags(I,:);

%% Analyse uniqueness

LH_tol = 0.05;  % this is an absolute tolerance
value_tol = 0.05;  % relative tolerance with respect to the largest value in the column (ToDo: Improve choice of group representative)
boundary_tol = 1;  % 0, no extra treatment of boundary, any position value is considered extra tolerance near the lower boundary (values below (1+boundary_tol)*min(rate_bounds) are set to (1+boundary_tol)*min(rate_bounds) for comparison)

nonunique_LHs_ind = zeros(size(LHs_mat));
nonunique_best_values_ind = zeros(size(values_con_mat,1),size(values_con_mat,3));
if size(LHs_mat,2)>1
    if max_complexity <= 6  % otherwise the loop blows up memory for complexity >= 7
        gcp;
        number_of_nodes = getParforArg();
    else
        number_of_nodes = 0;  % starts a normal for loop (which for some reason runs backwards)
    end
    %for i=1:size(LHs,1)  % for debugging
    parfor(i=1:size(LHs,1), number_of_nodes)
        show_progress(i,size(LHs,1))
        % Are there any nonunique LHs? (--> optimization did not always find the best LH value)
        [~, ~, temp_same_LHs] = uniquetol(-LHs_mat(i,:), LH_tol, 'DataScale', 1);
        nonunique_LHs_ind(i,:) = temp_same_LHs;  % same LH results get the same number, starting with 1 for the best LH value
        
        % Are there any nonunique values leading to the same best LH? (--> optimization did find more than one set of best input values)
        temp_values_reg = permute(values_reg_mat(i,:,:,temp_same_LHs==1),[4 3 2 1]);  % indeces are now: run_number, dataset, parameter, model_number (which has only one entry, so this dimension is omitted automatically)
        temp_values_reg = temp_values_reg(:,:,temp_values_reg(1,1,:)>0);  % discard NaN entries
        temp_values_reg = reshape(temp_values_reg,sum(temp_same_LHs==1),[]);  % move entries for different parameters next to each other
        temp_values_con = permute(values_con_mat(i,:,temp_same_LHs==1),[3 2 1]);   % indeces are now: run_number, parameter, model_number (which has only one entry, so this dimension is omitted automatically)
        temp_values_con = temp_values_con(:,temp_values_con(1,:)>0);
        temp_values = [temp_values_con temp_values_reg];  % indeces are now: run_number with same best LH, all parameter values (constitutive and regulated for each dataset)

        low_boundary_values = temp_values < (1+boundary_tol)*min(rate_bounds);
        temp_values(low_boundary_values) = (1+boundary_tol)*min(rate_bounds);

        [~, ~, temp_same_values] = uniquetol(temp_values, value_tol, 'ByRows', true);  % The groups around two rows that are not within tolerance are not unique. Representative rows of two different groups can be within tolerance.
        nonunique_best_values_ind_temp = zeros(1,size(values_con_mat,3));
        nonunique_best_values_ind_temp(temp_same_LHs==1) = temp_same_values';
        nonunique_best_values_ind(i,:) = nonunique_best_values_ind_temp;  % contains entries for each model and run: entry=0 -> run did not reach best LH value, entry>=1 -> entries are uniqueness group the run was asigned.
        % E.g. [0 1 0 2 2]: run 1 and 3 did not reach the best LH value, all others did and run 2 is in its own group and run 4 and 5 in another group
    end
end

exitflags_best_LH = exitflags_sorted;
exitflags_best_LH(nonunique_LHs_ind~=1) = NaN;

if size(LHs_mat,2)>1
    models_with_nonunique_LHs = find(max(nonunique_LHs_ind,[],2)>1);
    models_with_nonunique_best_values = find(max(nonunique_best_values_ind,[],2)>1);
    fraction_of_nonunique_LHs = length(models_with_nonunique_LHs) / size(LHs_mat,1)
    fraction_of_nonunique_best_values = length(models_with_nonunique_best_values) / size(LHs_mat,1)
    
    for i=1:size(LHs_mat,2)
        num_models_with_best_LH_i_times(i) = length(find(sum(nonunique_best_values_ind>0,2)==i));
    end
    num_models_with_best_LH_i_times;
end

%% Check models for reversible parameters and compare with nonuniqueness

if 1
    params_temp = zeros(1,size(params_con,2)+size(params_reg,2));
    error_param_missing = zeros(size(params_con,1),max_complexity+1);
    error_param_not_often_enough = zeros(size(params_con,1),max_complexity+1);
    models_all_reverse_params_temp = zeros(size(params_con,1),1);
    models_some_reverse_params_temp = zeros(size(params_con,1),1);

    parfor i = 1:size(params_con,1)
        show_progress(i,size(params_con,1))

        params_temp = [params_con(i,:) params_reg(i,:)];
        params_temp = params_temp(1,params_temp>0);
        [Wi, error_param_missing_temp, error_param_not_often_enough_temp] = calc_value_indeces_in_W( sort(params_temp), C, 8, 1 );
        error_param_missing(i,:) = [error_param_missing_temp zeros(1,max_complexity+1-length(error_param_missing_temp))];
        error_param_not_often_enough(i,:) = [error_param_not_often_enough_temp zeros(1,max_complexity+1-length(error_param_not_often_enough_temp))];

        is_revertible = 1;
        for j = 1:max(max(Wi))
            Wi_transpose = Wi';
            temp = Wi_transpose(Wi==j);
            if all(temp==temp(1)) && temp(1)>0
            else
                is_revertible = 0;
            end
        end
        if is_revertible
            models_all_reverse_params_temp(i) = 1;  % all params are revertible in this model
        end

        Wi2 = zeros(8);
        Wi2(Wi>2) = Wi(Wi>2);
        combs = nchoosek(3:max(max(Wi)),2);
        if size(combs,2)>1
            for j=1:size(combs,1)
                Wi3 = zeros(8);
                Wi3(Wi2==combs(j,1)) = 1;
                Wi3(Wi2==combs(j,2)) = -1;
                if ~any(any(Wi3+Wi3'))
                    models_some_reverse_params_temp(i) = 1;  % some pair of params (excluding global assembly and dissassembly) is a reverse pair
                    break
                end
            end
        end
    end

    models_param_missing = find(any(error_param_missing>0,2));
    models_param_not_often_enough = find(any(error_param_not_often_enough>0,2));
    num_models_param_missing = length(models_param_missing)
    num_models_param_not_often_enough = length(models_param_not_often_enough)

    models_all_reverse_params = find(models_all_reverse_params_temp);
    models_some_reverse_params = find(models_some_reverse_params_temp);

    if size(LHs_mat,2)>1
        [length(models_with_nonunique_best_values) length(models_all_reverse_params) length(models_some_reverse_params)]

        models_all_reverse_and_unique = setdiff(models_all_reverse_params, models_with_nonunique_best_values);  % not all completely revertible models have nonunique values. why?
        length(models_all_reverse_and_unique);

        models_nonunique_without_reverse_params = setdiff(models_with_nonunique_best_values, models_some_reverse_params);  % only possible revertible pair is A glo and D glo, which can cause nonunique results if both are regulated
        length(models_nonunique_without_reverse_params);

        models_AD_glo_ind = find(sum(params_reg==1,2) .* sum(params_reg==2,2));  % most models with A glo and D glo regulated have unique values
        length(intersect(models_nonunique_without_reverse_params,models_AD_glo_ind));
    end
else
    warning('Skipping reversible parameter check.')
end

%% Remove not needed data from workspace to save memory

if contains(pwd,"threshold_6")
    LH_threshold_factor = 100
elseif contains(pwd,"threshold_7")
    LH_threshold_factor = 100 * 10^(1/3)
elseif contains(pwd,"threshold_8")
    LH_threshold_factor = 100 * 100^(1/3)
else
    LH_threshold_factor = 100
end

max_index = find(LHs > sum(log10(LHs_data)) - n_data*log10(LH_threshold_factor),1,'last')

if 1
    clear params_con_setup params_reg_setup exitflags cluster_blocks error_param_missing error_param_not_often_enough % these have not changed and are still available in collected_results.mat

    %params_con = params_con(1:max_index,:);
    %params_reg = params_reg(1:max_index,:);
    complexity = complexity(1:max_index);
    sliding = sliding(1:max_index);
    values_con = values_con(1:max_index,:);
    values_reg = values_reg(1:max_index,:,:);
    values_con_mat = values_con_mat(1:max_index,:,:);
    values_reg_mat = values_reg_mat(1:max_index,:,:,:);
    LHs_mat = LHs_mat(1:max_index,:);
    exitflags_sorted = exitflags_sorted(1:max_index,:);
    exitflags_best_LH = exitflags_best_LH(1:max_index,:);
    nonunique_best_values_ind = nonunique_best_values_ind(1:max_index,:);
    nonunique_LHs_ind = nonunique_LHs_ind(1:max_index,:);
    %AICs = AICs(1:max_index,:);
    %BICs = BICs(1:max_index,:); 
end

%% Calc results table and observables and plot parameter analysis

clear params_con_text params_reg_text
models = 1:max_index;
for i=1:size(params_con,2)
    params_con_text(:,:,i) = C_text_mat(params_con(models,i)+1,:);
end
for i=1:size(params_reg,2)
    params_reg_text(:,:,i) = C_text_mat(params_reg(models,i)+1,:);
end
values_reg_reordered = NaN*ones(max_index,size(values_reg,2)*n_data);
for i=1:max_index  % reorder regulated values to (r1_act, r1_rep, r2_act, r2_rep,...)
    values_reg_temp = permute(values_reg(i,:,:),[3 2 1]);  % indeces are now dataset, parameter
    values_reg_temp = values_reg_temp(:,values_reg_temp(1,:)>0);
    values_reg_temp = reshape(values_reg_temp,1,[]);
    values_reg_reordered(i,:) = [values_reg_temp NaN*ones(1,size(values_reg_reordered,2)-length(values_reg_temp))];
end

[ LHs_single, ps, W_gains, fluxes, occupancy_probs, net_fluxes ] = MP_combined_fit_calc_obs( configuration_count_data, C, params_con, params_reg, values_con, values_reg, max_index, min_S_rate );

results = table(models', num2str(I(models)), AICs(models), BICs(models), LHs(models,:), LHs_single, [complexity(models) BIC_lower_bounds(complexity(models))], sum(nonunique_best_values_ind(models,:)>0,2), max(nonunique_best_values_ind(models,:),[],2), any(exitflags_best_LH(models,:)==1,2) ,'VariableNames',{'rank','orig_index','AIC','BIC','log10_LH','log10_single_LHs','complexity','LH_occ','alt_values','exitflags_1'});
for i=1:size(params_con,2)
    results = [results table(params_con_text(:,:,i), values_con(models,i),'VariableNames',{sprintf('con_%d',i),sprintf('con_%d_v',i)})];
end
for i=1:size(params_reg,2)
    results = [results table(params_reg_text(:,:,i), values_reg_reordered(:,n_data*(i-1)+1:n_data*i),'VariableNames',{sprintf('reg_%d',i),sprintf('reg_%d_v',i)})];
end


LHs_check = max(abs(sum(LHs_single,2)-LHs(1:max_index)));
if LHs_check > 1e-10
    LHs_check
    error('LHs_check too high!')
end

parula_reverse = parula;
parula_reverse = parula_reverse(64:-1:1, :);
c_map = [1 1 1; parula_reverse];

temp = table2array(results(:,'complexity'));  % low complexity    

if 0  % area to check models without net fluxes (detailed balance): run previous section with LH_threshold_factor = 10^100

    results(temp(:,1) < 5 & sliding == 0,:)  % low complexity and no sliding
    
    % test individual models:
    r = 9337;  % with one config param
    r = 6698;  % with one site param
    myfontsize = 10;
    rotation_mode = 'horizontal';

    draw_network_nice_2(C, C_text, params_con(r,:), params_reg(r,:), myfontsize-1,0, rotation_mode)
    plot_flux_networks(fluxes, NaN, ps, r, I, 1)
    plot_net_flux_networks(fluxes, NaN, ps, r, I, 1)
    net_fluxes{1}(:,:,r)
    
    % detailed balance models:
    net_flux_threshold = 1e-15;
    DB_models = find(max(max(net_fluxes{1},[],1),[],2) < net_flux_threshold & max(max(net_fluxes{2},[],1),[],2) < net_flux_threshold & max(max(net_fluxes{3},[],1),[],2) < net_flux_threshold);
    results(DB_models,:)
  
    % non-eq. models with close A-D-symmetry:
    close_symmetry_buddies = [9 10];
    results(close_symmetry_buddies,:)

    a = 1;
    subplot(1,2,a)
    draw_network_nice_2(C, C_text, params_con(close_symmetry_buddies(a),:), params_reg(close_symmetry_buddies(a),:), myfontsize-1,0, rotation_mode)

    a = 2;
    subplot(1,2,a)
    draw_network_nice_2(C, C_text, params_con(close_symmetry_buddies(a),:), params_reg(close_symmetry_buddies(a),:), myfontsize-1,0, rotation_mode)

    plot_matrix_value_analysis_2('net fluxes', 'log', fluxes, NaN, C, C_text, close_symmetry_buddies, I, close_symmetry_buddies, 'clustering manual', [1 2], NaN, c_map, data_text);
    plot_flux_networks(fluxes, NaN, ps, close_symmetry_buddies, I, [1 2])
    plot_net_flux_networks(fluxes, NaN, ps, close_symmetry_buddies, I, [1 2])
end

LH_offset = sum(log10(LHs_data));

plot_parameter_value_analysis(values_con, values_reg, params_con, params_reg, C_text, ...
    [I LH_offset-LHs], NaN, 1:30, NaN, 'clustering off', NaN, NaN, c_map, "s1_top_30");

plot_parameter_value_analysis(values_con, values_reg, params_con, params_reg, C_text, ...
    [I LH_offset-LHs], NaN, find(temp(:,1)==5), NaN, 'clustering off', NaN, NaN, c_map, "s1_complexity_5");

plot_parameter_value_analysis(values_con, values_reg, params_con, params_reg, C_text, ...
    [I LH_offset-LHs], NaN, find(temp(:,1)==6), NaN, 'clustering off', NaN, NaN, c_map, "s1_complexity_6");

%% Expand nonunique results

if size(LHs_mat,2)>1
    expand_models = 1:min(50000,max_index);
    
    first_columns = table(NaN,"",NaN,NaN,NaN,NaN*ones(1,n_data),[NaN NaN],NaN,NaN,NaN,'VariableNames',{'rank','orig_index','AIC','BIC','log10_LH','log10_single_LHs','complexity','LH_occ','alt_values','exitflags_1'});
    
    results_nonunique_started = 0;
    models_nonunique_LHs_single = [];
    for i=expand_models
        num_values = max(nonunique_best_values_ind(i,:));
        
        if num_values>1
            LHs_single_temp = zeros(num_values,n_data);
            for j=1:num_values
                run_indeces = find(nonunique_best_values_ind(i,:)==j);
                run_index = find(nonunique_best_values_ind(i,:)==j,1);

                new_entry = first_columns;
                new_entry.rank = results(i,:).rank;
                if j==1
                    new_entry.orig_index = results(i,:).orig_index;
                    new_entry.AIC = results(i,:).AIC;
                    new_entry.BIC = results(i,:).BIC;
                    new_entry.complexity = results(i,:).complexity;
                end
                new_entry.LH_occ = length(run_indeces);
                new_entry.alt_values = j;  % jth set of values
                new_entry.exitflags_1 = any(exitflags_sorted(i,run_indeces)==1);
                
                LHs_single_temp(j,:) = MP_combined_fit_calc_obs(configuration_count_data, C, params_con(i,:), params_reg(i,:), values_con_mat(i,:,run_index), values_reg_mat(i,:,:,run_index), 1, min_S_rate);
                
                new_entry.log10_LH = LHs_mat(i,run_index);
                new_entry.log10_single_LHs = LHs_single_temp(j,:);

                values_reg_temp = values_reg_mat(i,:,:,run_index);
                values_reg_temp = permute(values_reg_temp,[3 2 1]);
                values_reg_reordered = reshape(values_reg_temp,1,[]);
                
                for k=1:size(params_con,2)
                    temp_table = table([params_con_text(i,:,k); params_con_text(i,:,k)], [values_con_mat(i,k,run_index); values_con_mat(i,k,run_index)],'VariableNames',{sprintf('con_%d',k),sprintf('con_%d_v',k)});  % if matlab doesn't want to create tables of one line with character arrays...
                    new_entry = [new_entry temp_table(1,:)];  % ...we create one with two lines and just take the first
                end
                for k=1:size(params_reg,2)
                    temp_table = table([params_reg_text(i,:,k); params_reg_text(i,:,k)], [values_reg_reordered(n_data*(k-1)+1:n_data*k); (n_data*(k-1)+1:n_data*k)],'VariableNames',{sprintf('reg_%d',k),sprintf('reg_%d_v',k)});
                    new_entry = [new_entry temp_table(1,:)];
                end
                if results_nonunique_started
                    results_nonunique_values = [results_nonunique_values; new_entry];
                else
                    results_nonunique_values = new_entry;
                    results_nonunique_started = 1;
                end
            end
            [~, ~, temp_same_LHs_single] = uniquetol(-LHs_single_temp, LH_tol, 'ByRows', true, 'DataScale', 1);
            if any(temp_same_LHs_single~=1)
                models_nonunique_LHs_single = [models_nonunique_LHs_single i];
            end
        end
    end
    
    if results_nonunique_started
        results_nonunique_values = [results_nonunique_values table(results_nonunique_values.rank,'VariableNames',{'rank2'})];  % helps with large tables
        results_nonunique_values;

        %results_nonunique_values(myismember(results_nonunique_values.rank,models_nonunique_without_reverse_params(1:10)),:);
        results_nonunique_values(myismember(results_nonunique_values.rank,models_nonunique_LHs_single),:);  % these models show a tradeoff between the two maximum likelihood fits

        results_nonunique_LHs = results(models_with_nonunique_LHs(models_with_nonunique_LHs<max_index),:);  % only shows the values with the best seen LH
    end
else
    warning('Only one realisation in dataset. No statements about uniqueness possible.')
end

%% Exitflags (ToDo: check different exitflags in each model)

% For [x,fval,exitflag,output] = fmincon(___) exitflag = ...  means:
% 1: First-order optimality measure was less than options.OptimalityTolerance, and maximum constraint violation was less than options.ConstraintTolerance.
% 2: Change in x was less than options.StepTolerance and maximum constraint violation was less than options.ConstraintTolerance.
% 0: Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.
% In general, strictly positive exit flags are successful outcomes.

figure
subplot(2,1,1)
histogram(exitflags_sorted(:))
title('exitflags')

num_exitflags_zero = sum(sum(exitflags_sorted == 0))

subplot(2,1,2)
histogram(-LHs_mat(exitflags_sorted==0),0:0.5:80)
xlim([10 80])
xlabel('- log_{10}(LH)')
ylabel('Number of networks')
title('Combined fit, models with one exitflag=0')

%% Check/plot parameter occurences

switch 1
    case 1
        LH_bound = sum(log10(LHs_data)) - n_data*log10(LH_threshold_factor);
        good_s1 = find(LHs > LH_bound);
        fprintf('Found %d models with LH > %f.\n',length(good_s1),LH_bound)
    case 2
        BIC_bound = 104;
        good_s1 = find(BICs < BIC_bound);
        fprintf('Found %d models with BIC < %f.\n',length(good_s1),BIC_bound)
end

max_good_s1 = 100000;
if length(good_s1) > max_good_s1
    warning('Found more than %d good_s1. Reducing to the first %d.',max_good_s1,max_good_s1)
    good_s1 = good_s1(1:max_good_s1);
end

parameter_occurrences_histogram_pdf(params_con,params_reg,C_text,'param_occurrences_all_models',"",1,"among all models")
parameter_occurrences_histogram_pdf(params_con(good_s1,:),params_reg(good_s1,:),C_text,'param_occurrences_good_s1',"",0,"among good stage 1 models")

%% Plot LH and BIC histogram

LH_perfect = sum(log10(LHs_data));
x_max = 12;

plot_LH_histogram(LHs, LH_bound, LH_perfect, x_max, 'All models', '_s1_all_models', 'R_1')

no_sliding_models = find(~any(params_con' > 32) & ~any(params_reg' > 32));
plot_LH_histogram(LHs(no_sliding_models), LH_bound, LH_perfect, x_max, 'Models without sliding', '_s1_no_sliding', 'R_1')
10^(LHs(1) - LHs(no_sliding_models(1)))

config_params = [3 4 5 6 8 9 10 11 13 14 15 16 19 20 21 22 24 25 26 27 29 30 31 32 36 37 39 40 42 43 45 46];
no_config_params = setdiff(0:46, config_params);
no_config_param_models = find(all(ismember(params_con,no_config_params)') & all(ismember(params_reg,no_config_params)'));
plot_LH_histogram(LHs(no_config_param_models), LH_bound, LH_perfect, x_max, 'Models without configuration-specific processes', '_s1_no_config_params', 'R_1')
10^(LHs(1) - LHs(no_config_param_models(1)))

plot_LH_histogram(LHs(intersect(no_sliding_models,no_config_param_models)), LH_bound, LH_perfect, x_max, 'Models without sliding and config.-spec. proc.', '_s1_no_sliding_no_config_params', 'R_1')

lower_complexity_models = find(sum(params_con>0,2) + n_data*sum(params_reg>0,2) - 1 < max(sum(params_con>0,2) + n_data*sum(params_reg>0,2) - 1));
plot_LH_histogram(LHs(lower_complexity_models), LH_bound, LH_perfect, x_max, 'Models with lower complexity', '_s1_lower_complexity', 'R_1')
10^(LHs(1) - LHs(lower_complexity_models(1)))

save(sprintf("stage_1_results_LH_threshold_%d.mat", LH_threshold_factor))

%% Find topological equilibrium models in stage 1 (i.e. models with equal balance no matter what parameter values)

only_global = [1 17 33];
C_text(only_global)

no_sliding_no_config = [1 2 7 12 17 18 23 28];
C_text(no_sliding_no_config)

params = [params_con params_reg];

eq_models = sort([calc_models_with_params_from_param_set(params, only_global) calc_models_with_params_from_param_set(params, no_sliding_no_config)]);

length(eq_models)
min(LH_perfect - LHs(eq_models))  % logarithmic likelihood ratios

%% Chromatin opening/closing rate ratios after stage 1 (only Eigenvalues, fast)

opening_rates = [];
closing_rates = [];
A_reg_models = false(size(good_s1));
D_reg_models = false(size(good_s1));

for i = good_s1'
    params_reg_temp = params_reg(i,:);
    params_reg_temp = params_reg_temp(params_reg_temp > 0);
    if all(params_reg_temp <= 16)
        A_reg_models(i) = true;
    elseif all(params_reg_temp >= 17 & params_reg_temp <= 32)
        D_reg_models(i) = true;
    end
    opening_rates(i) = calc_chromatin_shift_rate(W_gains{1}(:,:,i));  % dynamics into the activated state
    closing_rates(i) = calc_chromatin_shift_rate(W_gains{2}(:,:,i));  % dynamics into the repressed state
end
shift_rate_ratios = closing_rates./opening_rates;

figure

subplot(2,1,1)
hist(shift_rate_ratios(A_reg_models))
title("A-reg.")

subplot(2,1,2)
hist(shift_rate_ratios(D_reg_models))
title("D-reg.")

%% Chromatin opening/closing rate ratios after stage 1 (with distance decay rates, slow)

A_reg_models = false(size(good_s1));
D_reg_models = false(size(good_s1));
opening_rates_ev = zeros(size(good_s1));
closing_rates_ev = zeros(size(good_s1));
opening_rates_dist = zeros(size(good_s1));
closing_rates_dist = zeros(size(good_s1));
tic
for i = good_s1'
    params_reg_temp = params_reg(i,:);
    params_reg_temp = params_reg_temp(params_reg_temp > 0);
    if all(params_reg_temp <= 16)
        A_reg_models(i) = true;
    elseif all(params_reg_temp >= 17 & params_reg_temp <= 32)
        D_reg_models(i) = true;
    end

    data_set = 1;  % starting from repressed state using the dynamics of the activated state
    [opening_rates_ev(i), opening_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains{data_set}(:,:,i), ps(i,:,setdiff([1 2], data_set))', false);

    data_set = 2;  % starting from repressed state using the dynamics of the activated state
    [closing_rates_ev(i), closing_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains{data_set}(:,:,i), ps(i,:,setdiff([1 2], data_set))', false);
end
toc
shift_rate_ratios_ev = closing_rates_ev./opening_rates_ev;
shift_rate_ratios_dist = closing_rates_dist./opening_rates_dist;

figure('Position',[1500,200,400,600],'PaperUnits', 'centimeters','PaperSize', [10 15],'PaperPosition', [0 0 10 15])
defaultaxesfontsize = get(0,'defaultaxesfontsize');
myfontsize = 5;
set(0,'defaultaxesfontsize',myfontsize);

subplot(3,2,1)
hist(shift_rate_ratios_ev(A_reg_models))
title('Eigenvalue rate ratios (A-reg. models)')

subplot(3,2,3)
hist(shift_rate_ratios_dist(A_reg_models))
title('Distance decay rate ratios (A-reg. models)')

subplot(3,2,5)
plot(shift_rate_ratios_ev(A_reg_models), shift_rate_ratios_dist(A_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (A-reg. models)')
ylabel('Distance decay rate ratios (A-reg. models)')

subplot(3,2,2)
hist(shift_rate_ratios_ev(D_reg_models))
title('Eigenvalue rate ratios (D-reg. models)')

subplot(3,2,4)
hist(shift_rate_ratios_dist(D_reg_models))
title('Distance decay rate ratios (D-reg. models)')

subplot(3,2,6)
plot(shift_rate_ratios_ev(D_reg_models), shift_rate_ratios_dist(D_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (D-reg. models)')
ylabel('Distance decay rate ratios (D-reg. models)')

pause(0.5)
print(sprintf('img/shift_rate_ratios_s1.pdf'), '-dpdf')
close gcf
open(sprintf('img/shift_rate_ratios_s1.pdf'))
set(0,'defaultaxesfontsize',defaultaxesfontsize);

%% Stage 2 run (with sticky N-3 data)

prefactor_bounds_log = log([1/5 5]);
num_prefactors = 4;  % 4 or 12

fixed_param  % to check what was used in the first stage

stage = 2;
models = good_s1;
dummy_prefactors = NaN*ones(length(models), 2*num_prefactors);  % needed for parfor to work

gcp; 
tic
[~, Fs_s2, F_terms_s2, log10_LHs_s2, ps_s2, W_gains_s2, squ_sticky_N3_errors_s2, occs_sticky_N3_s2, ...
    values_con_s2, values_reg_s2, prefactors_s2, exitflags_s2, outputs_s2, lambdas_s2] = ...
    MP_combined_fit_all_data(stage, LHs_data, configuration_count_data, C, params_con(models,:), params_reg(models,:,:), fixed_param, rate_bounds, min_S_rate, ...
        values_con(models,:), values_reg(models,:,:), dummy_prefactors, NaN, NaN, NaN, prefactor_bounds_log, num_prefactors);
toc

%% Stage 2 results

F_max = n_data*log10(LH_threshold_factor);  % = -(LH_bound-LH_perfect)
good_s2 = find(Fs_s2 <= F_max);

plot_LH_histogram(-Fs_s2, LH_bound-LH_perfect, 0, x_max, 'All models', '_s2_all_models', 'R_2')

no_sliding_models_s2 = no_sliding_models(no_sliding_models <= max(good_s1));
plot_LH_histogram(-Fs_s2(no_sliding_models_s2), LH_bound-LH_perfect, 0, x_max, 'Models without sliding', '_s2_no_sliding', 'R_2')

no_config_param_models_s2 = no_config_param_models(no_config_param_models <= max(good_s1));
plot_LH_histogram(-Fs_s2(no_config_param_models_s2), LH_bound-LH_perfect, 0, x_max, 'Models without configuration-specific processes', '_s2_no_config_params', 'R_2')

lower_complexity_models_s2 = lower_complexity_models(lower_complexity_models <= max(good_s1));
plot_LH_histogram(-Fs_s2(lower_complexity_models_s2), LH_bound-LH_perfect, 0, x_max, 'Models with lower complexity', '_s2_lower_complexity', 'R_2')

parameter_occurrences_histogram_pdf(params_con(good_s2,:),params_reg(good_s2,:),C_text,'param_occurrences_good_s2',"",0,"among good stage 2 models")

figure
subplot(4,1,1)
histogram(Fs_s2,0:0.1:1.5*F_max)
xlabel("stage 2 minimal objective function value")
ylabel("count")
subplot(4,1,2)
scatter(F_terms_s2(:,1), F_terms_s2(:,2),'.')
xlim([0 1.5*F_max]);
ylim([0 1.5*F_max]);
xlabel("likelihood term")
ylabel("sticky N-3 error term")
subplot(4,1,3)
diff_values_con = abs(log10(values_con_s2)-log10(values_con(good_s1,:)));
diff_values_reg = abs(log10(values_reg_s2)-log10(values_reg(good_s1,:,:)));
%diff_values_con(isnan(diff_values_con)) = 0;
%diff_values_reg(isnan(diff_values_reg)) = 0;
%histogram(max([max(diff_values_reg,[],3) diff_values_con],[],2))
histogram(mean([mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') diff_values_con],2,'omitnan'), 0:0.04:2)
xlim([0 2])
xlabel("mean abs. log10 parameter change from s1")
ylabel("count")
subplot(4,1,4)
diff_values_con = diff_values_con(good_s2,:);
diff_values_reg = diff_values_reg(good_s2,:,:);
histogram(mean([mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') diff_values_con],2,'omitnan'), 0:0.04:2)
xlim([0 2])
xlabel("mean abs. log10 parameter change from s1 of good s2 models")
ylabel("count")

mean(sum(log10_LHs_s2(good_s2, :), 2))
mean(LHs(good_s2))

x = exitflags_s2;
binc = [min(0,min(x)):max(x)];
counts = hist(x,binc);
occurrences = [binc; counts]

table_s2 = table(num2str(I(good_s1)), Fs_s2, F_terms_s2, exitflags_s2, values_con_s2, values_reg_s2(:,:,1), values_reg_s2(:,:,2), values_reg_s2(:,:,3));
table_s2_good = table_s2(good_s2, :);

table_s2 = sortrows(table_s2);
table_s2_good = sortrows(table_s2);

results_s2_good = results(good_s2,1:end);

length(good_s2)

%% Stage 2 optimization infos

for i = find(exitflags_s2<=0)'
    outputs_s2(i)
end

%% Stage 2 log prefactor plots and parameter analysis

stage = 2;
max_log_factor = prefactor_bounds_log(2);

log_prefactors_s2 = NaN * ones(length(good_s1),num_prefactors,2);
log_prefactors_s2(:,:,1) = log(prefactors_s2(:,1:num_prefactors));
log_prefactors_s2(:,:,2) = log(prefactors_s2(:,(num_prefactors+1):(2*num_prefactors)));

if num_prefactors == 4
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 9],'PaperPosition', [-1.5 -0.2 21 9.5])
    prefactor_names = ["a3", "s23", "d3", "s32"];

else
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 27],'PaperPosition', [-1.5 -2.5 21 31])
    prefactor_names = ["a4-1", "a6-2", "a5-3", "a8-7", "s4-3", "s6-7", "d1-4", "d2-6", "d3-5", "d7-8", "s3-4", "s7-6"];
end

i=0;
for p=1:2
    for j=1:num_prefactors
        i=i+1;
        if num_prefactors == 4
            subplot(2,4,i)
        else
            subplot(6,4,i)
        end
        histogram(log_prefactors_s2(good_s2,j,p),-max_log_factor-0.05:0.1:max_log_factor+0.05)
        set(gca, 'FontSize',6)
        xlabel("log_{10} prefactor")
        if j== 1
            ylabel("Count")
        end
        if j == 2
            title({sprintf("                                                                     Sticky mutant %d", p), sprintf("%s", prefactor_names(j))}, 'fontsize', 8)
        else
            title({"", sprintf("%s", prefactor_names(j))}, 'fontsize', 8)
        end
    end
end
pause(0.5)
print("img/sticky_N3_prefactor_histograms_stage_" + stage,'-dpdf')
close gcf
open(sprintf('img/sticky_N3_prefactor_histograms_stage_%d.pdf', stage))

if num_prefactors == 4  % correlation plots only when using 4 prefactors (would be 66 plots for 12 prefactors)
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 22],'PaperPosition', [-1.5 -2 21 25])
    i=0;
    for p=1:2
        for j=1:3
            for k=(j+1):4
                i=i+1;
                subplot(4,3,i)
                scatter(log_prefactors_s2(good_s2,j,p), log_prefactors_s2(good_s2,k,p), 3, 'filled')
                set(gca, 'FontSize',6)
                xlim([-max_log_factor max_log_factor])
                ylim([-max_log_factor max_log_factor])
                if j == 1 && k == 3
                    title({sprintf("Sticky mutant %d",p); ""}, 'fontsize', 8)
                end
                xlabel(prefactor_names(j))
                ylabel(prefactor_names(k))
            end
        end
    end

    pause(0.5)
    print("img/sticky_N3_prefactor_scatterplots_stage_" + stage,'-dpdf')
    close gcf
    open(sprintf('img/sticky_N3_prefactor_scatterplots_stage_%d.pdf', stage))
end

[~, I_clustering_s2] = sort(Fs_s2(good_s2));

plot_parameter_value_analysis(values_con_s2, values_reg_s2, params_con, params_reg, C_text, ...
    [I(good_s1) Fs_s2], NaN, good_s2, NaN, 'clustering manual', I_clustering_s2, NaN, c_map, "s2");

plot_parameter_value_analysis(values_con_s2, values_reg_s2, params_con, params_reg, C_text, ...
    [I(good_s1) Fs_s2], NaN, good_s2, NaN, 'clustering manual', I_clustering_s2(1:min(30,length(I_clustering_s2))), NaN, c_map, "s2_top_30");

save(sprintf('stage_2_results.mat'))

%% Stage 4 run (with sticky N-3 data and with Flag/Myc data, runtime possibly > 1h)

% Stage 4 ("_s4") in program code corresponds to stage 3 in the paper manuscript.

stage = 4;
s4_cluster = 0;
models = good_s2;  % also start directly from stage 2 results
length(models)  % has to be N in MP_stage_3_submit_job.sh
fit_time_scale = 1;
time_scale_bounds = rate_bounds;  % only if equal to rate_bounds, the boundary treatment in sensitivity analysis is correct

if ~s4_cluster  % local calculation
    gcp;
    tic
    [Fs_fmincon_s4, Fs_s4, F_terms_s4, log10_LHs_s4, ps_s4, W_gains_s4, squ_sticky_N3_errors_s4, occs_sticky_N3_s4, ...
        values_con_s4, values_reg_s4, prefactors_s4, exitflags_s4, outputs_s4, lambdas_s4, time_rescaling_s4, num_interations_s4, Flag_probs_s4, ...
        cond_Flag_probs_s4, fmincon_gradients_s4, fmincon_hessians_s4, derivest_gradients_s4, derivest_hessians_s4] = ...
        MP_combined_fit_all_data(stage, LHs_data, configuration_count_data, C, params_con(models,:), params_reg(models,:,:), fixed_param, rate_bounds, min_S_rate, ...
            values_con_s2(models,:), values_reg_s2(models,:,:), prefactors_s2(models,:), fit_time_scale, time_scale_bounds, F_max, prefactor_bounds_log, num_prefactors);
    toc
else
    clear JOB_ID SGE_TASK_ID SGE_TASK_LAST TMPDIR HOSTNAME 
    save 'stage_4/data_for_s4.mat' 'stage' 'LHs_data' 'configuration_count_data' 'C' 'params_con' 'params_reg' 'fixed_param' 'rate_bounds' 'min_S_rate' ...
            'values_con' 'values_reg' 'values_con_s2' 'values_reg_s2' 'prefactors_s2' 'n_data' 'models' 'fit_time_scale' 'time_scale_bounds' 'F_max'
    
    % NOW: run MP_stage_4_submit_job.sh in terminal, wait and collect results
end

%% Stage 4 save data

save(sprintf('stage_4_results.mat'))

%% Stage 4 results

complexity_threshold = 7;

below_complexity_threshold = table2array(results_s2_good(:,'complexity')) <= complexity_threshold;

good_s4 = good_s2(Fs_s4 <= F_max & below_complexity_threshold);
good_s4_in_good_s2 = find(Fs_s4 <= F_max & below_complexity_threshold);

if length(good_s4) > 0
    plot_LH_histogram(-Fs_s4, LH_bound-LH_perfect, 0, x_max, 'All models', '_s4_all_models', 'R_3')

    no_config_param_models_s4 = ismember(no_config_param_models_s2, good_s2);
    plot_LH_histogram(-Fs_s4(no_config_param_models_s4), LH_bound-LH_perfect, 0, x_max, 'Models without configuration-specific processes', '_s4_no_config_params', 'R_3')
    
    parameter_occurrences_histogram_pdf(params_con(good_s4,:),params_reg(good_s4,:),C_text,'param_occurrences_good_s4',"",0,"among good stage 3 models")
end

figure
subplot(4,2,1)
histogram(Fs_s4, 0:0.1:1.5*F_max)
xlabel("F s4 value")
ylabel("count")

subplot(4,2,2)
histogram((F_terms_s4(:,1) - F_terms_s2(good_s2,1))./(F_terms_s4(:,1) + F_terms_s2(good_s2,1))/2, -1:0.025:1)
xlim([-1 1]);
xlabel("rel. diff. 1st term of stage 4 - 2 ")
ylabel("count")

subplot(4,2,3)
diff_values_con = abs(log10(values_con_s4)-log10(values_con(good_s2,:)));
diff_values_reg = abs(log10(values_reg_s4)-log10(values_reg(good_s2,:,:)));
%diff_values_con(isnan(diff_values_con)) = 0;
%diff_values_reg(isnan(diff_values_reg)) = 0;
%histogram(max([max(diff_values_reg,[],3) diff_values_con],[],2), 0:0.04:2)
histogram(mean([mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') mean(diff_values_reg,3,'omitnan') diff_values_con],2,'omitnan'), 0:0.04:2)
xlabel("mean abs. log10 parameter change from stage 1")
ylabel("count")

subplot(4,2,[5:8])
scatter3(F_terms_s4(:,1), F_terms_s4(:,2), F_terms_s4(:,3), ".")
xlim([0 1.5*F_max]);
xlabel("s4 Boeger log10 LH")
ylabel("s4 sticky N-3 log10 LH")
zlabel("s4 Flag/Myc N-1 log10 LH")

length(good_s4)

table_s4 = table(num2str(I(good_s2)), Fs_s4, F_terms_s4, exitflags_s4);
table_s4_good = table_s4(good_s4_in_good_s2,:);

table_s4 = sortrows(table_s4);
table_s4_good = sortrows(table_s4_good);

x = exitflags_s4;
binc = [min(0,min(x)):max(x)];
counts = hist(x,binc);
occurrences = [binc; counts]

length(good_s4)

%% Stage 4 optimization infos

for i = find(exitflags_s4 <= 0)'
    outputs_s4(i)
end

%% Stage 4 log prefactor plots

stage = 4;

log_prefactors_s4 = NaN * ones(length(good_s2),num_prefactors,2);
log_prefactors_s4(:,:,1) = log(prefactors_s4(:,1:num_prefactors));
log_prefactors_s4(:,:,2) = log(prefactors_s4(:,(num_prefactors+1):(2*num_prefactors)));

if num_prefactors == 4
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 9],'PaperPosition', [-1.5 -0.2 21 9.5])
    prefactor_names = ["a3", "s23", "d3", "s32"];

else
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 27],'PaperPosition', [-1.5 -2.5 21 31])
    prefactor_names = ["a4-1", "a6-2", "a5-3", "a8-7", "s4-3", "s6-7", "d1-4", "d2-6", "d3-5", "d7-8", "s3-4", "s7-6"];
end

i=0;
for p=1:2
    for j=1:num_prefactors
        i=i+1;
        if num_prefactors == 4
            subplot(2,4,i)
        else
            subplot(6,4,i)
        end
        histogram(log_prefactors_s4(good_s4_in_good_s2,j,p),-max_log_factor-0.05:0.1:max_log_factor+0.05)
        set(gca, 'FontSize',6)
        xlabel("log_{10} prefactor")
        if j== 1
            ylabel("Count")
        end
        if j == 2
            title({sprintf("                                                                     Sticky mutant %d", p), sprintf("%s", prefactor_names(j))}, 'fontsize', 8)
        else
            title({"", sprintf("%s", prefactor_names(j))}, 'fontsize', 8)
        end
    end
end
pause(0.5)
print("img/sticky_N3_prefactor_histograms_stage_" + stage,'-dpdf')
close gcf
open(sprintf('img/sticky_N3_prefactor_histograms_stage_%d.pdf', stage))

if num_prefactors == 4  % correlation plots only when using 4 prefactors (would be 66 plots for 12 prefactors)
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 22],'PaperPosition', [-1.5 -2 21 25])
    i=0;
    for p=1:2
        for j=1:3
            for k=(j+1):4
                i=i+1;
                subplot(4,3,i)
                scatter(log_prefactors_s4(good_s4_in_good_s2,j,p), log_prefactors_s4(good_s4_in_good_s2,k,p), 3, 'filled')
                set(gca, 'FontSize',6)
                xlim([-max_log_factor max_log_factor])
                ylim([-max_log_factor max_log_factor])
                if j == 1 && k == 3
                    title({sprintf("Sticky mutant %d",p); ""}, 'fontsize', 8)
                end
                xlabel(prefactor_names(j))
                ylabel(prefactor_names(k))
            end
        end
    end

    pause(0.5)
    print("img/sticky_N3_prefactor_scatterplots_stage_" + stage,'-dpdf')
    close gcf
    open(sprintf('img/sticky_N3_prefactor_scatterplots_stage_%d.pdf', stage))
end

%% Stage 4 time rescaling

plot_log2_F_M_ratio_diffs(cond_Flag_probs_s4(:,:,below_complexity_threshold), Fs_s4(below_complexity_threshold), F_max)

figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 9],'PaperPosition', [0 0 18 9])
subplot(1,2,1)
histogram(log10(time_rescaling_s4),log10(time_scale_bounds(1)):0.25:log10(time_scale_bounds(2)))
xlabel('log10 time rescaling factor')
ylabel('Count')
subplot(1,2,2)
histogram(log10(time_rescaling_s4(good_s4_in_good_s2)),log10(time_scale_bounds(1)):0.25:log10(time_scale_bounds(2)))
xlabel('log10 time rescaling factor (good models)')
ylabel('Count')

pause(0.5)
print("img/rescaling_factors",'-dpdf')
close gcf
open(sprintf('img/rescaling_factors.pdf', stage))

%% Stage 4 N-1 over N-2 Flag prob. plots

if fit_time_scale
    %plot_log2_Flag_ratios_subsets(Flag_probs_s4, good_s2(1:size(Flag_probs_s4,3)), no_sliding_models, no_config_param_models, lower_complexity_models, NaN, occupancy_prob_data, Fs_s4);
    %plot_log2_Flag_ratios_subsets(Flag_probs_s4(:,:,good_s4_in_good_s2), good_s4, no_sliding_models, no_config_param_models, lower_complexity_models, NaN, occupancy_prob_data, Fs_s4(good_s4_in_good_s2));
    plot_log2_Flag_ratios_subsets(Flag_probs_s4(:,:,below_complexity_threshold), find(below_complexity_threshold), no_sliding_models, no_config_param_models, lower_complexity_models, NaN, occupancy_prob_data, Fs_s4(below_complexity_threshold), F_max);
else
    plot_H3_exchange_Flag_ratio_subsets(Flag_probs_s4, good_s2(1:size(Flag_probs_s4,3)), no_sliding_models, no_config_param_models, lower_complexity_models, ts, 1, 0, 0, 'ellipse', NaN, occupancy_prob_data);
    plot_H3_exchange_Flag_ratio_subsets(Flag_probs_s4, good_s2(1:size(Flag_probs_s4,3)), no_sliding_models, no_config_param_models, lower_complexity_models, ts, 1, 0, 0, 'ellipse', 0.5, occupancy_prob_data);
end

%% Stage 4 calc results table and observables

clear params_con_text params_reg_text
models = good_s2;
for i=1:size(params_con,2)
    params_con_text(:,:,i) = C_text_mat(params_con(models,i)+1,:);
end
for i=1:size(params_reg,2)
    params_reg_text(:,:,i) = C_text_mat(params_reg(models,i)+1,:);
end
values_reg_reordered = NaN*ones(length(models),size(values_reg_s4,2)*n_data);
for i=1:length(models)  % reorder regulated values to (r1_act, r1_rep, r2_act, r2_rep,...)
    values_reg_temp = permute(values_reg_s4(i,:,:),[3 2 1]);  % indeces are now dataset, parameter
    values_reg_temp = values_reg_temp(:,values_reg_temp(1,:)>0);
    values_reg_temp = reshape(values_reg_temp,1,[]);
    values_reg_reordered(i,:) = [values_reg_temp NaN*ones(1,size(values_reg_reordered,2)-length(values_reg_temp))];
end

[ LHs_single_s4, ps_s4, W_gains_s4, fluxes_s4, occupancy_probs_s4 ] = MP_combined_fit_calc_obs( configuration_count_data, C, params_con(good_s2,:), ...
    params_reg(good_s2,:,:), values_con_s4, values_reg_s4, length(models), min_S_rate );

results_s4 = table(models, num2str(I(models)), Fs_s4, F_terms_s4, LHs_single_s4, complexity(models), exitflags_s4, ...
    time_rescaling_s4, 'VariableNames',{'s1_rank','orig_index','Fs_s4','F_terms_s4','log10_single_LHs','complexity','exitflags_s4','time_rescaling'});
for i=1:size(params_con,2)
    results_s4 = [results_s4 table(params_con_text(:,:,i), values_con_s4(:,i),'VariableNames',{sprintf('con_%d',i),sprintf('con_%d_v',i)})];
end
for i=1:size(params_reg,2)
    results_s4 = [results_s4 table(params_reg_text(:,:,i), values_reg_reordered(:,n_data*(i-1)+1:n_data*i),'VariableNames',{sprintf('reg_%d',i),sprintf('reg_%d_v',i)})];
end

results_good_s4 = results_s4(good_s4_in_good_s2,:);

%% Stage 4 parameter value analysis

if num_prefactors == 4 && LH_threshold_factor == 100 && max_complexity == 6 && length(good_s4) == 7
    I_clustering = [1 4 5 6 7 2 3]  % order of networks used in the paper
    I_clustering_fluxes = [1 5 2 4 6 3 0 7 0];  % when plotting the flux networks, this helps to get networks of the same group on top of each other
else
    [~, I_clustering] = sort(table2array(results_good_s4(:, 3)));
    I_clustering_fluxes = I_clustering;
end

plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
    [I(good_s2) Fs_s4 log10(time_rescaling_s4)], NaN, good_s4_in_good_s2, NaN, 'clustering manual', I_clustering, NaN, c_map, "s4");

plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
    [I(good_s2) Fs_s4], NaN, good_s4_in_good_s2, time_rescaling_s4, 'clustering manual', I_clustering, NaN, c_map, "s4");

plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
    [I(good_s2) Fs_s4], NaN, good_s4_in_good_s2, time_rescaling_s4, 'clustering manual', I_clustering(1:min(30,length(I_clustering))), NaN, c_map, "s4_top30");

plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
    [I(good_s2) Fs_s2(good_s2) Fs_s4], NaN, 1:length(good_s2), NaN, 'clustering manual', I_clustering_s2, NaN, c_map, "s4_all_candidates");  
  
%% Stage 4 Chromatin opening/closing rate ratios with children (with distance decay rates, slow)

A_reg_models = false(size(good_s4));
D_reg_models = false(size(good_s4));
opening_rates_ev = zeros(size(good_s4));
closing_rates_ev = zeros(size(good_s4));
opening_rates_dist = zeros(size(good_s4));
closing_rates_dist = zeros(size(good_s4));
tic
parfor i = 1:length(good_s4)
    r = good_s4(i);
    s = good_s4_in_good_s2(i);  % index after stage 2 i.e. in stage 4
    model_time_rescaling = time_rescaling_s4(s);
    
    params_reg_temp = params_reg(r,:);
    params_reg_temp = params_reg_temp(params_reg_temp > 0);
    if all(params_reg_temp <= 16)
        A_reg_models(i) = true;
    elseif all(params_reg_temp >= 17 & params_reg_temp <= 32)
        D_reg_models(i) = true;
    end

    data_set = 1;  % starting from repressed state using the dynamics of the activated state
    [opening_rates_ev(i), opening_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', false);
    
    data_set = 2;  % starting from repressed state using the dynamics of the activated state
    [closing_rates_ev(i), closing_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', false);
end
toc
shift_rate_ratios_ev = closing_rates_ev./opening_rates_ev;
shift_rate_ratios_dist = closing_rates_dist./opening_rates_dist;


figure('Position',[1500,200,400,600],'PaperUnits', 'centimeters','PaperSize', [10 15],'PaperPosition', [0 0 10 15])
defaultaxesfontsize = get(0,'defaultaxesfontsize');
myfontsize = 5;
set(0,'defaultaxesfontsize',myfontsize);

subplot(3,2,1)
hist(shift_rate_ratios_ev(A_reg_models))
title('Eigenvalue rate ratios (A-reg. models)')

subplot(3,2,3)
hist(shift_rate_ratios_dist(A_reg_models))
title('Distance decay rate ratios (A-reg. models)')

subplot(3,2,5)
plot(shift_rate_ratios_ev(A_reg_models), shift_rate_ratios_dist(A_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (A-reg. models)')
ylabel('Distance decay rate ratios (A-reg. models)')

subplot(3,2,2)
hist(shift_rate_ratios_ev(D_reg_models))
title('Eigenvalue rate ratios (D-reg. models)')

subplot(3,2,4)
hist(shift_rate_ratios_dist(D_reg_models))
title('Distance decay rate ratios (D-reg. models)')

subplot(3,2,6)
plot(shift_rate_ratios_ev(D_reg_models), shift_rate_ratios_dist(D_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (D-reg. models)')
ylabel('Distance decay rate ratios (D-reg. models)')

pause(0.5)
print(sprintf('img/shift_rate_ratios_s4.pdf'), '-dpdf')
close gcf
open(sprintf('img/shift_rate_ratios_s4.pdf'))
set(0,'defaultaxesfontsize',defaultaxesfontsize);

%% Stage 4 find child models of lower complexity models

results_final = results_good_s4;

temp = table2array(results_final(:,'complexity'));
low_complex_ranks = table2array(results_final(temp(:,1) < max(temp(:,1)), 's1_rank'));
all_ranks = table2array(results_final(:, 's1_rank'));

is_child_of = zeros(size(results_final,1),1);

if length(low_complex_ranks) < 1
    warning("no lower complexity models found")
else
    for i = 1:length(low_complex_ranks)
        for j = 1:length(all_ranks)
            l = low_complex_ranks(i);
            h = all_ranks(j);
            is_child = 0;
            if all(params_reg(l,:) == params_reg(h,:)) && l~=h % regulated parameter (names) have to agree
                params_temp = [params_reg(l,:) params_con(l,:)];
                params_temp = params_temp(1,params_temp>0);
                Wi_low = calc_value_indeces_in_W( sort(params_temp), C, 8, 1 );
                
                params_temp = [params_reg(h,:) params_con(h,:)];
                params_temp = params_temp(1,params_temp>0);
                Wi_high = calc_value_indeces_in_W( sort(params_temp), C, 8, 1 );
                
                if all(all((Wi_high > 0) - (Wi_low > 0) >= 0))
                    is_child = 1;
                    for p = unique(Wi_low(Wi_low > 0))'
                        Wi_p_low = Wi_low == p;
                        Wi_p_low_in_high_1 = Wi_high(Wi_p_low);
                        if any(Wi_p_low_in_high_1 == 0)
                            is_child = 0;
                        end
                        Wi_p_low_in_high_0 = Wi_high(~Wi_p_low);
                        if any(myismember(Wi_p_low_in_high_1, Wi_p_low_in_high_0))
                            is_child = 0;
                        end
                    end
                else
                    is_child = 0;
                end
            end
            if is_child
                is_child_of(j) = str2num(table2array(results_final(table2array(results_final(:,'s1_rank')) == l,2)));
            end
        end
    end
end

table_of_children = [results_final(:,1:2) table(num2str(is_child_of), 'VariableNames', {'is_child_of'})];
sortrows(table_of_children, 2)

good_final = good_s4(is_child_of == 0);
good_final_in_good_s2 = good_s4_in_good_s2(is_child_of == 0);

results_final = results_final(is_child_of==0,:);

if size(results_final,1) < size(results_good_s4,1)
    [~, I_clustering] = sort(table2array(results_final(:, 3)));
    I_clustering_fluxes = I_clustering;

    plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
        [I(good_s2) Fs_s4], NaN, good_final_in_good_s2, time_rescaling_s4, 'clustering manual', I_clustering, NaN, c_map, "s4_no_children");
    plot_parameter_value_analysis(values_con_s4, values_reg_s4, params_con(good_s2,:), params_reg(good_s2,:), C_text, ...
        [I(good_s2) Fs_s4], NaN, good_final_in_good_s2, time_rescaling_s4, 'clustering manual', I_clustering(1:min(30,length(I_clustering))), NaN, c_map, "s4_no_children_top30");
end

%% If there are too many models, move to the last step now

%% Stage 4 fluxes (adult models only)

models_to_plot = good_final_in_good_s2;
scaling_factors = time_rescaling_s4;

% activated and repressed state values in one plot
plot_matrix_value_analysis('reaction rates', W_gains_s4, C, C_text, [I(good_s2) Fs_s4 log10(scaling_factors)], NaN, models_to_plot, NaN, 'clustering manual', I_clustering, NaN, c_map);
[max_ADS, min_ADS] = plot_matrix_value_analysis('reaction rates', W_gains_s4, C, C_text, [I(good_s2) Fs_s4 log10(scaling_factors)], NaN, models_to_plot, scaling_factors, 'clustering manual', I_clustering, NaN, c_map);

max_ADS
min_ADS

if strcmp(data_text{3}, 'half-act: pho4[85-99] pho80D TATAm')
    data_text{3} = 'weakly act: pho4[85-99] pho80D TATAm';
end

% activated and repressed state log values in separate plots
plot_matrix_value_analysis_2('all fluxes', 'log', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
plot_matrix_value_analysis_2('all fluxes', 'log', fluxes_s4, NaN, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_matrix_value_analysis_2('all fluxes', 'linear', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_matrix_value_analysis_2('all fluxes', 'linear_cutoff', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);

plot_matrix_value_analysis_2('net fluxes', 'log', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
plot_matrix_value_analysis_2('net fluxes', 'log', fluxes_s4, NaN, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_matrix_value_analysis_2('net fluxes', 'linear', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_matrix_value_analysis_2('net fluxes', 'linear_cutoff', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);

plot_net_binned_fluxes('all fluxes', 'log', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_net_binned_fluxes('all fluxes', 'linear', fluxes_s4, scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);


% cycle fluxes
plot_cycle_fluxes(fluxes_s4, 'log', scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);
%plot_cycle_fluxes(fluxes_s4, 'linear', scaling_factors, C, C_text, models_to_plot, [I(good_s2)], NaN, 'clustering manual', I_clustering, NaN, c_map, data_text);

% (net) flux networks
if length(models_to_plot) <= 24
    plot_flux_networks(fluxes_s4, scaling_factors, ps_s4, models_to_plot, I(good_s2), I_clustering_fluxes)

    plot_net_flux_networks(fluxes_s4, scaling_factors, ps_s4, models_to_plot, I(good_s2), I_clustering_fluxes)

    plot_net_binned_flux_networks_2(fluxes_s4, scaling_factors, C, models_to_plot, I(good_s2), I_clustering_fluxes)
end

%% Stage 4 parameter value tables (adult models only)

rate_columns = [8 10:2:size(results_final,2)];
num_columns = [3 4 5 rate_columns];
results_output = results_final(I_clustering, :);

results_output(:,rate_columns) = varfun(@(var) log10(var), results_output(:, rate_columns));

results_output(:,num_columns) = varfun(@(var) round(var, 2, 'significant'), results_output(:, num_columns));

results_output

for i=1:size(results_output,1)  % reorder regulated values to rep. half-act. and act.
    t = results_output(i,'log10_single_LHs');
    results_output(i,'log10_single_LHs') = table(t{1,1}([2 3 1]));
    
    t = results_output(i,'reg_1_v');
    results_output(i,'reg_1_v') = table(t{1,1}([2 3 1]));
    
    t = results_output(i,'reg_2_v');
    results_output(i,'reg_2_v') = table(t{1,1}([2 3 1]));
end

writetable(results_output(:, [2 8:16]), "parameter_values_constitutive.csv")
writetable(results_output(:, [2 8 size(results_output,2)-[3 2]]), "parameter_values_regulated.csv")
writetable(results_output(:, [2 3 4 5 6]), "model_fit_results.csv")

%% Stage 4 Chromatin opening/closing rate ratios without children (with distance decay rates, slow)

A_reg_models = false(size(good_final));
D_reg_models = false(size(good_final));
opening_rates_ev = zeros(size(good_final));
closing_rates_ev = zeros(size(good_final));
opening_rates_dist = zeros(size(good_final));
closing_rates_dist = zeros(size(good_final));
tic
parfor i = 1:length(good_final)
    r = good_final(i);
    s = good_final_in_good_s2(i);  % index after stage 2 i.e. in stage 4
    model_time_rescaling = time_rescaling_s4(s);
    
    params_reg_temp = params_reg(r,:);
    params_reg_temp = params_reg_temp(params_reg_temp > 0);
    if all(params_reg_temp <= 16)
        A_reg_models(i) = true;
    elseif all(params_reg_temp >= 17 & params_reg_temp <= 32)
        D_reg_models(i) = true;
    end

    data_set = 1;  % starting from repressed state using the dynamics of the activated state
    [opening_rates_ev(i), opening_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', false);
    
    data_set = 2;  % starting from repressed state using the dynamics of the activated state
    [closing_rates_ev(i), closing_rates_dist(i), ~] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', false);
end
toc
shift_rate_ratios_ev = closing_rates_ev./opening_rates_ev;
shift_rate_ratios_dist = closing_rates_dist./opening_rates_dist;


figure('Position',[1500,200,400,600],'PaperUnits', 'centimeters','PaperSize', [10 15],'PaperPosition', [0 0 10 15])
defaultaxesfontsize = get(0,'defaultaxesfontsize');
myfontsize = 5;
set(0,'defaultaxesfontsize',myfontsize);

subplot(3,2,1)
hist(shift_rate_ratios_ev(A_reg_models))
title('Eigenvalue rate ratios (A-reg. models)')

subplot(3,2,3)
hist(shift_rate_ratios_dist(A_reg_models))
title('Distance decay rate ratios (A-reg. models)')

subplot(3,2,5)
plot(shift_rate_ratios_ev(A_reg_models), shift_rate_ratios_dist(A_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (A-reg. models)')
ylabel('Distance decay rate ratios (A-reg. models)')

subplot(3,2,2)
hist(shift_rate_ratios_ev(D_reg_models))
title('Eigenvalue rate ratios (D-reg. models)')

subplot(3,2,4)
hist(shift_rate_ratios_dist(D_reg_models))
title('Distance decay rate ratios (D-reg. models)')

subplot(3,2,6)
plot(shift_rate_ratios_ev(D_reg_models), shift_rate_ratios_dist(D_reg_models), '.', 'MarkerSize', 0.1)
xlabel('Eigenvalue rate ratios (D-reg. models)')
ylabel('Distance decay rate ratios (D-reg. models)')

pause(0.5)
print(sprintf('img/shift_rate_ratios_s4_no_children.pdf'), '-dpdf')
close gcf
open(sprintf('img/shift_rate_ratios_s4_no_children.pdf'))
set(0,'defaultaxesfontsize',defaultaxesfontsize);

%% Stage 4 model overview (adult models only)

for i = 1:length(good_final)
%for i=1

    for part = [1 2]

    r = good_final(i);  % rank index and index in stage 2
    s = good_final_in_good_s2(i);  % index after stage 2 i.e. in stage 4
    model_time_rescaling = time_rescaling_s4(s);

    defaultaxesfontsize = get(0,'defaultaxesfontsize');
    myfontsize = 5;
    set(0,'defaultaxesfontsize',myfontsize);
    rotation_mode = 'horizontal';
   
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [16*3/4 25],'PaperPosition', [-1 -1.8 19*3/4 28])
    
    barcolor = [0.9 0.9 0.9];

    ax1=subplot(5,2,1);
    draw_network_nice_2(C, C_text, params_con(r,:), params_reg(r,:), myfontsize-1,0, rotation_mode)
    ht=title(sprintf('Model %d, F = %2.2f, F terms = (%2.2f, %2.2f, %2.2f, %2.2f)',I(r),Fs_s4(s),F_terms_s4(s,1),F_terms_s4(s,2),F_terms_s4(s,3),F_terms_s4(s,4) ));
    ax1.Legend.Position = ax1.Legend.Position + [0.0075 0 0 0];
    if strcmp(rotation_mode,'vertical')
        ht.Position = ht.Position + [9 0 0];
        xlim_factor1 = 1.14;  % to enlarge the box to look similar to the ones below
        ax1.XLim = (ax1.XLim + 5) * xlim_factor1;
    else
        ht.Position = ht.Position + [0 9 0];
        ylim([-9 21]*1.07)
    end  
    
    ax3=subplot(5,2,3);
    data_set = 2;
    if part == 1
        plot_ps_with_data(ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    elseif part == 2
        plot_ps_with_data_all_stages(ps(r,:,data_set),ps_s2(:,data_set,r)',ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    end
    ht=title(sprintf('%s, log10 rel. LH = %2.2f -> %2.2f',data_text{data_set},LHs_single(r,data_set)-log10(LHs_data(data_set)),LHs_single_s4(s,data_set)-log10(LHs_data(data_set))));
    ht.Position = ht.Position + [0.5 0 0];
    
    ax2=subplot(5,2,5);
    data_set = 3;    
    if part == 1
        plot_ps_with_data(ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    elseif part == 2
        plot_ps_with_data_all_stages(ps(r,:,data_set),ps_s2(:,data_set,r)',ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    end
    ht=title(sprintf('%s, log10 rel. LH = %2.2f -> %2.2f',data_text{data_set},LHs_single(r,data_set)-log10(LHs_data(data_set)),LHs_single_s4(s,data_set)-log10(LHs_data(data_set))));
    ht.Position = ht.Position + [1 0 0];

    ax4=subplot(5,2,7);
    data_set = 1;
    if part == 1
        plot_ps_with_data(ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    elseif part == 2
        plot_ps_with_data_all_stages(ps(r,:,data_set),ps_s2(:,data_set,r)',ps_s4(s,:,data_set),configuration_count_data(:,data_set),'bars')
    end
    ht=title(sprintf('%s, log10 rel. LH = %2.2f -> %2.2f',data_text{data_set},LHs_single(r,data_set)-log10(LHs_data(data_set)),LHs_single_s4(s,data_set)-log10(LHs_data(data_set))));
    ht.Position = ht.Position + [0.5 0 0];
    
    ax5=subplot(5,2,9);
    marker_size = 3;
    nucl_order_N1 = [3 2 1];
    nucl_order_N2_N3 = [6 5 4 3 2 1];
    sticky_mut_occ = [[calc_occs_after_rate_manipulation(W_gains_s4{1}(:,:,s), log(prefactors_s4(s,1:4)), [-inf inf]); ...
                       calc_occs_after_rate_manipulation(W_gains_s4{2}(:,:,s), log(prefactors_s4(s,1:4)), [-inf inf])] ...
                      [calc_occs_after_rate_manipulation(W_gains_s4{1}(:,:,s), log(prefactors_s4(s,5:8)), [-inf inf]); ...
                       calc_occs_after_rate_manipulation(W_gains_s4{2}(:,:,s), log(prefactors_s4(s,5:8)), [-inf inf])]];

    b = bar(1-[occupancy_probs_s4(s,nucl_order_N1,2) occupancy_probs_s4(s,nucl_order_N1,1); sticky_mut_occ(nucl_order_N2_N3,1)'; sticky_mut_occ(nucl_order_N2_N3,2)';]');
    b(1).FaceColor = [0.2 0.7 0.1];
    b(2).FaceColor = [0.8 0.4 0.1];
    b(3).FaceColor = [0.95 0.95 0.05];
    hold on
    grid on
    base_fold_change_on_Boeger_data = 0;  % see MP_calc_objective_function_all_data
    if base_fold_change_on_Boeger_data
        [fix_occ_N2_N3_1, sd_fix_occ_N2_N3_1] = get_fixation_occupancy_data(occupancy_prob_data([2 3], 1), 1);
        [fix_occ_N2_N3_2, sd_fix_occ_N2_N3_2] = get_fixation_occupancy_data(occupancy_prob_data([2 3], 1), 2);
    else
        [fix_occ_N2_N3_1, sd_fix_occ_N2_N3_1] = get_fixation_occupancy_data(occupancy_probs_s4(s,[2 3],1)', 1);
        [fix_occ_N2_N3_2, sd_fix_occ_N2_N3_2] = get_fixation_occupancy_data(occupancy_probs_s4(s,[2 3],1)', 2);
    end
    % plot([0.775; 1.775; 2.775; 3.775; 4.775; 5.775],1-[occupancy_prob_data(:,1)' occupancy_prob_data(:,2)']','xk','MarkerSize',4)  % Boeger data
    errorbar([5 4 5.224 4.224]', 1-[fix_occ_N2_N3_1; fix_occ_N2_N3_2], [sd_fix_occ_N2_N3_1; sd_fix_occ_N2_N3_2]','.k','MarkerSize',8,'CapSize',2);
    xlim([0.5 6.5])
    ylim([0 1])
    legend('Unchanged model','Sticky mut. 1 prefactors','Sticky mut. 2 prefactors','Rescaled RE acc. data','Location','NorthEast')
    xlabel('           Repressed                                Activated')
    ylabel('Accessibility')
    ht=title('Accessibility change in sticky N-3 mutants');
    ht.Position = ht.Position + [0.25 0 0];
    ax5.XTickLabel = {"N-3" "N-2" "N-1" "N-3" "N-2" "N-1"};
    ax5.Legend.Position = ax5.Legend.Position + [-0.11 +0.01 0 0];
    
    [ts, Flag_assembly_prob_fun, ts_data, log2_F_M_ratios_N1_data, ~, ~, lambda, nucl, log2_F_M_ratio_sd] = get_Flag_Myc_data_Dion(0);
    t_subset=2:10:length(ts);
    t_indeces = ismember(ts(2:end), ts_data);
    
    % perturbed Flag parameters for both Flag plots
%     mycolor = [0.7 0.7 0.7];    
%     for Dion_param_index=1:2
%         [~, ~, Flag_probs_temp, cond_Flag_probs_temp, log2_F_M_ratio_error_temp, log2_F_N1_N2_ratio_error_temp, log2_F_M_ratio_diffs_error_temp] = ...
%             MP_calc_FM_error_and_time_scale_deterministic(W_gains_s4{2}(:,:,s), "", model_time_rescaling, Dion_param_index);
%         subplot(5,2,2)
%         log2_F_M_ratios_temp = squeeze(log2( cond_Flag_probs_temp(:,nucl) ./ (1 - cond_Flag_probs_temp(:,nucl)) ));
%         p_pert_2 = plot(ts(t_subset), log2_F_M_ratios_temp(t_subset) - mean(log2_F_M_ratios_temp(t_indeces)),'color',mycolor,'LineWidth',0.1); hold on
%         subplot(5,2,4)
%         p_pert_4 = plot(ts(2:10:end)', log2(Flag_probs_temp(1:10:end,1)./Flag_probs_temp(1:10:end,2)),'Color',mycolor,'LineWidth',0.1); hold on
%     end
    
    subplot(5,2,2)
    log2_F_M_ratios_temp = squeeze(log2( cond_Flag_probs_s4(:,nucl,s) ./ (1 - cond_Flag_probs_s4(:,nucl,s)) ));
    p1 = plot(ts(t_subset), log2_F_M_ratios_temp(t_subset) - mean(log2_F_M_ratios_temp(t_indeces)),'color',[0.1 0.7 0.1],'LineWidth',0.5); hold on
    p2 = errorbar(ts_data,log2_F_M_ratios_N1_data - mean(log2_F_M_ratios_N1_data),0.75*log2_F_M_ratio_sd*ones(4,1),0.75*log2_F_M_ratio_sd*ones(4,1),'.k','MarkerSize',8);
    ylim([-4 2])
    grid on
    xlabel('Time in h after lag')
    ylabel('Norm. log2 Flag / Myc (N-1)')
    %legend([p2,p1,p_pert_2],{'Dion et al. data','Model','Model with pert. Flag param.'},'location','southeast')
    legend([p2,p1],{'Dion et al. data','Model'},'location','southeast')
    ht=title('Flag over Myc ratio of N-1');
    ht.Position = ht.Position + [0.1 0 0];
    
    ax8 = subplot(5,2,4);
    h_p0 = plot(ts(2:10:end)', log2(Flag_probs_s4(1:10:end,1,s)./Flag_probs_s4(1:10:end,2,s)),'Color',[0.1 0.7 0.1],'LineWidth',0.5); hold on
    [t_data, ~, error_t_data, ~, log2_N1_N2_F_ratio_mean, log2_N1_N2_F_ratio_sd] = get_Flag_Myc_data_Rufiange();
    h_data = errorbar(t_data,log2_N1_N2_F_ratio_mean,log2_N1_N2_F_ratio_sd,log2_N1_N2_F_ratio_sd,'.k','MarkerSize',8);
    h_limit = plot([ts(1); ts(end)], ones(2,1) * log2(occupancy_probs_s4(s,1,2)/occupancy_probs_s4(s,2,2)), 'k:');
    xlabel('Time after lag in h')
    ylabel('log2 Flag ratio (N-1 / N-2)')
    t_max = max(ts);
    ylim([-2 1])
    xlim([0 t_max])  
    grid on
    %legend([h_data,h_p0(1),p_pert_4,h_limit],{'Rufiange et al. data','Model','Model with pert. Flag param.','Model steady state ratio'},'Location','NorthEast');
    legend([h_data,h_p0(1),h_limit],{'Rufiange et al. data','Model','Model steady state ratio'},'Location','NorthEast');
    ht=title('Flag ratio of N-1 over N-2');
    ht.Position = ht.Position + [0.1 0 0];
    
    ax6=subplot(5,2,6);
    if part == 1
        data_set = 2;
        net_fluxes_temp = calc_net_fluxes(fluxes_s4{data_set}(:,:,s)*model_time_rescaling);
        fmax_temp = max(max(net_fluxes_temp));
        draw_flux_network(net_fluxes_temp, 0, fmax_temp, ps_s4(s,:,data_set), myfontsize-1, 0, rotation_mode)
        ht=title(sprintf('Net fluxes, %s, max = %.2f / h',data_text{data_set}, fmax_temp));
        if strcmp(rotation_mode,'vertical')
            ht.Position = ht.Position + [3 0 0];
            ax6.XLim = (ax6.XLim + 0) * xlim_factor1;
        else
            ht.Position = ht.Position + [0 3 0];
            ylim([-15 15]*1.07)
        end
    elseif part == 2
        data_set = 1;  % starting from repressed state using the dynamics of the activated state
        [shift_rate, shift_rate_with_distance_measure, shift_probs] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', true);
        title(sprintf('rep. to act. distance decay rate %f/h\nslowest eigenrate %f/h', shift_rate_with_distance_measure, shift_rate))
    end
    

    ax7=subplot(5,2,8);
    if part == 1
        data_set = 1;
        net_fluxes_temp = calc_net_fluxes(fluxes_s4{data_set}(:,:,s)*model_time_rescaling);
        fmax_temp = max(max(net_fluxes_temp));
        draw_flux_network(net_fluxes_temp, 0, fmax_temp, ps_s4(s,:,data_set), myfontsize-1, 0, rotation_mode)
        ht=title(sprintf('Net fluxes, %s, max = %.2f / h',data_text{data_set}, fmax_temp));
        if strcmp(rotation_mode,'vertical')
            ht.Position = ht.Position + [3 0 0];
            ax6.XLim = (ax6.XLim + 0) * xlim_factor1;
        else
            ht.Position = ht.Position + [0 3 0];
            ylim([-15 15]*1.07)
        end
    elseif part == 2
        data_set = 2;
        [shift_rate, shift_rate_with_distance_measure, shift_probs] = calc_chromatin_shift_rate_with_distance_measure(W_gains_s4{data_set}(:,:,s) * model_time_rescaling, ps_s4(s,:,setdiff([1 2], data_set))', true);
        title(sprintf('act. to rep. distance decay rate %f/h\nslowest eigenrate %f/h', shift_rate_with_distance_measure, shift_rate))
    end
    
    ax10 = subplot(5,2,10);
    fluxes_mat = zeros(2,32);  % number of possible reactions
    for data_set = 1:2
        n = 1;
        for p=1:length(C)
            if size(C{p},1)==1  % single reaction
                fluxes_mat(data_set,n) = fluxes_s4{data_set}(C{p}(1,1),C{p}(1,2),s) * model_time_rescaling;
                n = n+1;
            end
        end
    end
    binned_fluxes_mat = zeros(2,10);
    for n=1:6
        binned_fluxes_mat(:,n) = sum(fluxes_mat(:,4*(n-1)+(1:4)),2);
    end
    for n=7:10
        binned_fluxes_mat(:,n) = sum(fluxes_mat(:,12+2*(n-1)+(1:2)),2);
    end
    % entries of binned_fluxes: A1  A2  A3  D1  D2  D3  S2>1 S2>3 S1>2 S3>2
    net_binned_fluxes_mat = binned_fluxes_mat - binned_fluxes_mat(:,[4:6 1:3 9 10 7 8]);
    fmax = max(max(abs(net_binned_fluxes_mat)));
    net_binned_fluxes_mat(net_binned_fluxes_mat<1e-2*fmax) = 0;
    draw_nucl_pov_fluxes_2(net_binned_fluxes_mat, myfontsize)
    title(sprintf('Site-centric net fluxes (rep./act.), max = %.2f / h', fmax))

    pause(0.5)
    if part == 1
        print(sprintf('img/models/model_overview_%d_part_1',I(r)),'-dpdf')
        close gcf
        %open(sprintf('img/models/model_overview_%d_part_1.pdf',I(r)))
    elseif part == 2
        print(sprintf('img/models/model_overview_%d_part_2',I(r)),'-dpdf')
        close gcf
        %open(sprintf('img/models/model_overview_%d_part_2.pdf',I(r)))
    end
    
    set(0,'defaultaxesfontsize',defaultaxesfontsize);
    
    end  % part 1 and 2
end

%% Find models used in Boeger's paper

model_D_params = [1 17 34];
model_F_params = [1 17 19 21 34];

C_text(model_D_params)
C_text(model_F_params)

model_D = [];  % should be more than one, because of different regulation parameters
model_F = [];
model_D_relatives = [];
model_F_relatives = [];

for i=1:size(params,1)
    if all(myismember(model_D_params, params(i,:)))
        model_D_relatives = [model_D_relatives i];
        if all(myismember(params(i,:), [0 model_D_params]))
            model_D = [model_D i];
        end
    end
    if all(myismember(model_F_params, params(i,:)))
        model_F_relatives = [model_F_relatives i];
        if all(myismember(params(i,:), [0 model_F_params]))
            model_F = [model_F i];
        end
    end
end

params(model_D,:)
params(model_F,:)

LH_bound

LHs(model_D)
LHs(model_F)

max(LHs(model_D_relatives))
max(LHs(model_F_relatives))

model_D_relatives_s1_good = model_D_relatives(LHs(model_D_relatives) >= LH_bound);
model_F_relatives_s1_good = model_F_relatives(LHs(model_F_relatives) >= LH_bound);

model_D_relatives_s2_good = model_D_relatives_s1_good(myismember(model_D_relatives_s1_good, good_s2));
model_F_relatives_s2_good = model_F_relatives_s1_good(myismember(model_F_relatives_s1_good, good_s2));

model_D_relatives_s4_good = model_D_relatives_s2_good(myismember(model_D_relatives_s2_good, good_s4));
model_F_relatives_s4_good = model_F_relatives_s2_good(myismember(model_F_relatives_s2_good, good_s4));

%% Sensitivity analysis (using eigenvectors of the Hessian)

% uses DERIVESTsuite for gradient and hessian calculation

max_complexity_with_timescale = max_complexity+1;

N_models = length(good_final);

err_minus = NaN * zeros(max_complexity_with_timescale+2,max_complexity_with_timescale,N_models);
err_plus = NaN * zeros(max_complexity_with_timescale+2,max_complexity_with_timescale,N_models);
err_minus_final = NaN * zeros(N_models,max_complexity_with_timescale);
err_plus_final = NaN * zeros(N_models,max_complexity_with_timescale);

err_minus_still_ok = NaN * zeros(max_complexity_with_timescale+2,max_complexity_with_timescale,N_models);
err_plus_still_ok = NaN * zeros(max_complexity_with_timescale+2,max_complexity_with_timescale,N_models);
err_minus_final_still_ok = NaN * zeros(N_models,max_complexity_with_timescale);
err_plus_final_still_ok = NaN * zeros(N_models,max_complexity_with_timescale);

boundary_values = NaN * zeros(N_models,1);
boundary_prefactors = NaN * zeros(N_models,1);
grad_norms = NaN * zeros(N_models,1);
low_evs = NaN * zeros(N_models,4);
log10_time_scale_plus = NaN * zeros(N_models,1);  % scalar product of corresponding eigenvector with unit vector (if in log space)
log10_time_scale_minus = NaN * zeros(N_models,1);

sloppy_boundary_params = NaN * zeros(N_models,1);
sloppy_boundary_params_out = NaN * zeros(N_models,1);
sloppy_boundary_params_in = NaN * zeros(N_models,1);
sloppy_internal_params = NaN * zeros(N_models,1);
sloppy_internal_params_plus = NaN * zeros(N_models,1);
sloppy_internal_params_minus = NaN * zeros(N_models,1);

better_minimum_found = zeros(N_models,1);
better_values_s4 = NaN * zeros(N_models,max_complexity_with_timescale);

M = zeros(n_data,1);
for k = 1:n_data
    M(k) = multinomial(configuration_count_data(:,k));
end
log10_M = log10(M);
log10_LH_part2 = -sum(log10_M) + sum(log10(LHs_data));

%for i = 1:5  % for testing
for i = 1:length(good_final)

    tic
    
    r = good_final(i);  % rank index and index in stage 2
    s = good_final_in_good_s2(i);  % index after stage 2 i.e. in stage 4
    
    disp(sprintf("\n%d  %d\n",i, I(r)))
    
    values_con_temp = values_con_s4(s,:);
    values_con_temp = values_con_temp(values_con_temp > 0);
    values_reg_temp = values_reg_s4(s,:,:);
    values_reg_temp = values_reg_temp(values_reg_temp > 0);
    values = [values_con_temp values_reg_temp(:)'];

    params_con_temp = params_con(r,:);
    params_con_temp = params_con_temp(params_con_temp > 0);
    params_reg_temp = params_reg(r,:);
    params_reg_temp = params_reg_temp(params_reg_temp > 0);
    params_all = [params_con_temp params_reg_temp params_reg_temp params_reg_temp];
    
    log_values = [log(values) log(time_rescaling_s4(s))];
    this_complexity_total = length(values) + 8;
    log_rate_bounds = log(rate_bounds);

    p = find(params_all == fixed_param,1);
    log_values = log_values([1:p-1, p+1:end]);

    value_indeces_in_W = calc_value_indeces_in_W([params_con_temp params_reg_temp], C, size(configuration_count_data,1), 0);
 
    f = @(log_values) MP_calc_objective_function_all_data(4, p, length(values_con_temp), length(params_reg_temp), value_indeces_in_W, ...
        [log_values(1:end-1) log(prefactors_s4(s,:))], 0, configuration_count_data, min_S_rate, log10_LH_part2, num_prefactors, exp(log_values(end)));

    
    if 0  % calc hessian
        [grad,err,finaldelta] = gradest(f, log_values);
        [hess,err] = hessian(f, log_values);

        gradients{i} = grad;
        hessians{i} = hess;
    else
        grad_all = derivest_gradients_s4(s, 1:this_complexity_total);
        hess_all = derivest_hessians_s4(1:this_complexity_total,1:this_complexity_total,s);
        
        [V_all, D_all] = eig(hess_all);
        
        grad = grad_all([1:length(log_values)-1 end]);
        hess = hess_all([1:length(log_values)-1 end],[1:length(log_values)-1 end]);
    end
    
    C_text(params_all)
    log_values
    grad
    grad_norm = norm(grad)
    [V, D] = eig(hess)
    
    
    EVs = diag(D);
            
    % find bounds for parameter values that change the objective function by dF in direction of eigenvectors
    %dF = log10(5/4);  % log10(5/4) corresponds to a likelihood decrease with factor 4/5
    %sloppy_dx = log(1.5);  % when a value error is called sloppy
    dF = log10(2);  % log10(2) corresponds to a likelihood decrease with factor 0.5
    sloppy_dx = log(2);  % when a value error is called sloppy
    
    % eps will be the scalar factor in front of the normalized eigenvector. 
    % eps needs to be high enough to reach dF even for directions that change almost only one very sloppy parameter, but find the finite error of the other parameters.
    % eps needs to be fine grained enough in the first steps to find the right scalar factor that barely leads to dF for non-sloppy directions 
    % (if we not interested in errors but only if a parameter is sloppy, then we can start eps with sloppy_dx)
    %eps = sloppy_dx/50 * 10.^([0:0.1:2 2.2:0.2:4 4.4:0.4:7]);  % test steps in log rate space
    eps = 10.^(-2:0.005:4);

    
    %eps_high = diff(log_rate_bounds)*sqrt(length(log_values));  % then we definitely cross a rate bound in any unit direction
    %eps = [sloppy_dx/10:sloppy_dx/10:sloppy_dx*5 (sqrt(sloppy_dx*6):sloppy_dx:sqrt(eps_high)).^2 2*eps_high*(1:6).^3]; 
    
    max_boundary_crossing = log(100);  % to keep rates from exploding: allow crossing the boundaries by a factor of 100, but not more (can be very slow to solve the ode with extremly high rates)
    
    V_test = [eye(length(log_values)) V grad'/grad_norm];
    
    F_min = f(log_values);
    for v = 1:size(V_test,2)
        dF_1 = zeros(1,length(eps));
        dF_2 = zeros(1,length(eps));
        vec = V_test(:,v)';
        
        eps_1 = Inf;
        eps_1_still_ok = 0;
        for j = 1:length(eps)
            log_values_1 = log_values + eps(j) * vec;
            log_values_1 = min([log_values_1; (log_rate_bounds(2)+max_boundary_crossing) * ones(1,length(log_values_1))]);
            log_values_1 = max([log_values_1; (log_rate_bounds(1)-max_boundary_crossing) * ones(1,length(log_values_1))]);
            dF_1(j) = f(log_values_1) - F_min;
            if dF_1(j) < better_minimum_found(i) && all(log_values_1 >= log_rate_bounds(1)) && all(log_values_1 <= log_rate_bounds(2))  % found a better minimum inside boundary!
                better_minimum_found(i) = dF_1(j);
                better_values_s4(i,1:length(log_values_1)) = exp(log_values_1);
            end
             
            if dF_1(j) > dF
                eps_1 = eps(j);
                break
            end
            eps_1_still_ok = eps(j);
        end
        
        eps_2 = Inf;
        eps_2_still_ok = 0;
        for j = 1:length(eps)
            log_values_2 = log_values - eps(j) * vec;
            log_values_2 = min([log_values_2; (log_rate_bounds(2)+max_boundary_crossing) * ones(1,length(log_values_2))]);
            log_values_2 = max([log_values_2; (log_rate_bounds(1)-max_boundary_crossing) * ones(1,length(log_values_2))]);
            dF_2(j) = f(log_values_2) - F_min;
            if dF_2(j) < better_minimum_found(i) && all(log_values_2 >= log_rate_bounds(1)) && all(log_values_2 <= log_rate_bounds(2))  % found a better minimum inside boundary!
                better_minimum_found(i) = dF_2(j);
                better_values_s4(i,1:length(log_values_2)) = exp(log_values_2);
            end
            if dF_2(j) > dF
                eps_2 = eps(j);
                break
            end
            eps_2_still_ok = eps(j);
        end
   
        err_1 = eps_1 * vec;
        err_2 = -eps_2 * vec;
               
        err_plus(v,vec >= 0, i) = err_1(vec >= 0);
        err_plus(v,vec < 0, i) = err_2(vec < 0);
        
        err_minus(v,vec >= 0, i) = err_2(vec >= 0);
        err_minus(v,vec < 0, i) = err_1(vec < 0);
        
        err_1_still_ok = eps_1_still_ok * vec;
        err_2_still_ok = -eps_2_still_ok * vec;
        
        err_plus_still_ok(v,vec >= 0, i) = err_1_still_ok(vec >= 0);
        err_plus_still_ok(v,vec < 0, i) = err_2_still_ok(vec < 0);
        
        err_minus_still_ok(v,vec >= 0, i) = err_2_still_ok(vec >= 0);
        err_minus_still_ok(v,vec < 0, i) = err_1_still_ok(vec < 0);
    end
   
    err_minus_temp = err_minus(:,:,i);
    err_plus_temp = err_plus(:,:,i);
    
    [err_minus_final_temp, I_minus] = min(err_minus_temp,[],1);
    [err_plus_final_temp, I_plus] = max(err_plus_temp,[],1);
    
    err_minus_final(i,:) = err_minus_final_temp;
    err_plus_final(i,:) = err_plus_final_temp;
    
    err_minus_temp_still_ok = err_minus_still_ok(:,:,i);
    err_plus_temp_still_ok = err_plus_still_ok(:,:,i);
    
    [err_minus_final_temp_still_ok, I_minus_still_ok] = min(err_minus_temp_still_ok,[],1);
    [err_plus_final_temp_still_ok, I_plus_still_ok] = max(err_plus_temp_still_ok,[],1);
    
    err_minus_final_still_ok(i,:) = err_minus_final_temp_still_ok;
    err_plus_final_still_ok(i,:) = err_plus_final_temp_still_ok;
    
    log_values
    err_minus_final_temp
    err_plus_final_temp
    
    log_boundary_area = log(1.1);
       
    params_at_lower_boundary = log_values < log_rate_bounds(1) + log_boundary_area & [ones(1,length(log_values)-1) 0];  % do not include time scale here
    params_at_upper_boundary = log_values > log_rate_bounds(2) - log_boundary_area & [ones(1,length(log_values)-1) 0];
    params_at_boundary = params_at_lower_boundary | params_at_upper_boundary;
    
    sloppy_boundary_params(i) = sum(abs(err_minus_final_temp(params_at_boundary)) >= sloppy_dx | abs(err_plus_final_temp(params_at_boundary)) >= sloppy_dx); 
    sloppy_boundary_params_out(i) = sum(abs([err_minus_final_temp(params_at_lower_boundary) err_plus_final_temp(params_at_upper_boundary)]) >= sloppy_dx);
    sloppy_boundary_params_in(i) = sum(abs([err_minus_final_temp(params_at_upper_boundary) err_plus_final_temp(params_at_lower_boundary)]) >= sloppy_dx);
    
    sloppy_internal_params(i) = sum(abs(err_minus_final_temp(~params_at_boundary)) >= sloppy_dx | abs(err_plus_final_temp(~params_at_boundary)) >= sloppy_dx);
    sloppy_internal_params_plus(i) = sum(abs(err_plus_final_temp(~params_at_boundary)) >= sloppy_dx);
    sloppy_internal_params_minus(i) = sum(abs(err_minus_final_temp(~params_at_boundary)) >= sloppy_dx);
    
    log10_time_scale_plus(i) = log10(exp(1)) * err_plus_final_temp(end);
    log10_time_scale_minus(i) = log10(exp(1)) * err_minus_final_temp(end);
    
    EVs_all = diag(D);
    
    boundary_values(i) = sum(params_at_lower_boundary) + sum(params_at_upper_boundary);
    boundary_prefactors(i) = sum(log(prefactors_s4(s,:)) < log(1/5) + log_boundary_area) + sum(log(prefactors_s4(s,:)) > log(5) - log_boundary_area);
    grad_norms(i) = grad_norm;
    low_evs(i,:) = EVs_all(1:4);
    
    toc
end
log10_time_scale = log10(table2array(results_final(:,8)));
results_sensitivity_analysis = [results_final(:,[2 3]) table(better_minimum_found, grad_norms, low_evs(:,1:3), sloppy_internal_params, boundary_values, sloppy_boundary_params, boundary_prefactors, ...
    log10_time_scale, log10_time_scale_plus, log10_time_scale_minus)];

results_sensitivity_analysis(I_clustering,:)

mean(log10(grad_norms))
mean(sloppy_internal_params)
mean(log10_time_scale_plus)

log10_err_plus_final = err_plus_final * log10(exp(1));
log10_err_minus_final = err_minus_final * log10(exp(1));

csvwrite("parameter_errors_plus.csv", round(log10_err_plus_final(I_clustering,:), 4, 'significant'))
csvwrite("parameter_errors_minus.csv", round(log10_err_minus_final(I_clustering,:), 4, 'significant'))

% Save final data

save('final_results.mat', '-regexp', '^(?!(values_reg_mat|values_con_mat|exitflags|exitflags_best_LH|exitflags_sorted|LHs_mat|nonunique_best_values_ind|nonunique_LHs_ind|error_param_missing|error_param_not_often_enough|params_con_setup|params_reg_setup)$).')
