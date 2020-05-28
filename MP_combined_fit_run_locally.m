% This script sets up the regulated on-off-slide models (i.e. their process/parameter combinations) and runs the first stage of the fit procedure
% MP -> multi process/parameter
% combined -> simultaneous fit to three configuration data sets

clear

set(0,'defaultfigureposition',[1700,200,800,1200]')
set(0,'defaultfigurecolor','white')
set(0,'defaultaxesfontsize',10)
warning('WARNING TEST OK, DEFAULTS SET')

addpath data_functions

%% Setup models and data

% load experimental data: 1 'rep: pho4D pho80D TATAm', 2 'half-act: pho4[85-99] pho80D TATAm', 3 'act: pho80D TATAm', 4 'rep: wt', 5 'act: pho80D', 6 'rep: pho2D'
strain_indeces = [5 4 2];

preset_params = [1 17];  % these parameters have to be in the model; they can be regulated (useful options are [1 17] and [1 17 33]); global assembly: 1, disassembly: 17, sliding: 33
params_never_reg = [];
fixed_param = 1;  % has to be in preset_params (default: 1 for global assembly in activated state)
max_complexity = 6;  % total number of real free parameters without time scale (constitutive params + 2x regulated params - 1)
max_num_params_con = [];
max_num_params_reg = [];
rate_bounds = [1e-2 1e2];
min_S_rate = 0;
%min_S_rate = min(rate_bounds);

n_data = length(strain_indeces);
[ C, C_text, C_text_mat, params_con_setup, params_reg_setup ] = MP_combined_fit_param_setup( preset_params, params_never_reg, max_complexity, max_num_params_con, max_num_params_reg, n_data );

[ configuration_count_data, configuration_prob_data, occupancy_prob_data, LHs_data, data_text ] = get_configuration_data( strain_indeces );

if 0  % limit number of models for test purposes
    warning('Number of models manually LIMITED!')
    params_con_setup = params_con_setup(1:10000,:);
    params_reg_setup = params_reg_setup(1:10000,:);
end

%% Check models

% error_param_missing indicates models where a parameter is not used at all
% error_param_not_often_enough indicates models where a parameter can be replaced by another (effective model duplicates)
% both are zero if MP_combined_fit_param_setup sets up the models correctly

params_temp = zeros(size(params_con_setup,1),size(params_con_setup,2)+size(params_reg_setup,2));
value_indeces_in_W = zeros(8,8,size(params_con_setup,1));
error_param_missing = zeros(size(params_con_setup,1),max_complexity+1);
error_param_not_often_enough = zeros(size(params_con_setup,1),max_complexity+1);

for i = 1:size(params_con_setup,1)
    
    params_temp(i,:) = [params_con_setup(i,:) params_reg_setup(i,:)];
    params_temp_2 = params_temp(i,params_temp(i,:)>0);
    [value_indeces_in_W(:,:,i), error_param_missing_temp, error_param_not_often_enough_temp] = calc_value_indeces_in_W( params_temp_2, C, 8, 1 );
    error_param_missing(i,1:length(error_param_missing_temp)) = error_param_missing_temp;
    error_param_not_often_enough(i,1:length(error_param_not_often_enough_temp)) = error_param_not_often_enough_temp;
    
end

models_param_missing = find(any(error_param_missing>0,2));
models_param_not_often_enough = find(any(error_param_not_often_enough>0,2));

num_models_param_missing = length(models_param_missing)
num_models_param_not_often_enough = length(models_param_not_often_enough)

%% Run combined optimization multiple times

run_repeats = 2;
iteration_factor = 1e2;  % 10 is usually enough to avoid exitflag 0 when tolerace_factor = 1, 1e3 increases the needed time too much when tolerance_factor = 1e-4
tolerance_factor = 1e-4;  % lower means lower tolerances

%fmincon_options = optimoptions('fmincon');  % uses default algorithm 'interior-point', takes at least 1.5 to 2 times longer without being more robust (in default setting even less accurate), however sometimes finds better LH-values
fmincon_options = optimoptions('fmincon','Algorithm','sqp');

fmincon_options.Display = 'none';  % Display set to iter gives progress output

fmincon_options.MaxIterations = iteration_factor * fmincon_options.MaxIterations;
if strcmp(fmincon_options.Algorithm,'sqp')
    fmincon_options.MaxFunctionEvaluations = iteration_factor * 100*max_complexity;
else
    fmincon_options.MaxFunctionEvaluations = iteration_factor * fmincon_options.MaxFunctionEvaluations;
end

fmincon_options.ConstraintTolerance = tolerance_factor * fmincon_options.ConstraintTolerance;
fmincon_options.OptimalityTolerance = tolerance_factor * fmincon_options.OptimalityTolerance;
fmincon_options.StepTolerance = tolerance_factor * fmincon_options.StepTolerance;

options = fmincon_options;

gcp;  % start parpool for parfor loop in MP_combined_fit if not running

tic
for run_number=1:run_repeats  % this loop can be parallized on a cluster to speed things up when doing many repeats
    run_number
    [ LHs_mat(:,run_number), values_con_mat(:,:,run_number), values_reg_mat(:,:,:,run_number), exitflags(:,run_number) ] = MP_combined_fit( configuration_count_data, C, params_con_setup, params_reg_setup, fixed_param, rate_bounds, min_S_rate, options);
end
calculation_time = toc
calculation_time_h = calculation_time / 3600;

if ~isreal(LHs_mat)
    complex_indeces = ~arrayfun(@(x) isreal(x),LHs_mat);
    warning('Complex LH values found!')
end

save temp_results.mat

%% Continue with MP_combined_fit_analysis_all_data

% run MP_combined_fit_analysis_all_data.m
