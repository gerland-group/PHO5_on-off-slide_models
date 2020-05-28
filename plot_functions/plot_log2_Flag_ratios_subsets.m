function plot_log2_Flag_ratios_subsets(Flag_probs, models, no_sliding_models, no_config_param_models, lower_complexity_models, t_max, occupancy_prob_data, Fs_s4, F_max)

    ts = get_Flag_Myc_data_Dion();

    if ts(1) == 0
        ts = ts(2:end);
    end
    
    if isfinite(t_max) && t_max < max(max(ts))
        i_max = find(ts >= t_max,1);
        ts = ts(1:i_max);
        Flag_probs = Flag_probs(1:i_max,:,:);
    else
        t_max = max(max(ts));
    end     

    if length(ts) > 250
        step = round(length(ts) / 250);
        t_subset = 1:step:length(ts);
        ts = ts(t_subset);
        Flag_probs = Flag_probs(t_subset,:,:);
    end
    
    if size(ts,1) == 1
        ts = ones(size(Flag_probs,3),1)*ts;
    end
    
    defaultaxesfontsize = get(0,'defaultaxesfontsize');
    myfontsize = 6;
    set(0,'defaultaxesfontsize',myfontsize);
    
    no_sliding_models = find(ismember(models,no_sliding_models));
    no_config_param_models = find(ismember(models,no_config_param_models));
    lower_complexity_models = find(ismember(models,lower_complexity_models));
    
    [t_data, ~, ~, ~, log2_N1_N2_F_ratio_mean, log2_N1_N2_F_ratio_sd] = get_Flag_Myc_data_Rufiange();

    log2_N1_N2_F_ratios = log2(squeeze(Flag_probs(:,1,:) ./ Flag_probs(:,2,:)));
    N = size(log2_N1_N2_F_ratios,2);
  
    figure('Position',[2000,100,320,320]','color','white','PaperUnits','centimeters','PaperSize',[8.5 7],'PaperPosition',[-0.2 0 9.3 7.3])

    good_models = find(Fs_s4<=F_max);
    
    models = setdiff(1:N,good_models');
    h_rest = plot(ts(models,:)', log2_N1_N2_F_ratios(:,models),'Color',[0.5 0.5 0.5],'LineWidth',0.1); hold on;
    
    models = good_models;
    h_good = plot(ts(models,:)', log2_N1_N2_F_ratios(:,models),'Color',[0.1 0.7 0.1],'LineWidth',0.1);
    
    %models = no_sliding_models;
    %h_no_sliding = plot(ts(models,:)', log2_N1_N2_F_ratios(:,models),'Color',[0.9 0.1 0.1],'LineWidth',0.1);

    %models = no_config_param_models;
    %h_no_config_param = plot(ts(models,:)', log2_N1_N2_F_ratios(:,models),'Color',[0.1 0.1 0.9],'LineWidth',0.1);

    %models = lower_complexity_models;
    %h_lower_complexity = plot(ts(models,:)', log2_N1_N2_F_ratios(:,models),'Color',[0.3 0.7 0.3],'LineWidth',0.1);
    
    h_limit = plot([ts(1,1); ts(1,end)], ones(2,1) * log2(occupancy_prob_data(1,2)/occupancy_prob_data(2,2)), 'k:');
   
    h_data = errorbar(t_data,log2_N1_N2_F_ratio_mean,log2_N1_N2_F_ratio_sd,log2_N1_N2_F_ratio_sd,'b.','CapSize',3,'MarkerSize',8);

    xlabel('Time after lag in h')
    ylabel('log2 Flag ratio (N-1 over N-2)')
    ylim([-2.5 1])
    xlim([0 t_max])  
    grid on

    if 0 && length(h_no_sliding) > 0 && length(h_no_config_param) > 0 && length(h_lower_complexity) > 0
        h = legend([h_data,h_no_sliding(1),h_no_config_param(1),h_lower_complexity(1),h_rest(1),h_limit],...
            'Rufiange et al. data','No sliding','No config.-spec. param.','Lower complexity','Other models','Boeger et al. nucl. occ. ratio','Location','SouthEast');
    else
        h = legend([h_data,h_good(1),h_rest(1),h_limit],...
            'Rufiange et al. data','Good models','Other models','Boeger et al. nucl. occ. ratio','Location','SouthEast');
    end

    pause(1)
    filename = sprintf('img/log2_Flag_ratio_subsets_t_max_%2.0f.pdf',10*t_max);
    print('-dpdf','-painters',filename)
    close gcf
    open(filename)

    set(0,'defaultaxesfontsize',defaultaxesfontsize);
    
end
