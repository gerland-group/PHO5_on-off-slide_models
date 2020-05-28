function [] = plot_log2_F_M_ratio_diffs(cond_Flag_probs, Fs_s4, F_max)
    
    [ts, Flag_assembly_prob_fun, ts_data, log2_F_M_ratios_N1_data, err_plus, err_minus, lambda, nucl, log2_F_M_ratio_sd] = get_Flag_Myc_data_Dion();

    if ts(1) == 0
        ts = ts(2:end);
    end
    
    if length(ts) > 250
        step = round(length(ts) / 250);
        t_subset = 1:step:length(ts);
    end
    
    t_indeces = find(ismember(ts, ts_data));
    
    i_norm = [1 2 3 4];
    is_rest = [1 2 3 4];
 
    log2_F_M_ratios_N1 = squeeze(log2( cond_Flag_probs(:,nucl,:) ./ (1 - cond_Flag_probs(:,nucl,:)) ));

    defaultaxesfontsize = get(0,'defaultaxesfontsize');
    myfontsize = 6;
    set(0,'defaultaxesfontsize',myfontsize);

    figure('Position',[2000,100,320,320]','color','white','PaperUnits','centimeters','PaperSize',[8.5 7],'PaperPosition',[-0.2 0 9.3 7.3])

    for i = 1:size(log2_F_M_ratios_N1,2)
        if Fs_s4(i) <= F_max
            mycolor = [0.1 0.7 0.1];
            p1=plot(ts(t_subset), log2_F_M_ratios_N1(t_subset,i) - mean(log2_F_M_ratios_N1(t_indeces(i_norm),i)),'color',mycolor,'LineWidth',0.1); hold on
        else
            mycolor = [0.5 0.5 0.5];
            p0=plot(ts(t_subset), log2_F_M_ratios_N1(t_subset,i) - mean(log2_F_M_ratios_N1(t_indeces(i_norm),i)),'color',mycolor,'LineWidth',0.1); hold on
        end
    end

    p2=errorbar(ts_data(is_rest),log2_F_M_ratios_N1_data(is_rest) - mean(log2_F_M_ratios_N1_data(i_norm)),log2_F_M_ratio_sd*ones(length(is_rest),1),log2_F_M_ratio_sd*ones(length(is_rest),1),'.b','CapSize',3,'MarkerSize',8);
    ylim([-4 2])
    grid on
    xlabel('Time in h after lag')
    ylabel('log2 Flag / Myc (N-1) - mean over the four measurement times')
    legend([p1,p0,p2],{'Good models', 'Other models', 'Dion et al. data'},'location','southeast')

    pause(1)
    filename = sprintf('img/log2_F_M_ratio_diffs_N%d', nucl);
    print('-dpdf','-painters',filename)
    close gcf
    open([filename '.pdf'])

    set(0,'defaultaxesfontsize',defaultaxesfontsize);

end

