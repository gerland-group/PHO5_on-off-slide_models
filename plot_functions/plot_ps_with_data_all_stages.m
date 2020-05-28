function [ ] = plot_ps_with_data( ps_s1, ps_s2, ps_s4, data, mode )
    
    ps_data = data ./ sum(data);
    stds = sqrt(data.*(sum(data)-data)/(sum(data)^2*(sum(data)+1)));  % std of Dirichlet distribution 
    
    barcolor = [0.9 0.9 0.9];
    hbar = bar(ps_data(:,1),'FaceColor',barcolor); hold on

    if strcmp(mode,'bars')
        herror = errorbar(ps_data,stds,'.k','MarkerSize',0.001);
    elseif strcmp(mode,'nobars')
    else
        error('Unknown mode')
    end

    hline_1 = plot(ps_s1','b');
    hdots_1 = plot(ps_s1','.b','MarkerSize',6);
    
    hline_2 = plot(ps_s2','g');
    hdots_2 = plot(ps_s2','.g','MarkerSize',6);
    
    hline_3 = plot(ps_s4','r');
    hdots_3 = plot(ps_s4','.r','MarkerSize',6);
    
    xlim([0.5 8.5])
    ylim([0 0.7])
    grid on
    ylabel('Probability')
    
    legend([hbar, herror, hdots_1, hdots_2, hdots_3],'Brown et al. data','Statistical error','Stage 1 fit','Stage 2 fit','Stage 3 fit')
    
end

