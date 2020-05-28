function [ ] = plot_ps_with_data( ps, data, mode )
    
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

    hline = plot(ps','b');
    hdots = plot(ps','.b','MarkerSize',8);
    xlim([0.5 8.5])
    ylim([0 0.7])
    grid on
    ylabel('Probability')
    
    legend([hbar, herror, hdots],'Brown et al. data','Statistical error','Model fit')
    
end

