function [] = plot_LH_histogram(LHs, LH_threshold, offset, x_max, my_title, name_add, my_xlab)
    
    figure('Position',[2000,100,320,160],'PaperUnits', 'centimeters','PaperSize', [8.5 4.25],'PaperPosition',[-0.05 0 9.3 4.25])
    histogram(-LHs+offset,0:0.25:x_max)
    hold on
    set(gca, 'FontSize',6)
    xlabel(my_xlab)
    ylabel('Number of models')
    xlim([0 x_max])
    y_lims = ylim;
    LH_best = -max(LHs)+offset;
    line([LH_best; LH_best],[0; y_lims(2)],'LineStyle','--','Color',[0.2 0.8 0.2])
    line([-LH_threshold+offset; -LH_threshold+offset],[0; y_lims(2)],'Color',[1 0.1 0.1])
    %text(20, 0.85*y_lims(2), my_title, 'fontsize', 7)
    title(my_title)
    
    pause(0.5)
    print("img/LH_histogram" + name_add,'-dpdf')
    close gcf
    s = convertStringsToChars("img/LH_histogram" + name_add + ".pdf");
    open(s)
end

