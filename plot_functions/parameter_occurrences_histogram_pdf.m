function [  ] = parameter_occurrences_histogram_pdf( params_con, params_reg, C_text, filename, my_title, with_legend, my_ylab)

    params_con_vec = reshape(params_con,[],1);
    params_reg_vec = reshape(params_reg,[],1);
    
    figure('Position',[2000,100,360,180],'PaperUnits', 'centimeters','PaperSize', [8.5 4.25],'PaperPosition',[0 0.25 9.5 4])
    histogram(params_con_vec); hold on
    set(gca, 'FontSize',6)
    histogram(params_reg_vec)
    xlim([0.5 length(C_text)+0.5])
    ylim([0, size(params_con,1)])
    %xticklabel_rotate(1:length(C_text),90,[-0.014 -0.24],C_text,'fontsize',4.5)
    xticklabel_rotate(1:length(C_text),90,[-0.014 -0.10],C_text,'fontsize',4.5)
    title(my_title);
    if with_legend
        legend('Constitutive processes','Regulated processes','Location','northeast')
    end
    h=ylabel(sprintf('Number of occurrences\n%s', my_ylab));
    %h.Position = [-0.122 0.5 0];
    h=xlabel('Process');
    %h.Position = h.Position + [0 -0.04 0];
    pause(0.5)
    print(['img/' filename],'-dpdf')
    close gcf
    open(['img/' filename '.pdf'])

end
