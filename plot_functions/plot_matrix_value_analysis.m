function [ max_ADS, min_ADS, I_clustering ] = plot_matrix_value_analysis( analysis_mode, array3d_cell, C, C_text, model_infos_all, model_infos_subset, ...
    model_subset, normalization_factors, cluster_mode, I_manual, cluster_number, c_map )
% activated and repressed state values in one plot
% array3d_cell can be fluxes or W_gains
    
    N = length(model_subset);
    if ~exist('cluster_number','var') || isnan(cluster_number)
        cluster_number = round(N/2);
    end
    
    if strcmp(analysis_mode,'all fluxes')
        caxis_lims = [-3 1];
        file_name = 'flux_analysis';
    elseif strcmp(analysis_mode,'net fluxes')
        caxis_lims = [-3 1];
        file_name = 'net_flux_analysis';
    elseif strcmp(analysis_mode,'reaction rates')
        caxis_lims = [-2 2];
        file_name = 'reaction_rate_analysis';
    else
        error('Unknown mode given.')
    end
    
    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
        scale_string = "without_time_scale";
    else
        normalization_factors = normalization_factors(model_subset);
        scale_string = "with_time_scale";
    end
    
    values_image = zeros(N,2*32);  % number of possible reactions for activated and repressed state
    n = 1;
    for p=1:length(C)
        if size(C{p},1)==1  % single reaction
            for m=1:N
                if strcmp(analysis_mode,'reaction rates') || strcmp(analysis_mode,'all fluxes')
                    values_image(m,n) = array3d_cell{2}(C{p}(1,1),C{p}(1,2),model_subset(m));
                    values_image(m,n+1) = array3d_cell{1}(C{p}(1,1),C{p}(1,2),model_subset(m));
                elseif strcmp(analysis_mode,'net fluxes')
                    values_image(m,n) = max([0 array3d_cell{2}(C{p}(1,1),C{p}(1,2),model_subset(m))-array3d_cell{1}(C{p}(1,2),C{p}(1,1),model_subset(m))]);
                    values_image(m,n+1) = max([0 array3d_cell{1}(C{p}(1,1),C{p}(1,2),model_subset(m))-array3d_cell{2}(C{p}(1,2),C{p}(1,1),model_subset(m))]);
                end
                x_label_text{(n+1)/2} = C_text{p};
                
                values_image(m,[n n+1]) = values_image(m,[n n+1]) * normalization_factors(m);
            end
            n = n+2;
        end
    end
    
    max_A = max(values_image(:,1:24),[],2);
    max_D = max(values_image(:,25:48),[],2);
    max_S = max(values_image(:,49:64),[],2);
    
    max_ADS = [max_A, max_D, max_S];
    
    min_A = min(values_image(:,1:24),[],2);
    min_D = min(values_image(:,25:48),[],2);
    min_S = min(values_image(:,49:64),[],2);
    
    min_ADS = [min_A, min_D, min_S];

    values_image_log = log10(values_image);
    
    indeces_not_inf = ~isinf(values_image_log);
    
    caxis_lims = [min(values_image_log(indeces_not_inf)) max(values_image_log(indeces_not_inf))];
    
    if any(any(values_image_log(indeces_not_inf)<caxis_lims(1)))
        warning('Some datapoints are increased to the lower color axis limit.')
        values_image_log(values_image_log<caxis_lims(1) & ~isinf(values_image_log))
        values_image_log(values_image_log<caxis_lims(1) & ~isinf(values_image_log)) = caxis_lims(1);
    end
    if any(any(values_image_log(indeces_not_inf)>caxis_lims(2)))
        warning('Some datapoints are decreased to the upper color axis limit.')
        values_image_log(values_image_log>caxis_lims(2) & ~isinf(values_image_log))
        values_image_log(values_image_log>caxis_lims(2) & ~isinf(values_image_log)) = caxis_lims(2);
    end
    
    if strcmp('clustering on',cluster_mode)
        T = clusterdata(values_image,cluster_number);  % some standard clustering method (ToDo: look into the details)
        [T ,I_clustering] = sort(T);
    elseif strcmp('clustering manual',cluster_mode)
        I_clustering = I_manual;
    else
        I_clustering = 1:N;
    end
    
    if contains(pwd,"reduced_threshold")
        h_factor = 0.45;
        h_offset = 0.3;
    else
        h_factor = 1;
        h_offset = 0;
    end
    
    figure('Position',[1700,200,765,h_factor*1000],'PaperUnits', 'centimeters','PaperSize', [17.7 h_factor*20],'PaperPosition', [-0.5 h_factor*-1.65+h_offset 20.3 h_factor*23.25-h_offset])

    myfontsize = 6;

    colormap(c_map)
    
    imagesc(values_image_log(I_clustering,:))
    max_ADS = max_ADS(I_clustering,:);
    min_ADS = min_ADS(I_clustering,:);
    
    set(gca,'FontSize',myfontsize)
    
    N_plot = length(I_clustering);
    
    if N_plot<=30
        for i=1:N_plot
            if size(model_infos_all,2) == 2
                ylabels_cell{i} = sprintf('%d\n%1.2f', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2));
            elseif size(model_infos_all,2) == 3
                ylabels_cell{i} = sprintf('%d\n%1.2f\n%1.2f', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2),  model_infos_all(model_subset(I_clustering(i)),3));
            end
        end
        format_yticks(gca,[],ylabels_cell,[],1:length(ylabels_cell),[],[],[])
        model_separation_line_width = 0.5;
    elseif N_plot<=100
        for i=1:N_plot
            ylabels_cell{i} = sprintf('%d (%1.2f)', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2));
        end
        format_yticks(gca,[],ylabels_cell,[],1:length(ylabels_cell),[],[],-0.01)
        model_separation_line_width = 0.2;
    else
        model_separation_line_width = 0;
    end
    
    caxis([caxis_lims(1)-1.001*(caxis_lims(2) - caxis_lims(1))/(size(c_map,1)-1) caxis_lims(2)])  % to correct for the lowest colour being white, which shouldn't be on the colorbar and only used for real zero values

    h = colorbar;
    if scale_string == "with_time_scale"
        h.Label.String = 'log_{10} values (in 1/h)';
    else
        h.Label.String = 'log_{10} relative rate value';
    end
    h.Label.FontSize = myfontsize*1.1;
    set(h, 'ylim', caxis_lims,'FontSize',myfontsize)
    
    if contains(pwd,"reduced_threshold")
        xticklabel_rotate(1.5:2:2*length(x_label_text),90,[-0.011 -0.055],x_label_text,'fontsize',myfontsize)
    else
        xticklabel_rotate(1.5:2:2*length(x_label_text),90,[-0.011 -0.025],x_label_text,'fontsize',myfontsize)
    end
    
    %xticklabel_rotate(1.5:2:2*length(x_label_text),90,[-0.011 -0.025],x_label_text,'fontsize',myfontsize)
    %xticklabel_rotate(1.5:2:2*length(x_label_text),90,[-0.011 -0.065],x_label_text,'fontsize',myfontsize)

    set(gca,'TickLength',[ 0 0 ])
    h.Position(3) = 0.015;
    
    % manual grid lines between data points ('grid on' gives the grid at the ticks, not between them)
    hold on
    x = xlim;
    y = ylim;
    dark_grey = [0 0 0];
    for xline = min(x):2:max(x)
        plot([xline; xline],y','Color',dark_grey,'LineWidth',0.5)
    end
    for xline = min(x)+1:2:max(x)
        plot([xline; xline],y','Color',[0.8,0.8,0.8],'LineWidth',0.5)
    end
    for xline = min(x)+8:8:8*8
        plot([xline; xline],y','Color',dark_grey,'LineWidth',1)
    end
    if model_separation_line_width > 0
        for yline = min(y):1:max(y)
            plot(x',[yline;yline],'Color',[0.8,0.8,0.8],'LineWidth',model_separation_line_width)
        end
    end
    
    h=xlabel('Reactions');
    h.Position = h.Position - [0 0.01 0];
    %h.Position = h.Position - [0 0.065 0];
    h=ylabel('Models');
    if N <= 30
        h.Position = h.Position - [0.10 0 0];
    end
    hold off
    
    pause(0.5)
    
    print(['img/' file_name sprintf('_%s',scale_string)],'-dpdf')
    close gcf
    open(['img/' file_name sprintf('_%s',scale_string) '.pdf'])
end
