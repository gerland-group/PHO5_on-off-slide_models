function [ I_clustering ] = plot_parameter_value_analysis( values_con, values_reg, params_con, params_reg, C_text, model_infos_all, model_infos_subset, ...
    model_subset, normalization_factors, cluster_mode, I_manual, cluster_number, c_map, name_add )
% plots an overview of the used parameters/processes and their values

    N = length(model_subset);
    if ~exist('cluster_number','var') || isnan(cluster_number)
        cluster_number = round(N/2);
    end
    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
        name_add = name_add + "_not_rescaled";
    else
        normalization_factors = normalization_factors(model_subset);
    end
    
    num_params = length(C_text);

    values_con = values_con(model_subset,:);
    values_reg = values_reg(model_subset,:,:);
    params_con = params_con(model_subset,:);
    params_reg = params_reg(model_subset,:);  
    
    values_image = zeros(N,2*num_params);
    for m=1:N
        for p=1:size(params_con(m,params_con(m,:)>0),2)
            values_image(m,2*params_con(m,p)+[-1 0]) = values_con(m,p)*normalization_factors(m);
        end
        num_ind_params = size(params_reg(m,params_reg(m,:)>0),2);
        for p=1:num_ind_params
            values_image(m,2*params_reg(m,p)-1) = values_reg(m,p,2)*normalization_factors(m);   % regulated value: repressed
            values_image(m,2*params_reg(m,p)) = values_reg(m,p,1)*normalization_factors(m);   % regulated value: activated
        end
    end

    values_image_log = log10(values_image);
     
    indeces_not_inf = ~isinf(values_image_log);
    
    %caxis_lims = [-2.5 2.5];
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
    
    values_image_clustering = values_image_log;
    values_image_clustering(values_image_clustering ==-Inf) = -5;
    
    if strcmp('clustering on',cluster_mode)
        T = clusterdata(values_image_clustering,cluster_number);  % some standard clustering method (ToDo: look into the details)
        [T ,I_clustering] = sort(T);
    elseif strcmp('clustering manual',cluster_mode)
        I_clustering = I_manual;
    else
        I_clustering = 1:N;
    end
    
    if contains(pwd,"reduced_threshold") && N <= 10
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
            if size(model_infos_all,2) == 2
                ylabels_cell{i} = sprintf('%d (%1.2f)', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2));
            elseif size(model_infos_all,2) == 3
                ylabels_cell{i} = sprintf('%d\n(%1.2f, %1.2f)', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2),  model_infos_all(model_subset(I_clustering(i)),3));
            end 
        end
        format_yticks(gca,[],ylabels_cell,[],1:length(ylabels_cell),[],[],0)
        model_separation_line_width = 0.2;
    else
        model_separation_line_width = 0;
    end
    
    caxis([caxis_lims(1)-1.001*(caxis_lims(2)-caxis_lims(1))/(size(c_map,1)-1) caxis_lims(2)])  % to correct for the lowest colour being white, which shouldn't be on the colorbar and only used for real zero values
    
    h = colorbar;
    %h.TickLabels = {'10^{-3}','10^{-2}','10^{-1}','1','10','10^{2}','10^{3}'};
    %h.TickLabels = {'10^{-2}','10^{-1}','1','10','10^{2}'};
    if contains(name_add,"not_rescaled")
        h.Label.String = 'log_{10} relative rate value';
    else
        h.Label.String = 'log_{10} values (in 1/h)';
    end
    h.Label.FontSize = myfontsize*1.1;
    set(h, 'ylim', caxis_lims,'FontSize',myfontsize)
    
    if contains(pwd,"reduced_threshold") && N <= 10
        xticklabel_rotate(1.5:2:2*length(C_text),90,[-0.011 -0.055],C_text,'fontsize',myfontsize)
    else
        xticklabel_rotate(1.5:2:2*length(C_text),90,[-0.011 -0.025],C_text,'fontsize',myfontsize)
    end
    set(gca,'TickLength',[ 0 0 ])
    h.Position(3) = 0.015;
    
    % manual grid lines between data points ('grid on' gives the grid at the ticks, not between them)
    line_width = 0.5;
    hold on
    x = xlim;
    y = ylim;
    for xline = min(x):2:max(x)
        if any(xline==[ 2.5 + (0:10:30) 34.5 + (0:10:30) 66.5 68.5 + (0:6:18) ])
            if any(xline == [32.5 64.5])
                plot([xline; xline],y','Color',[0.4,0.4,0.4],'LineWidth',1)
            else
                plot([xline; xline],y','Color',[0.4,0.4,0.4],'LineWidth',line_width)
            end
        else
            plot([xline; xline],y','Color',[0.8,0.8,0.8],'LineWidth',line_width)
        end
    end
    if model_separation_line_width > 0
        for yline = min(y):1:max(y)
            plot(x',[yline;yline],'Color',[0.8,0.8,0.8],'LineWidth',model_separation_line_width)
        end
    end
    
    % seperation lines for regulated params
    for m=1:N_plot
        for p=params_reg(I_clustering(m),params_reg(I_clustering(m),:)>0)
            plot([2*p-0.5; 2*p-0.5],[m-0.5;m+0.5],'Color',[1,0.2,1],'LineWidth',1)
        end
    end
    h=xlabel('Processes');
    h.Position = h.Position - [0 0.01 0];
    h=ylabel('Models');
    if N_plot<=100
        h.Position = h.Position - [0.10 0 0];
    end
    hold off
    
    pause(0.5)
    print("img/parameter_value_analysis_" + name_add + "_" + length(I_clustering) + "_models",'-dpdf')
    close gcf
    s = convertStringsToChars("img/parameter_value_analysis_" + name_add + "_" + length(I_clustering) + "_models.pdf");
    open(s)
end
