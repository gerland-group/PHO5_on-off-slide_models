function [ I_clustering ] = plot_matrix_value_analysis_2( analysis_mode, plot_mode, array3d_cell, normalization_factors, C, C_text, model_subset, model_infos_all, model_infos_subset, cluster_mode, I_manual, cluster_number, c_map, data_text)
% activated and repressed state values in separate plots
% array3d_cell can be fluxes or W_gains
    
    N = length(model_subset);
    
    if ~exist('cluster_number','var') || isnan(cluster_number)
        cluster_number = round(N/2);
    end
    
    if strcmp(analysis_mode,'all fluxes')
        caxis_log_lims = [-4.5 1.5];
        caxis_linear_lims = [0 6];
        file_name = 'flux_analysis';
        colorbar_string = 'flux values (in 1/h)';
        binning = 1;
    elseif strcmp(analysis_mode,'net fluxes')
        caxis_log_lims = [-3 1];
        caxis_linear_lims = [0 2];
        file_name = 'net_flux_analysis';
        colorbar_string = 'net flux values (in 1/h)';
        binning = 1;
    elseif strcmp(analysis_mode,'reaction rates')
        caxis_log_lims = [-3 3];
        caxis_linear_lims = [0 100];
        file_name = 'reaction_rate_analysis';
        colorbar_string = 'reaction rate values (in 1/h)';
        binning = 0;
    else
        error('Unknown analysis mode given.')
    end
    
    if strcmp(plot_mode,'log')
    elseif strcmp(plot_mode,'linear')
    elseif strcmp(plot_mode,'linear_cutoff')
    else
        error('Unknown plot mode given.')
    end
    
    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
        scale_string = "without_time_scale";
    else
        normalization_factors = normalization_factors(model_subset);
        scale_string = "with_time_scale";
    end

    for data_set=1:length(array3d_cell)
    
        values_image = zeros(N,32);  % number of possible reactions
        x_label_text = {};
        n = 1;
        for p=1:length(C)
            if size(C{p},1)==1  % single reaction
                for m=1:N
                    if strcmp(analysis_mode,'reaction rates') || strcmp(analysis_mode,'all fluxes')
                        values_image(m,n) = array3d_cell{data_set}(C{p}(1,1),C{p}(1,2),model_subset(m));
                    elseif strcmp(analysis_mode,'net fluxes')
                        values_image(m,n) = max([0 array3d_cell{data_set}(C{p}(1,1),C{p}(1,2),model_subset(m))-array3d_cell{data_set}(C{p}(1,2),C{p}(1,1),model_subset(m))]);
                    end
                    x_label_text{n} = C_text{p};
                    values_image(m,n) = values_image(m,n) * normalization_factors(m);
                end
                n = n+1;
            end
        end
        
        if strcmp('clustering on',cluster_mode)
            T = clusterdata(values_image,cluster_number);  % some standard clustering method (ToDo: look into the details)
            [T ,I_clustering] = sort(T);
        elseif strcmp('clustering manual',cluster_mode)
            I_clustering = I_manual;
        else
            I_clustering = 1:N;
        end
        
        values_image_binned = zeros(N,10);
        for n=1:6
            values_image_binned(:,n) = sum(values_image(:,4*(n-1)+(1:4)),2);
        end
        for n=7:10
            values_image_binned(:,n) = sum(values_image(:,12+2*(n-1)+(1:2)),2);
        end
        
        if binning == 1
            values_image = [values_image zeros(N,1) values_image_binned];
        end
        
        if strcmp(plot_mode,'log')
            values_image = log10(values_image);
        
            indeces_not_inf = ~isinf(values_image);
             
            caxis_log_lims = [min(values_image(indeces_not_inf)) max(values_image(indeces_not_inf))];
            
            if any(any(values_image(indeces_not_inf)<caxis_log_lims(1)))
                warning('Some datapoints are increased to the lower color axis limit.')
                values_image(values_image<caxis_log_lims(1) & ~isinf(values_image))
                values_image(values_image<caxis_log_lims(1) & ~isinf(values_image)) = caxis_log_lims(1);  % do not replace -inf values
            end
            if any(any(values_image(indeces_not_inf)>caxis_log_lims(2)))
                warning('Some datapoints are decreased to the upper color axis limit.')
                values_image(values_image>caxis_log_lims(2) & ~isinf(values_image))
                values_image(values_image>caxis_log_lims(2) & ~isinf(values_image)) = caxis_log_lims(2);  % do not replace -inf values
            end
        elseif strcmp(plot_mode,'linear')
            caxis_linear_lims = [min(min(values_image(values_image>0))) max(max(values_image))];
            
            %values_image(values_image<caxis_linear_lims(1)) = 0;
            
%             if any(any(values_image<caxis_linear_lims(1) & values_image>0))
%                 warning('Some datapoints are increased to the lower color axis limit.')
%                 values_image(values_image<caxis_linear_lims(1) & values_image>0)
%                 values_image(values_image<caxis_linear_lims(1) & values_image>0) = caxis_linear_lims(1);  % do not replace zero values
%             end
%             if any(any(values_image>caxis_linear_lims(2)))
%                 warning('Some datapoints are decreased to the upper color axis limit.')
%                 values_image(values_image>caxis_linear_lims(2))
%                 values_image(values_image>caxis_linear_lims(2)) = caxis_linear_lims(2);
%             end
        else  % linear_cutoff
            model_maxima = max(values_image(:,1:32), [], 2);
            for m = 1:length(model_maxima)
                below_cutoff = values_image(m,:) < 0.1 * model_maxima(m);
                values_image(m,below_cutoff) = 0;
            end
            caxis_linear_lims = [min(min(values_image(values_image>0))) max(max(values_image))];
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

        imagesc(values_image(I_clustering,:))
        set(gca,'FontSize',myfontsize)

        if N<=30
            for i=1:N
                %ylabels_cell{i} = sprintf('%d\ndLH %1.2f', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2));
                ylabels_cell{i} = sprintf('%d', model_infos_all(model_subset(I_clustering(i)),1));
            end
            format_yticks(gca,[],ylabels_cell,[],1:length(ylabels_cell),[],[],[])
            model_separation_line_width = 0.5;
        else
            model_separation_line_width = 0.2;
        end

        h = colorbar;
        h.Label.FontSize = myfontsize*1.1;

        if strcmp(plot_mode,'log')
            h.Label.String = ['log_{10} ' colorbar_string];
            caxis([caxis_log_lims(1)-1.001*(caxis_log_lims(2) - caxis_log_lims(1))/(size(c_map,1)-1) caxis_log_lims(2)])  % to correct for the lowest colour being white, which shouldn't be on the colorbar and only used for real zero values
            set(h, 'ylim', caxis_log_lims,'FontSize',myfontsize)
        else
            h.Label.String = colorbar_string;
            caxis([caxis_linear_lims(1)-1.001*(caxis_linear_lims(2) - caxis_linear_lims(1))/(size(c_map,1)-1) caxis_linear_lims(2)])  % to correct for the lowest colour being white, which shouldn't be on the colorbar and only used for real zero values
            set(h, 'ylim', caxis_linear_lims,'FontSize',myfontsize)
        end
        
        if contains(pwd,"reduced_threshold")
            x_label_offset = -0.055;
        else
            x_label_offset = -0.025;
        end
        
        if binning == 1
            x_label_text = [x_label_text {'         ', C_text{2}, C_text{7}, C_text{12}, C_text{18}, C_text{23}, C_text{28}, C_text{35}, C_text{38}, C_text{41}, C_text{44}}];
            xticklabel_rotate(0.5:1:length(x_label_text),90,[0.002 x_label_offset],x_label_text,'fontsize',myfontsize)
        else
            xticklabel_rotate(0.5:1:length(x_label_text),90,[0.006 x_label_offset],x_label_text,'fontsize',myfontsize)
        end

        set(gca,'TickLength',[ 0 0 ])
        h.Position(3) = 0.015;

        % manual grid lines between data points ('grid on' gives the grid at the ticks, not between them)
        hold on
        x = xlim;
        y = ylim;
        dark_grey = [0 0 0];
        for xline = min(x):1:max(x)
            plot([xline; xline],y','Color',[0.8,0.8,0.8],'LineWidth',0.5)
        end
        for yline = min(y):1:max(y)
            plot(x',[yline;yline],'Color',[0.8,0.8,0.8],'LineWidth',model_separation_line_width)
            plot(min(x) + [32; 33],[yline;yline],'Color',[1 1 1],'LineWidth',1)
        end
        for xline = min(x)+[0:4:8*4 33 36 39 41 43] 
            plot([xline; xline],y','Color',dark_grey,'LineWidth',1)
        end
        
        h=xlabel('Reactions');
        h.Position = h.Position - [0 0.01 0];
        h=ylabel(['Models (' data_text{data_set} ')']);
        if N <= 30
            h.Position = h.Position - [0.10 0 0];
        end
        hold off
        
        pause(0.5)

        print(['img/' file_name '_' plot_mode sprintf('_%s_%d',scale_string,data_set)],'-dpdf')
        close gcf
        open(['img/' file_name '_' plot_mode sprintf('_%s_%d',scale_string,data_set) '.pdf'])
    end
end
