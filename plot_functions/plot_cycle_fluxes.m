function [ I_clustering ] = plot_cycle_fluxes( fluxes3d_cell, plot_mode, normalization_factors, C, C_text, model_subset, model_infos_all, model_infos_subset, cluster_mode, I_manual, cluster_number, c_map, data_text)
% activated and repressed state values in separate plots
    
    N = length(model_subset);
    
    if ~exist('cluster_number','var') || isnan(cluster_number)
        cluster_number = round(N/2);
    end
    
    if strcmp(plot_mode,'log')
    elseif strcmp(plot_mode,'linear')
    else
        error('Unknown plot mode given.')
    end
    
    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
    else
        normalization_factors = normalization_factors(model_subset);
    end
    
    caxis_log_lims = [-4 0];
    colorbar_string = 'cycle flux values (in 1/h)';

    for data_set=1:length(fluxes3d_cell)
    
        values_image = zeros(N,18);  % number of possible cycles including different directions
        x_label_text = {};

        for m=1:N
            net_fluxes = fluxes3d_cell{data_set}(:,:,model_subset(m)) - fluxes3d_cell{data_set}(:,:,model_subset(m))';
            
            values_image(m,[1 2]) = [1 -1] * net_fluxes(2,1);
            x_label_text{1} = '123';
            values_image(m,[3 4]) = [1 -1] * net_fluxes(1,4);
            x_label_text{2} = '134';
            values_image(m,[5 6]) = [1 -1] * net_fluxes(6,2);
            x_label_text{3} = '267';
            values_image(m,[7 8]) = [1 -1] * net_fluxes(4,6);
            x_label_text{4} = '456';
            
            values_image(m,[9 10]) = [1 -1] * (net_fluxes(2,1) + net_fluxes(2,3));
            x_label_text{5} = '273';
            values_image(m,[11 12]) = [1 -1] * (net_fluxes(1,4) + net_fluxes(3,4));
            x_label_text{6} = '354';
            values_image(m,[13 14]) = [1 -1] * (net_fluxes(5,6) + net_fluxes(4,6));
            x_label_text{7} = '586';
            values_image(m,[15 16]) = [1 -1] * (net_fluxes(6,7) + net_fluxes(6,2));
            x_label_text{8} = '687';

            values_image(m,[17 18]) = [1 -1] * (net_fluxes(3,5) + net_fluxes(3,4) + net_fluxes(1,4));
            x_label_text{9} = '3785';
            
            %test = net_fluxes(3,5) + net_fluxes(3,4) + net_fluxes(1,4) - (net_fluxes(7,3) + net_fluxes(2,3) + net_fluxes(2,1))

            values_image(m,:) = values_image(m,:) * normalization_factors(m);
        end
        values_image(values_image<0) = 0;


        file_name = 'cycle_flux_analysis';
        if strcmp('clustering on',cluster_mode)
            T = clusterdata(values_image,cluster_number);  % some standard clustering method (ToDo: look into the details)
            [T ,I_clustering] = sort(T);
            file_name = 'cycle_flux_analysis_clustered';
        elseif strcmp('clustering manual',cluster_mode)
            I_clustering = I_manual;
        else
            I_clustering = 1:N;
        end  
        
        if strcmp(plot_mode,'log')
            values_image = log10(values_image);

            indeces_not_inf = ~isinf(values_image);
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
        else
            values_image(values_image==0) = -Inf;
            indeces_not_inf = ~isinf(values_image);
            if any(any(values_image(indeces_not_inf)<10^caxis_log_lims(1)))
                warning('Some datapoints are increased to the lower color axis limit.')
                values_image(values_image<10^caxis_log_lims(1) & ~isinf(values_image))
                values_image(values_image<10^caxis_log_lims(1) & ~isinf(values_image)) = 10^caxis_log_lims(1);  % do not replace -inf values
            end
            if any(any(values_image(indeces_not_inf)>10^caxis_log_lims(2)))
                warning('Some datapoints are decreased to the upper color axis limit.')
                values_image(values_image>10^caxis_log_lims(2) & ~isinf(values_image))
                values_image(values_image>10^caxis_log_lims(2) & ~isinf(values_image)) = 10^caxis_log_lims(2);  % do not replace -inf values
            end
        end
        
        if contains(pwd,"reduced_threshold")
            h_factor = 0.45;
            h_offset = 0.55;
        else
            h_factor = 1;
            h_offset = 0;
        end

        figure('Position',[1700,200,765,h_factor*1000],'PaperUnits', 'centimeters','PaperSize', [17.7 h_factor*20],'PaperPosition', [-0.5 h_factor*-1.65+h_offset 20.3 h_factor*23.25-h_offset-0.2])

        %figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [17.7 10],'PaperPosition', [0.5 -0.3 19.5 10.2])
        myfontsize = 6;

        colormap(c_map)

        imagesc(values_image(I_clustering,:))
        set(gca,'FontSize',myfontsize)

        if N <= 30
            for i=1:N
                %ylabels_cell{i} = sprintf('%d\n%1.2f', model_infos_all(model_subset(I_clustering(i)),1), model_infos_all(model_subset(I_clustering(i)),2));
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
            caxis([10^caxis_log_lims(1)-1.001*(10^caxis_log_lims(2) - 10^caxis_log_lims(1))/(size(c_map,1)-1) 10^caxis_log_lims(2)])  % to correct for the lowest colour being white, which shouldn't be on the colorbar and only used for real zero values
            set(h, 'ylim', 10.^caxis_log_lims,'FontSize',myfontsize)
        end
        
        %xticks(1.5:2:18)
        %xticklabels(x_label_text)
        xticklabel_rotate(1.5:2:2*length(x_label_text),0,[-0.015 0],x_label_text,'fontsize',myfontsize) 

        set(gca,'TickLength',[ 0 0 ])
        h.Position(3) = 0.015;

        % manual grid lines between data points ('grid on' gives the grid at the ticks, not between them)
        hold on
        x = xlim;
        y = ylim;
        dark_grey = [0 0 0];
        for xline = min(x):1:max(x)
            plot([xline; xline],y','Color',[0.8 0.8 0.8],'LineWidth',0.5)
        end
        for xline = min(x):2:max(x)
            plot([xline; xline],y','Color',dark_grey,'LineWidth',1)
        end
        for yline = min(y):1:max(y)
            plot(x',[yline;yline],'Color',[0.8 0.8 0.8],'LineWidth',model_separation_line_width)
        end
    
        h=xlabel('Cycles (using state numbers)');
        h.Position = h.Position - [0 0.01 0];
        %h.Position = h.Position - [0 0.04 0];

        h=title(data_text{data_set});
        h.Position = h.Position + [0 0.025 0];
        
        h=ylabel('Models');
        if N <= 30
            h.Position = h.Position - [0.10 0 0];
        end
        hold off
        
        pause(0.5)

        print(['img/' file_name '_' plot_mode '_' sprintf('%d',data_set)],'-dpdf')

        close gcf
        open(['img/' file_name '_' plot_mode '_' sprintf('%d',data_set) '.pdf'])
    end
end
