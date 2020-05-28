function plot_net_flux_networks( fluxes, normalization_factors, ps, model_subset, I, I_clustering )

    N = length(model_subset);

    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
    else
        normalization_factors = normalization_factors(model_subset);
    end

    % find maximal net flux
    fmax = [];
    for data_set = 1:3
        for i=1:length(model_subset)
            net_fluxes_temp = calc_net_fluxes(fluxes{data_set}(:,:,model_subset(i))*normalization_factors(i));
            fmax = max([fmax max(max(net_fluxes_temp))]);
        end
    end

    defaultaxesfontsize = get(0,'defaultaxesfontsize');
    myfontsize = 6;
    set(0,'defaultaxesfontsize',myfontsize);
    
    if length(model_subset) <= 12
        rows = 4;
        cols = 3;
    elseif length(model_subset) <= 15
        rows = 5;
        cols = 3;
    elseif length(model_subset) <= 20
        rows = 5;
        cols = 4;
    else
        rows = 6;
        cols = 4;
    end
            
    for data_set = 1:3
        figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 26],'PaperPosition', [-2.5 -3 22.5 30.5])
        
        for i=1:length(model_subset)
            plot_pos = find(I_clustering==i);
            if isempty(plot_pos)
                continue
            end
            subplot(rows,cols,plot_pos);
            net_fluxes_temp = calc_net_fluxes(fluxes{data_set}(:,:,model_subset(i))*normalization_factors(i));
            fmax_temp = max(max(net_fluxes_temp));
            draw_flux_network(net_fluxes_temp, 0, fmax_temp, ps(model_subset(i),:,data_set), myfontsize-1, 0, 'horizontal')
            % draw_flux_network(net_fluxes_temp, 0, fmax_temp, ps(model_subset(i),:,data_set), myfontsize-1, 0.15, 'horizontal')
            title(sprintf('%d, max = %.3f / h',I(model_subset(i)), fmax_temp))
        end
    
        pause(0.5)
        print(sprintf('img/net_flux_networks_%d',data_set),'-dpdf')
        close gcf
        open(sprintf('img/net_flux_networks_%d.pdf',data_set))
    end
    set(0,'defaultaxesfontsize',defaultaxesfontsize);
end
