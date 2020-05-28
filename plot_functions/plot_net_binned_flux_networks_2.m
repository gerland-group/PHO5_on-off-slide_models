function plot_net_binned_flux_networks( fluxes, normalization_factors, C, model_subset, I, I_clustering )

    defaultaxesfontsize = get(0,'defaultaxesfontsize');
    myfontsize = 5;
    set(0,'defaultaxesfontsize',myfontsize);

    N = length(model_subset);
    
    if isnan(normalization_factors)
        normalization_factors = ones(N,1);
    else
        normalization_factors = normalization_factors(model_subset);
    end
    
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
    
    figure('Position',[1700,200,765,1000],'PaperUnits', 'centimeters','PaperSize', [18 26],'PaperPosition', [-2.5 -3 22.5 30.5])
    for m=1:length(model_subset)
        plot_pos = find(I_clustering==m);
        if isempty(plot_pos)
            continue
        end
            subplot(rows,cols,plot_pos);
        fluxes_mat = zeros(2,32);  % number of possible reactions
        for data_set = 1:2
            n = 1;
            for p=1:length(C)
                if size(C{p},1)==1  % single reaction
                    fluxes_mat(data_set,n) = fluxes{data_set}(C{p}(1,1),C{p}(1,2),model_subset(m))* normalization_factors(m);
                    n = n+1;
                end
            end
        end

        binned_fluxes_mat = zeros(2,10);
        for n=1:6
            binned_fluxes_mat(:,n) = sum(fluxes_mat(:,4*(n-1)+(1:4)),2);
        end
        for n=7:10
            binned_fluxes_mat(:,n) = sum(fluxes_mat(:,12+2*(n-1)+(1:2)),2);
        end

        % entries of binned_fluxes: A1  A2  A3  D1  D2  D3  S2>1 S2>3 S1>2 S3>2
        net_binned_fluxes_mat = binned_fluxes_mat - binned_fluxes_mat(:,[4:6 1:3 9 10 7 8]);
        fmax = max(max(abs(net_binned_fluxes_mat)));
        net_binned_fluxes_mat(net_binned_fluxes_mat<1e-2*fmax) = 0;

        draw_nucl_pov_fluxes_2(net_binned_fluxes_mat, myfontsize)
        title(sprintf('%d, max = %.3f / h',I(model_subset(m)), fmax))
    end

    pause(0.5)
    print(sprintf('img/net_summed_fluxes_switched'),'-dpdf')
    close gcf
    open(sprintf('img/net_summed_fluxes_switched.pdf'))

    set(0,'defaultaxesfontsize',defaultaxesfontsize);

end
