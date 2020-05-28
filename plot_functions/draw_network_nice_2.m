function  draw_network_nice_2( C, C_text, params_con, params_reg, fontsize, font_offset, rotation_mode )
% draws the network including constitutive and regulated parameters    

    shift = 0.25;
    linewidth = 0.5;
    maxheadsize = 0.002;
    
    params_con = params_con(params_con>0);
    params_reg = params_reg(params_reg>0);
    [params, I] = sort([params_con params_reg]);
    
    W_ind = calc_value_indeces_in_W(params, C, 8, 0);
    for i=1:max(max(W_ind))
        net_fluxes = W_ind==i;
        if i==1
            [x,y,u,v] = prepare_network_plot(net_fluxes, NaN, NaN, 'new', fontsize, font_offset, rotation_mode );
        else
            [x,y,u,v] = prepare_network_plot(net_fluxes, NaN, NaN, 'add', fontsize, font_offset, rotation_mode );
        end
        x = x + shift.*(-v./sqrt(u.*u+v.*v));  % shift arrows perpendicular to their direction to avoid arrow overlapping
        y = y + shift.*(u./sqrt(u.*u+v.*v));

        % This plots the same network outside the viewing area and prevents a bug where the arrowheads of single reaction parameters are too
        % small:
        x = [x; x];
        y = [y; y+1000];
        u = [u; u];
        v = [v; v];
        
        if ismember(params(i),params_reg)  % regulated parameter
            quiver(x,y,u,v,0,'LineWidth',2*linewidth,'MaxHeadSize',maxheadsize);
        else  % constitutive parameter
            quiver(x,y,u,v,0,'LineWidth',linewidth,'MaxHeadSize',maxheadsize);
        end
    end

    set(gca,'xtick',[],'ytick',[]);
    legend(C_text{params})
  
end
