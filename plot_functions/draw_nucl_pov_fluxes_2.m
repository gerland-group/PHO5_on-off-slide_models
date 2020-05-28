function  draw_nucl_pov_fluxes_2( net_binned_fluxes_mat, myfontsize )
% net summed fluxes from the nucleosome point of view (plot can be missleading, always compare with actual network flux plot)

    maxarrowwith = 1;
   
    % Act
    % 3 2 1     A1  A2  A3  D1  D2  D3  S2>1 S2>3 S1>2 S3>2
    x =        [7.5 4.5 1.5 7.5 4.5 1.5 5.5  3.5  6.5  2.5];
    y =        [3   3   3   2   2   2   1    1    1    1  ];
    alpha = pi*[3/2 3/2 3/2 1/2 1/2 1/2 0    1    1    0  ];
    
    x = [x x];
    y = [y(1:6) y(7:10)+0.75 -1 -1 -1 0 0 0 y(7:10)-0.75];
    alpha = [alpha alpha(1:6)-pi alpha(7:10)];
    
    % switch repressed and activated arrow positions
    y = [y(11:20) y(1:10)];
    alpha = [alpha(11:20) alpha(1:10)];
    
    max_flux = max(max(net_binned_fluxes_mat));
    
    rectangle('Position',[1 0.5 1 1],'LineStyle','-','Curvature',[1 1])
    rectangle('Position',[4 0.5 1 1],'LineStyle','-','Curvature',[1 1])
    rectangle('Position',[7 0.5 1 1],'LineStyle','-','Curvature',[1 1])
    
    text(1.2,1,'N-3','FontSize',myfontsize)
    text(4.2,1,'N-2','FontSize',myfontsize)
    text(7.2,1,'N-1','FontSize',myfontsize)
    
    net_binned_fluxes_vec = [net_binned_fluxes_mat(1,:) net_binned_fluxes_mat(2,:)];
    for i=1:length(net_binned_fluxes_vec)
        if i<=10
            mycolor = [0.2 0.8 0.1];
        else
            mycolor = [0.8 0.2 0.1];
        end
        if net_binned_fluxes_vec(i)>0
            draw_arrow(x(i),y(i),alpha(i),net_binned_fluxes_vec(i)/max_flux*maxarrowwith, mycolor)
        end
    end
    set(gca,'xtick',[],'ytick',[]);
    xlim([0.5 8.5])
    ylim([-1.5 3.5])
    axis equal
    
end


