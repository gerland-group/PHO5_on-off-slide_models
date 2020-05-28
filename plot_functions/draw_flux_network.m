function [  ] = draw_flux_network( fluxes, shift, scaling, ps, fontsize, font_offset, rotation_mode )
    
    linewidth = 0.5;
    maxheadsize = 0.3;

    [x,y,u,v] = prepare_network_plot(fluxes, scaling, ps, 'new', fontsize, font_offset, rotation_mode, 1);
    x = x + shift.*(-v./sqrt(u.*u+v.*v));  % shift arrows perpendicular to their direction to avoid arrow overlapping
    y = y + shift.*(u./sqrt(u.*u+v.*v));
    quiver(x,y,u,v,0,'LineWidth',linewidth,'MaxHeadSize',maxheadsize)
    set(gca,'xtick',[],'ytick',[])
    
end
