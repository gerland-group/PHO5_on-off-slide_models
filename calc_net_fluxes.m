function [ net_fluxes ] = calc_net_fluxes( fluxes )

    net_fluxes = fluxes - permute(fluxes, [2 1 3]);
    net_fluxes(net_fluxes<0) = 0;

end

