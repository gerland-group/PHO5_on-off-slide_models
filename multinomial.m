function coefficient = multinomial( ks )

    warning('off','MATLAB:nchoosek:LargeCoefficient');

    binomials = zeros(size(ks));
    for i = 1:length(ks)
        try  % to catch warning for to large coefficients
            binomials(i) = nchoosek(sum(ks(1:i)),ks(i));
        catch err
        end
    end
    coefficient =  prod(binomials);

end

