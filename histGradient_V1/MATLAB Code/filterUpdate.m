function [ Z ] = filterUpdate( X, Y )

    Z = krylovSubspace( X, Y ); 
end

