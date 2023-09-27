function [ Z ] = krylovSubspace( X, Y )
    
    Z = BiCGSTAB( X, Y );
    
end