function [ Z ] = verifyFilter( X, Y )
    
    Z = krylovSubspace( X, Y );   
    
end