function [ RA_NN, RA_GMR ] = filterCreation( A_ )
    
    [ RA_NN ]  = CNN( A_ );
    
    [ RA_GMR ] = GMRES( A_ );
end