function [ RA ] = filterCreation( A )

    RA = convGradDecent(A);

    % RA = krylovSubspace([],A);
end