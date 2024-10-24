function [ Y_ ] = addNoise( X_ )
     
    Y_ = X_;

    % We add AGWN as a random percentage of the individual time samples
    % contained in the series of frames.
    
    MU = 0; SD = 0.02;

    P = MU + SD * randn([1,1]);
    
    for i = 1:1:size(X_,1)

        P = MU + SD * randn([1,1]);
    
        Y_(i,1) = Y_(i,1) + floor(X_(i,1)*P);
    end