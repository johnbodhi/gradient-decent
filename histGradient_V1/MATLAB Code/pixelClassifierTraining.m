function pixelClassifierTraining( dataSet )

global RA Q CFLAG

    CFLAG = 0;

    % Allocate for cyclic weighting / infinite parameter gain.
    
    RA = zeros( size( classType, 2), 30, size( classGroups, 2) ); 
    
    [ RA ] = runningAverage( dataSet );

    Q = RA; % This is a redundant RA for resetting the average during classification.
end