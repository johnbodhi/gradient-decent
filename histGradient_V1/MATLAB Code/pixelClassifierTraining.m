function pixelClassifierTraining( dataSet )

    global RA Q 
    
    RA = zeros( 8, 10, 3 ); 
    
    [ RA ] = runningAverage( dataSet );

    Q = RA; % This is a redundant RA for resetting the average during classification.
end