function pixelClassifierTraining( dataSet )

    global RA Q 
    
    RA = zeros( 8, 10, 3 );  
    
    [ RA ] = runningAverage( dataSet );

    Q = RA; % Resetting RA during classification.
end