function pixelClassifierTraining( dataSet, verObservation )

    global RA Q  
    
    [ RA ] = runningAverage( dataSet, verObservation );

    Q = RA; % Resetting RA during classification.
end