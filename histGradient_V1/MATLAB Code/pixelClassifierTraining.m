function pixelClassifierTraining( dataSet, Sup, verObservation )

    global RA Q  
    
    [ RA ] = runningAverage( dataSet, Sup, verObservation );

    Q = RA; % Resetting RA during classification.
end