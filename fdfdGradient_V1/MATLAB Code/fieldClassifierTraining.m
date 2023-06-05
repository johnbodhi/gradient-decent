function fieldClassifierTraining( dataSetRandomized, fieldObservation, trainingN, totalN )

global classType i n C RA R CFLAG

    CFLAG = 0;
    
    N = size( classType, 2 );

    R = zeros( N, totalN );
    
    RA = zeros( N, C-1 );

    for i = 1:trainingN                    

        n = fieldObservation( i ); 
        
        S_ = dataSetRandomized( i, 1:C-1 );
        
        runningAverage( S_ );
    end
end