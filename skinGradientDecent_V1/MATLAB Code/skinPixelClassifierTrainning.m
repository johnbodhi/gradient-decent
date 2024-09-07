function skinPixelClassifierTrainning( dataSetRandomized, skinObservation, trainingN, totalN, C )

global i n RA R G B CFLAG

    CFLAG = 0;

    N = 6;

    R = zeros( N, totalN ); G = zeros( N, totalN ); B = zeros( N, totalN ); 
    
    RA = zeros( N, C-1 );

    for i = 1:trainingN                    

        n = skinObservation( i ); 
        
        RGB = dataSetRandomized( i, 1:C-1 );
        
        runningAverage( RGB );
    end
end