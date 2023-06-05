function skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN, totalN, C )

global i n N W RA R G B classGroups classType CFLAG

    CFLAG = 0;

    R = zeros( N, totalN ); G = zeros( N, totalN ); B = zeros( N, totalN );
    
    RA = zeros( 2, 3, size( classGroups, 2)+1 );

    W = zeros( 2, 3, size( classGroups, 2) );

    for kk = 1:1:size( classGroups, 2 )
        for ii = 1:1:size(classType, 2 )
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e2;
            elseif ( ii == 2 )
    
                W(ii,:,kk) = 1e3;
            end
        end
    end

    for i = 1:trainingN                    

        n = skinObservation( i ); 
        
        RGB = dataSetRandomized( i, 1:C-1 );
        
        runningAverage( RGB );
    end
end