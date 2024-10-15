function skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN )

global i n N W RA R G B C classGroups classType CFLAG

    CFLAG = 0;

    R = zeros( N, trainingN ); 
    
    G = zeros( N, trainingN ); 
    
    B = zeros( N, trainingN );
    
    RA = zeros( size( classType, 2), 3, size( classGroups, 2) ); % Allocate for cyclic weighting / infinite parameter gain.

    W  = zeros( size( classType, 2), 3, size( classGroups, 2) );

    for kk = 1:1:size( classGroups, 2 )
        for ii = 1:1:size(classType, 2 )
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e1;
            elseif ( ii == 2 )
    
                W(ii,:,kk) = 1e3;

            elseif ( ii == 3 )
    
                W(ii,:,kk) = 0;
            end

        end
    end

    for i = 1:trainingN                    

        n = skinObservation( i ); 
        
        RGB = dataSetRandomized( i, 1:C-1 );
        
        [ RA ] = runningAverage( RGB );
    end
end