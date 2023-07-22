function pixelClassifierTraining( dataSetRandomized, Observation, trainingN )

global i n R G B ZA W C classGroups classType numRA CFLAG

    CFLAG = 0;

    R = zeros( 2, trainingN ); 
    
    G = zeros( 2, trainingN ); 
    
    B = zeros( 2, trainingN );

    % Allocate for cyclic weighting / infinite parameter gain.

    ZA = zeros( size( classType, 2), 3*numRA, size( classGroups, 2) ); 

    W  = zeros( size( classType, 2), 3*numRA, size( classGroups, 2) );

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

        n = Observation( i ); 
        
        RGB = dataSetRandomized( i, 1:C );
        
        [ ZA ] = runningAverage( RGB );
    end    

    % ZA = [ RA11 RA12 RA13 RA14 ]; These are block updates...
end