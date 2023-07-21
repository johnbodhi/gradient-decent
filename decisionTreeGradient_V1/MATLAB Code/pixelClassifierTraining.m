function pixelClassifierTraining( dataSetRandomized, Observation, trainingN )

global i n W RA R G B C M classGroups classType CFLAG ZA RA11 RA12 RA13 RA14

    CFLAG = 0;

    R = zeros( N, trainingN ); 
    
    G = zeros( N, trainingN ); 
    
    B = zeros( N, trainingN );

    % Allocate for cyclic weighting / infinite parameter gain.

    ZA = zeros( size( classType, 2), size(ZA,2), size( classGroups, 2) ); 

    W  = zeros( size( classType, 2), size(ZA,2), size( classGroups, 2) );

%     RA = zeros( size( classType, 2), 3, size( classGroups, 2) );  
%     W  = zeros( size( classType, 2), 3, size( classGroups, 2) );

%     RA11 = zeros( size( classType, 2), 3, size( classGroups, 2) ); 
% 
%     RA12 = zeros( size( classType, 2), 3, size( classGroups, 2) ); 
% 
%     RA13 = zeros( size( classType, 2), 3, size( classGroups, 2) ); 
% 
%     RA14 = zeros( size( classType, 2), 3, size( classGroups, 2) ); 
%
%     W    = zeros( size( classType, 2), 3, size( classGroups, 2) );

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