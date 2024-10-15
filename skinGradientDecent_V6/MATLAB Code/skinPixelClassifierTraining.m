function skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN )

global RA W Q R G B Nl classGroups classType CFLAG

    CFLAG = 0;

    R = zeros( size(Nl,2), trainingN ); 
    
    %G = zeros( size(Nl,2), trainingN );
    
    %B = zeros( size(Nl,2), trainingN );

    % Allocate for cyclic weighting / infinite parameter gain.
    
    RA = zeros( size( classType, 2), 3, size( classGroups, 2) ); 

    W  = zeros( size( classType, 2), 3, size( classGroups, 2) );

    for kk = 1:1:size( classGroups, 2 )
        for ii = 1:1:size(classType, 2 )
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e1;
            elseif ( ii == 2 )
    
                W(ii,:,kk) = 1e2;

            elseif ( ii == 3 )
    
                W(ii,:,kk) = 0;
            end
        end
    end        
        
    RGB = dataSetRandomized;
    
    [ RA ] = runningAverage( RGB, [], skinObservation );

    Q = RA; % This is a redundant RA for resetting the average during classification.
end