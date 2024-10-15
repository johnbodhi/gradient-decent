function pixelClassifierTraining( dataSetRandomized, Observation, trainingN )

global RA W Q R G B Nl classGroups classType CFLAG

    CFLAG = 0;

    R = zeros( size(Nl,2), trainingN ); 
    
%     G = zeros( size(Nl,2), trainingN );
%     
%     B = zeros( size(Nl,2), trainingN );

    % Allocate for cyclic weighting / infinite parameter gain.
    
    RA = zeros( size( classType, 2), 3, size( classGroups, 2) ); 

    W  = zeros( size( classType, 2), 3, size( classGroups, 2) );

    for kk = 1:1:size( classGroups, 2 )
        for ii = 1:1:size(classType, 2 )
        
            if ( ii == 1 )
    
                W(ii,:,kk) = 1e1;
            elseif ( ii == 2 )
    
                W(ii,:,kk) = 1e2;
            end
        end
    end        
        
    rgbData = dataSetRandomized;
    
    [ RA ] = runningAverage( rgbData, Observation );

    Q = RA; % This is a redundant RA for resetting the average during classification.
end