function [ D, E ] = fieldClassification( dataSet, fieldObservation, trainingN, testN, C )
    
    global i n vectorLength targetDecision CFLAG

    CFLAG = 1;

    D = 0; E = 0; 

    ii = 1;

    for i = trainingN+1:trainingN+testN

        n = fieldObservation( i ); 

        S_ = dataSet( i, 1:C ); 
        
        S( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1; % Store an RGB pixels of contained in the image of length imageLength to pass into Gradient Dencent.
        
        if ( size( S, 1 ) == vectorLength )

            S = kmeans( S, fieldObservation ); % k-means image data set.

            targetDecision = imageDecision( S ); D = D + 1; % Take image to classify.

                if ( targetDecision ~= fieldObservation( i ) )
                    E = E + 1;
                end

            S = 0; ii = 1;
            A = [ i n targetDecision D E ];
            disp( A )
        end
        % runningAverage( S_ );
    end
end