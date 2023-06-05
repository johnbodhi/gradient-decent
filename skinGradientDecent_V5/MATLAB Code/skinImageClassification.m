function [ D, E ] = skinImageClassification( dataSet, skinObservation, trainingN, testN, C )
    
    global i n imageLength imgDecision CFLAG

    CFLAG = 1;

    D = 0; E = 0; 

    ii = 1;

    Nr = 25; Mr = 25;

    imageLength = Nr * Mr;

    for i = trainingN+1:trainingN+testN
        
        n = skinObservation( i );

        RGB = dataSet( i, 1:C ); 
        
        rgbData( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1; % Store an RGB pixels of contained in the image of length imageLength to pass into Gradient Dencent.
        
        if ( size( rgbData, 1 ) == imageLength )

            rgbData = kmeans( rgbData, skinObservation ); % k-means image data set.

            imgDecision = imageDecision( rgbData ); D = D + 1; % Take image to classify.

            if ( imgDecision ~= skinObservation( i ) )
                E = E + 1;
            end

            rgbData = 0; ii = 1;
            A = [ n imgDecision D E ];
            disp( A )
        end

        % runningAverage( RGB ); % We can utilize non-stationary RA during classification.
    end

end