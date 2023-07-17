function [ D, E ] = skinImageClassification( dataSet, skinObservation, trainingN, testN )
    
    global i n A C imageLength imgDecision CFLAG

    CFLAG = 1;

    D = 0; E = 0; 

    ii = 1;

    for i = trainingN+1:1:trainingN+testN
        
        n = skinObservation( i, 1 );

        RGB = dataSet( i, 1:C);
        
        % Store an RGB pixels of contained in the image of length 
        % imageLength to pass into the Gradient.

        rgbData( ii, 1:C) = dataSet( i, 1:C); ii = ii + 1; 
        
        if ( size( rgbData, 1 ) == imageLength )

            rgbData = frameSieve(rgbData); % Duplicate and re-label each frame.

            skinObservation_ = rgbData(:,C);

            rgbData = kmeans( rgbData, skinObservation_ ); % k-means image data set.

            imgDecision = imageDecision( rgbData ); D = D + 1; % Take image to classify in the gradient.

%             if ( imgDecision ~= skinObservation_( 1, 1 ) )
% 
%                 E = E + 1;
%             end

            rgbData = 0; ii = 1;

            A = [ 0 imgDecision D E ]; disp( A )
        end

        runningAverage( RGB ); % We can utilize non-stationary RA during classification.
    end
end