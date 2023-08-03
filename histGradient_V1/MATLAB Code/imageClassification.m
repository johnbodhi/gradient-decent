function [ D, E ] = imageClassification( dataSet, testObservation, verObservation )
    
    global RA Q pp

    D = 0; E = 0; 

    ii = 1;

    for i = 1:1:size(dataSet,1)

        % Store an RGB pixels of contained in the image of length 
        % imageLength to pass into the Gradient.

        rgbData( ii, 1:size(dataSet,2)) = dataSet( i, 1:size(dataSet,2)); ii = ii + 1; 
        
        if ( size( rgbData, 1 ) == size( dataSet, 2 ) )

            % rgbData = frameSieve(rgbData); % Duplicate and re-label each frame.

            % Observation_ = rgbData(:,C);

            % RA = Q; % We need to reset RA between classes...

            % We can utilize non-stationary RA during classification to
            % monitor dissimilarity between objects...

            % runningAverage( rgbData, Observation_ ); 

            % rgbData = kmeans( rgbData, Observation_ ); % k-means image data set.

            imgDecision = imageDecision( rgbData ); D = D + 1; % Take image to classify in the gradient.

            % Supervised Error...
            
%             if ( imgDecision ~= testObservation( i, 1 ) )
% 
%                 E = E + 1;
%             end
% 
%             % Display observation type, classifier decision, cumulative decision per
%             % class, and cumulative error per class...
% 
%             J = [ testObservation(i,1) imgDecision D E ]; disp( J )

            % Unsupervised Error...

            if ( imgDecision ~= verObservation( pp, 1 ) )

                E = E + 1;
            end
            pp = pp + 1;
            
            J = [ 0 imgDecision D E ]; disp( J )

            rgbData = 0; ii = 1;
        end
    end
end