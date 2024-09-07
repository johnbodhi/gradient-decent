function [ D, E ] = classifier( dataSet, l )
    
    global RA Q Supervision

    D = 0; E = 0;

    for i = 1:1:size(dataSet,1)

        % Store an RGB pixels of contained in the image of length 
        % imageLength to pass into the Gradient.

        histData(1,1:size(dataSet,2)) = dataSet(i,1:size(dataSet,2));  

        histData = frameSieve(histData); % Duplicate and re-label each frame.

        l_ = histData(:,end,:);

        RA = Q; % We need to reset RA between classes...

        % We can utilize non-stationary RA during classification to
        % monitor dissimilarity between objects...

        histData = kmeans(histData, l_ ); % k-means image data set.

        frameD = frameDecision( histData ); D = D + 1; % Take image to classify in the gradient.
        
        if( Supervision )

            % Supervised Error...
            
            if ( frameD ~= dataSet(i,size(dataSet,2)) )
    
                E = E + 1;
            end
    
            % Display observation type, classifier decision, 
            % cumulative decision per
            % class, and cumulative error per class...
    
            J = [ dataSet(i,size(dataSet,2)) frameD D E ]; % disp( J )            
        else

            % Unsupervised Error...
    
            if ( frameD ~= l( i, 1 ) )
            
                E = E + 1;
            end

            % Display observation type, classifier decision, 
            % cumulative decision per
            % class, and cumulative error per class...
               
            J = [ l(i,1) frameD D E ]; disp( J )
        end

        histData = zeros(1,size(dataSet,2),size(dataSet,3));
    end
end