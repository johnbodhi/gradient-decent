function [ D, E ] = classifier( dataSet )
    
    global RA Q Train

    Train = 0;

    D = 0; E = 0;

    for i = 1:1:size(dataSet,1)

        % Store an RGB pixels of contained in the image of length 
        % imageLength to pass into the Gradient.

        histData(1,1:size(dataSet,2),:) = dataSet(i,1:size(dataSet,2),:);  

        histData = frameSieve(histData); % Duplicate and re-label each frame.

        Observation_ = histData(:,size(dataSet,2));

%         RA = Q; % We need to reset RA between classes...

        % We can utilize non-stationary RA during classification to
        % monitor dissimilarity between objects...

%         SVM( histData, [], Observation_ ); 

%         histData = kmeans( histData, Observation_ ); % k-means image data set.

        imgDecision = imageDecision( histData ); D = D + 1; % Take image to classify in the gradient.

        % Supervised Error...
        
        if ( imgDecision ~= dataSet(i,size(dataSet,2)) )

            E = E + 1;
        end

        % Display observation type, classifier decision, cumulative decision per
        % class, and cumulative error per class...

        J = [ dataSet(i,size(dataSet,2)) imgDecision D E ]; disp( J )

        % Unsupervised Error...

%             if ( imgDecision ~= Observation( pp, 1 ) )
% 
%                 E = E + 1;
%             end
%             pp = pp + 1;
%             
%             J = [ 0 imgDecision D E ]; disp( J )

        histData = zeros(1,size(dataSet,2),size(dataSet,3));
    end
end