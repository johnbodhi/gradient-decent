clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global L C X N Nl RA GFLAG classGroups classType imageLength uu vv pp

classType = [ 1 2 ]; % Number of column-wise designations.

classGroups = zeros( 1, 0.50 * size( numImages, 2 ) ); % Groupings for cyclic weight. Over allocate at classGroups = 1

N = size(numImages,2);

Nr = 25; Mr = 25; imageLength = Nr * Mr; % Photo length, and width. 

% Number of images per class to classify.

Nl = zeros(size(numImages,2),1);
for i = 1:1:size(Nl,1)

    Nl(i,1) = numImages(i); % Number of objects per class.
    %Nl(i,1) = 5; 
end
totalN = sum(Nl); 

% We can generate an objective label vector to keep track of our errors
% with unsupervised data...

ii = 1; jj = 1; kk = 1; pp = 1;
for i = 1:Nl(ii,1):2*totalN

    while ( jj <= Nl(ii,1) )

        verObservation(kk,1) = ii;

        jj = jj + 1; kk = kk + 1;
    end

    ii = ii + 1; jj = 1;

    if ( size(verObservation,1) >= totalN )

        break;
    end
end

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\Data\Excel Data");

dataSet = readmatrix( 'trainRGB.csv' ); % Samples not in test.

L = size( dataSet, 1 ); C = size( dataSet, 2 );

% dataSetRandomized = dataSetRandomized( dataSet, L, C );  % Randomize all pixels.

dataSetRandomized = readmatrix( 'dataSetRandomized.csv' ); % Only samples in train.

M = 1; % Class training epochs 1-total photos. (Trains RA on a percentage of the pixels in M photos)

dataSetRandomized = dataSetRandomized( 1:M * Nr * Mr, 1:C );  % Assign randomized pixels over the length of images.

skinObservation = dataSetRandomized( :, C ); % Extract randmized observations for training.

totalN = size( dataSetRandomized, 1 ); 

trainingN = floor( 0.30 * totalN ); 

testN = floor( 0.0 * totalN );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\MATLAB Code");

skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We can feed the gradient with organized images, or randomized images,
% and then select out which observations to feed. Otherwise, we can
% feed the gradient a totally random distribution of images.

% In this version the classifier is fed in an indefinite number of batches.

% Supervised Verification / Test...

uu = 1; vv = 2; hh = 1; cc = 1; T = 0;

for k = 1:1:size(RA,3)

    if ( hh > size(Nl,1) )

        break; % End test sequence.
    end

    X = [ uu, vv ]; % Batch indexes.

    GFLAG = 1;

    for j = 1:1:size(X,2)        

        cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\Data\Excel Data");

        % Test sequences not included in training data...

        % dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.
    
        % dataSet = readmatrix( 'testRGB (1).csv' ); % Supervised test sequence.

        dataSet = readmatrix( 'testRGB (2).csv' ); % Unsupervised test sequence.        

        % dataSet = randomizedPhotos( dataSet, Nr, Mr, L, C, Nl ); % Randomize all photos.

        L = size(dataSet,1);
    
        % We can grab all observations in order, no matter the order of the
        % data content. The scoop.

        % Supervised verification scoop.

%         ii = 1;      
%         for i = 1:1:L % We can choose more than one photo per class.
%             if ( dataSet( i, C ) == X(1,j) )
% 
%                 dataSet_( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1;
%             end
%         end    
        
        ii = 1;      
        for i = ( 1 + T ):1:( Nl(cc,1)*imageLength + T ) % We can choose more than one photo per class.
            if ( dataSet( i, C ) == 0 )

                dataSet_( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1;
            end
        end    
        T = T + Nl(cc,1)*imageLength; cc = cc + 1;
    
        % dataSet = readmatrix( 'randomizedPhotos.csv' ); % Random assortment of images.
        % dataSet_ = dataSet; 
        clear dataSet
        % dataSet = dataSet_( 1:N * Nr * Mr, 1:C );

        % We can iterate through all class images we desire.
        
        % dataSet = dataSet_( 1:Nl( 1, hh ) * Nr * Mr, 1:C ); 

        dataSet = dataSet_( :, 1:C ); 
        
        clear dataSet_

        cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V6\MATLAB Code");
        
        skinObservation = dataSet( :, C );
        
        totalN = size( dataSet, 1 );     
    
        trainingN = floor( 0.0 * totalN ); 
        
        testN = floor( 1.0 * totalN );          
        
        [ D, E ] = skinImageClassification( dataSet, testN, verObservation );

        if( hh <= size( Nl, 1 ) )
        
            [ PREC( hh ), REC( hh ), ACC( hh ), F1( hh ) ] = fMeasure( D, E ); 
        end

        hh = hh + 1; GFLAG = 0;
    end
    
    uu = uu + 2; vv = vv + 2; 
end

AVE = [ mean(PREC) mean(REC) mean(ACC) mean(F1) ];

toc;