clc; clear all; close all; tic

% Classifications: Non-Cancerous = 1, Cancerous = 2, Unknown = 3

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global L C

classType = [ 1 2 3 ]; % Skin class designations.

Nr = 58; Mr = 58; % Photo length, and width. 

% Number of images per class to classify.

N = [ 4 4 4 ]; % Total photos.

M = 1; % Class training epochs 1-Total Photos. (Trains RA on a percentage of the pixels in M photos)

for j = 1:1:size( classType, 2 )

    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\Data\Excel Data");
    
    dataSet = readmatrix( 'trainRGB.csv' ); 
    
    L = size( dataSet, 1 ); C = size( dataSet, 2 );

    % photoToArray(); % Pre-process all images.

    % randomizedPhotos = randomizedPhotos( dataSet, Nr, Mr, L, C, N ); % Randomize all photos.
    
    % dataSetRandomized = dataSetRandomized( dataSet, L, C );  % Randomize all pixels.
    
    dataSetRandomized = readmatrix('dataSetRandomized_0.csv');
    
    dataSetRandomized = dataSetRandomized( 1:M * Nr * Mr, 1:C );  % Assign randomized pixels over the length of images.
    
    skinObservation = dataSetRandomized( :, C ); % Extract randmized observations for training.
    
    totalN = size( dataSetRandomized, 1 ); 
    
    trainingN = floor( 0.05 * totalN ); 
    
    testN = floor( 0.0 * totalN );
    
    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\MATLAB Code");
    
    skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN, totalN, C );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\Data\Excel Data");
    
    % We can feed the gradient with organized images, or randomized images,
    % and then select out which observations to feed.1143 Otherwise, we can
    % feed the gradient a totally random distribution of images.
    
    dataSet = readmatrix( 'testRGB.csv' ); % Organized.

    % We can grab all observations in order, no matter the order of the
    % data content.
    ii = 1;
    for i = 1:L
        if ( dataSet( i, C ) == classType( 1, j ) )
            dataSet_( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1;
        end
    end

    % dataSet = readmatrix( 'randomizedPhotos.csv' ); dataSet_ = dataSet; % Random assortment of images.
    clear dataSet
    % dataSet = dataSet_( 1:N * Nr * Mr, 1:C );
    
    dataSet = dataSet_( 1:N( 1, j ) * Nr * Mr, 1:C ); % We can iterate through all class images we desire.

    skinObservation = dataSet( :, C );
    
    totalN = size( dataSet, 1 );     

    trainingN = floor( 0.0 * totalN ); 
    
    testN = floor( 1.0 * totalN );
     
    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V3\MATLAB Code");
    
    [ D, E ] = skinImageClassification( dataSet, skinObservation, trainingN, testN, C );

    if( j < 3 )
    
        [ PREC( j ), REC( j ), ACC( j ), F1( j ) ] = fMeasure( D, E ); 
    end
    toc;
end

AVE = [ sum(PREC) sum(REC) sum(ACC) sum(F1) ] ./ (size(classType, 2) - 1);