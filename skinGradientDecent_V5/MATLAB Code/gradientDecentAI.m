% clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global L C X N Nl RA classGroups classType uu vv

classType = [ 1 2 ]; % Number of column-wise designations.

classGroups = zeros( 1, size(numImages, 2)); % Groupings for cyclic weight.

classGroups(1,1:end) = size(classType,2); 

Nr = 25; Mr = 25; % Photo length, and width. 

% Number of images per class to classify.

N = size(classType,2)*size(classGroups,2); % Total number of classes.

Nl = zeros(1,size(numImages,2));
for j = 1:1:size(numImages,2)

    % Nl(1,j) = numImages(j); % Number of objects per class.
    Nl(1,j) = 1; 
end
totalN = sum(Nl); 

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Excel Data");

dataSet = readmatrix( 'trainRGB.csv' ); 

L = size( dataSet, 1 ); C = size( dataSet, 2 );

% randomizedPhotos = randomizedPhotos( dataSet, Nr, Mr, L, C, Nl ); % Randomize all photos.

% dataSetRandomized = dataSetRandomized( dataSet, L, C );  % Randomize all pixels.

dataSetRandomized = readmatrix('dataSetRandomized.csv');

M = 1; % Class training epochs 1-Total Photos. (Trains RA on a percentage of the pixels in M photos)

dataSetRandomized = dataSetRandomized( 1:M * Nr * Mr, 1:C );  % Assign randomized pixels over the length of images.

skinObservation = dataSetRandomized( :, C ); % Extract randmized observations for training.

totalN = size( dataSetRandomized, 1 ); 

trainingN = floor( 0.3 * totalN ); 

testN = floor( 0.0 * totalN );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\MATLAB Code");

skinPixelClassifierTraining( dataSetRandomized, skinObservation, trainingN );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We can feed the gradient with organized images, or randomized images,
% and then select out which observations to feed. Otherwise, we can
% feed the gradient a totally random distribution of images.

uu = 1; vv = 2; hh = 1;

for k = 1:1:size(RA,3)

    if ( hh > size(Nl,2) )

        break; % End test sequence.
    end

    X = [ uu, vv ];    

    for j = 1:1:size( X, 2 )

        cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\Data\Excel Data");
    
        dataSet = readmatrix( 'testRGB.csv' ); % Organized test sequence.
    
        L = size(dataSet,1);
    
        % We can grab all observations in order, no matter the order of the
        % data content.
        
        ii = 1;
        for i = 1:L
            if ( dataSet( i, C ) == X( 1, j ) )

                dataSet_( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1;
            end
        end
    
        % dataSet = readmatrix( 'randomizedPhotos.csv' ); dataSet_ = dataSet; % Random assortment of images.
        clear dataSet
        % dataSet = dataSet_( 1:N * Nr * Mr, 1:C );
        
        dataSet = dataSet_( 1:Nl( 1, hh ) * Nr * Mr, 1:C ); % We can iterate through all class images we desire.
    
        skinObservation = dataSet( :, C );
        
        totalN = size( dataSet, 1 );     
    
        trainingN = floor( 0.0 * totalN ); 
        
        testN = floor( 1.0 * totalN );
         
        cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\skinGradientDecent_V5\MATLAB Code");
        
        [ D, E ] = skinImageClassification( dataSet, skinObservation, trainingN, testN, C );

        if( hh <= size(Nl,2) )
        
            [ PREC( hh ), REC( hh ), ACC( hh ), F1( hh ) ] = fMeasure( D, E ); 
        end

        hh = hh + 1;        
    end

    uu = uu + 2; vv = vv + 2;    
end

AVE = [ mean(PREC) mean(REC) mean(ACC) mean(F1) ];

toc;