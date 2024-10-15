clc; clear all; close all; tic

% We can use K-Means clustering, and a Poynting Vector running average
% within a Gradient.

global vectorLength classType L C

vectorLength = 181;

classType = [ 1 2 3 4 ]; % Target class designations.

valClasses = 2;

% We can rondomize the class type distibution;

% classType = randomizeClassType( classType );

% Number of images per class to classify.

N = [ 4 4 4 4 ]; % Total targets.

M = 1; % Class training epochs 1-Total Targets.

for j = 1:1:size( classType, 2 )

    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\fdfdGradient_V1\Data\Excel Data\Sphere");
    
    dataSet = readmatrix( 'Co-PolarizedRCS_Sphere_Train.csv' ); 
    
    L = size( dataSet, 1 ); C = size( dataSet, 2 );

    % randomizedFields = randomizedTargets( dataSet, Nr, Mr, L, C, N ); % Randomize all targets.
    
    % dataSetRandomized = dataSetRandomized( dataSet, L, C );  % Randomize all fields values.
    
    dataSetRandomized = readmatrix('dataSetRandomized.csv');
    
    dataSetRandomized = dataSetRandomized( 1:M * vectorLength, 1:C );  % Assign randomized pixels over the length of target vector.
    
    fieldObservation = dataSetRandomized( :, C ); % Extract randmized observations for training.
    
    totalN = size( dataSetRandomized, 1 ); 
    
    trainingN = floor( 0.05 * totalN ); 
    
    testN = floor( 0.0 * totalN );
    
    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\fdfdGradient_V1\MATLAB Code");
    
    fieldClassifierTraining( dataSetRandomized, fieldObservation, trainingN, totalN );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\fdfdGradient_V1\Data\Excel Data\Sphere");
    
    dataSet = readmatrix( 'Co-PolarizedRCS_Sphere_Test.csv' ); % Organized.

    L = size( dataSet, 1 ); C = size( dataSet, 2 );

    % We can grab all observations in order, no matter the order of the
    % data content.

    ii = 1;
    for i = 1:L
        if ( dataSet( i, C ) == classType( 1, j ) )
            dataSet_( ii, 1:C ) = dataSet( i, 1:C ); ii = ii + 1;
        end
    end

    % Random assortment of targets.

    % dataSet = readmatrix( 'randomizedTargets.csv' ); dataSet_ = dataSet; 

    clear dataSet

    % dataSet = dataSet_( 1:N * Nr * Mr, 1:C );
    
    % We can iterate through all class images we desire.

    dataSet = dataSet_( 1:N( 1, j ) * vectorLength, 1:C ); 

    fieldObservation = dataSet( :, C );
    
    totalN = size( dataSet, 1 );     

    trainingN = floor( 0.0 * totalN ); 
    
    testN = floor( 1.0 * totalN );
     
    cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\fdfdGradient_V1\MATLAB Code");
    
    [ D, E ] = fieldClassification( dataSet, fieldObservation, trainingN, testN, C );

    if( j <= valClasses )

        [ PREC( j ), REC( j ), ACC( j ), F1( j ) ] = fMeasure( D, E ); 
    end

    toc;
end

AVE = [ sum(PREC) sum(REC) sum(ACC) sum(F1) ] ./ valClasses;