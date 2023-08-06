clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global classType classGroups imageLength Supervision

Classes = 8; % Classes per group...

classType = zeros( 1, Classes ); 

Groups  = 1; % Groups per classification...

classGroups = zeros( 1, Groups );

Np = 100; Mp = 100; imageLength = Np * Mp; % Photo length, and width. 

% Number of images per class to classify.

N = zeros(size(numImages,2),1);

for i = 1:1:size(N,1)

    N(i,1) = numImages(i); % Number of objects per class. 
end
totalN = sum(N); 

% We can generate an objective label vector to keep track of our errors
% with unsupervised data...

verObservation = verificationList( N, totalN );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

Supervision = 0; % Supervision

if ( Supervision )

    dataSet = readmatrix( 'trainRGB (1).csv' );   % Supervised training data.

elseif( ~Supervision )

    dataSet = readmatrix( 'trainRGB (2).csv' );   % Unsupervised training data.
end

% dataSet = randomizeAll( dataSet, N ); % Randomize all frames.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

pixelClassifierTraining( dataSet, verObservation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test...      

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

% Test sequences not included in training data...

% dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.

dataSet = readmatrix( 'testRGB.csv' ); % Unsupervised test sequence. 

% dataSet = randomizeAll( dataSet, N ); % Randomize all photos.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

dataSet_ = histogramization( dataSet );

clear dataSet

dataSet = dataSet_( :, 1:size(dataSet,2)-1 ); 

testObservation = dataSet( :, size(dataSet,2) ); % Supervised observations.

clear dataSet_

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");    

[ D, E ] = imageClassification( dataSet, testObservation, verObservation );

if( hh <= size( N, 1 ) )

    [ PREC( hh ), REC( hh ), ACC( hh ), F1( hh ) ] = fMeasure( D, E ); 
end  

AVE = [ mean(PREC) mean(REC) mean(ACC) mean(F1) ];

toc;