% clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global classType classGroups imageLength BINS Supervision

Classes = 8; % Classes per group...

classType = zeros( 1, Classes ); 

Groups  = 3; % Groups per classification...

classGroups = zeros( 1, Groups );

Np = 100; Mp = 100; imageLength = Np * Mp; % Photo length, and width. 

BINS = 20; % Histogram bins...

% Number of images per class to classify.

N = zeros(size(numImages,2),1);

for i = 1:1:size(N,1)

    N(i,1) = numImages(i); % Number of objects per class. 
end
totalN = sum(N); 

% We can generate an objective label vector to keep track of our errors
% with unsupervised data...

Observation = verificationList( N, totalN );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

Supervision = 1; % Supervision...

if ( Supervision )

    dataSet = readmatrix( 'trainRGB (1).csv' );   % Supervised training data.

elseif( ~Supervision )

    dataSet = readmatrix( 'trainRGB (2).csv' );   % Unsupervised training data.
end

% dataSet = randomizeAll( dataSet, N ); % Randomize all frames.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

classifierTraining( dataSet, Observation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test...      

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

% Test sequences not included in training data...

% dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.

dataSet = readmatrix( 'testRGB.csv' ); % Unsupervised test sequence. 

% dataSet = randomizeAll( dataSet, N ); % Randomize all photos.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

dataSet = histogramization( dataSet, Observation ); 

Observation = dataSet( :, size(dataSet,2) ); % Supervised observations.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");    

[ D, E ] = classifier( dataSet, Observation );

[ PREC REC ACC F1 ] = fMeasure( D, E );

toc;