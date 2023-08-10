clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global classType classGroups imageLength BINS Supervision Randomized

Classes = 8; % Classes per group...

classType = zeros( 1, Classes ); 

Groups  = 3; % Groups per classification (RGB)...

classGroups = zeros( 1, Groups );

Np = 100; Mp = 100; imageLength = Np * Mp; % Photo length, and width. 

BINS = 12; % Histogram bins...

% Number of images per class to classify.

for i = 1:1:size(numImages,2)

    N(i,1) = numImages(i); % Number of objects per class.     
end
 
% We can generate an objective label vector to keep track of our errors
% with unsupervised data...

Observation = verificationList( N );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

Supervision = 0; % Supervision...

Randomized  = 0;

if ( Supervision )

    dataSet = readmatrix( 'trainRGB (1).csv' );   % Supervised training data.

elseif( ~Supervision )

    dataSet = readmatrix( 'trainRGB (2).csv' );   % Unsupervised training data.
end

% dataSet = randomizeAll( dataSet, N ); % Randomize all frames.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

classifierTraining( dataSet, Observation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test...      

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

% Test sequences not included in training data...

% dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.

dataSet = readmatrix( 'testRGB.csv' ); % Unsupervised test sequence. 

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

dataSet = histogramization( dataSet, N, Observation ); 

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");    

[ D, E ] = classifier( dataSet );

[ PREC REC ACC F1 ] = fMeasure( D, E );

AVE = [ PREC REC ACC F1 ];

toc;