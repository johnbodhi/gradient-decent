clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code"); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global classType classGroups frameLength DATARANGE BINS Supervision Randomized

Classes = 8; % Classes per group...

classType = zeros( 1, Classes ); 

Groups  = 3; % Groups per classification (RGB)...

classGroups = zeros( 1, Groups );

Np = 100; Mp = 100; frameLength = Np * Mp; % Photo length, and width. 

DATARANGE = 256;

ii = 1;
for i = 1:1:DATARANGE   
    if( mod(DATARANGE,i) == 0 )
        
        BINS_(ii,1) = i;ii = ii + 1;
    end
end

BINS = BINS_(9,1); % Histogram bins...

Supervision = 0; % Supervision...

Randomized  = 0; % Frame randomization...

% Number of images per class to classify.

for i = 1:1:size(numImages,2)

    N(i,1) = numImages(i); % Number of objects per class.     
end
 
% We can generate an objective label vector to keep track of our errors
% with supervised / unsupervised data...

Observation = verificationList( N );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

if ( Supervision )

    dataSet = readmatrix( 'trainRGB (1).csv' );   % Supervised training data.

elseif( ~Supervision )

    dataSet = readmatrix( 'trainRGB (2).csv' );   % Unsupervised training data.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

classifierTraining( dataSet, N, Observation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test... 


cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

% Test sequences not included in training data...

% dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.

dataSet = readmatrix( 'testRGB.csv' ); % Unsupervised test sequence. 

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

dataSet = histogramization( dataSet, N, Observation ); 

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");    

[ D, E ] = classifier( dataSet, Observation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ PREC REC ACC F1 ] = fMeasure( D, E );

AVE = [ PREC REC ACC F1 ]; 

toc;