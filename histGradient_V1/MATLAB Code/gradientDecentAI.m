clear all; close all; clc; 

tic;

cd("C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code"); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global classType classGroups frameLength DATARANGE BINS Supervision Randomized Optimized Noise

Classes = 8; % Classes per group...

classType = zeros( 1, Classes ); 

Groups  = 3; % Groups per classification (RGB)...

classGroups = zeros( 1, Groups );

Np = 25; Mp = 25; frameLength = Np * Mp; % Photo length, and width. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATARANGE = 2^8;

ii = 1;
for i = 1:1:DATARANGE   
    if( mod(DATARANGE,i) == 0 )
        
        BINS_(ii,1) = i; ii = ii + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BINS        = BINS_(6,1); % Histogram bins...

Supervision = 0; % Supervision...

Randomized  = 0; % Frame randomization...

Optimized   = 0; % RA optimization...

Noise       = 0; % RA / input noisy...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of images per class to classify.

for i = 1:1:size(numImages,2)

    N(i,1) = numImages(i); % Number of objects per class.     
end
 
% We can generate an objective label vector to keep track of our errors
% with supervised / unsupervised data...

Observation = verificationList( N );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd("C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");

if ( Supervision )

    dataSet = readmatrix( 'trainRGB (1).csv' ); % Supervised training data.

    % dataSet = readmatrix( 'trainRGB (3).csv' ); % Supervised training data.

elseif( ~Supervision )

    dataSet = readmatrix( 'trainRGB (2).csv' );   % Unsupervised training data.

    % dataSet = readmatrix( 'trainRGB (4).csv' );   % Unsupervised training data.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd("C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");

classifierTraining( dataSet, N, Observation );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test... 

PERSISTENCE = 1; T = 100;

while( PERSISTENCE < T )

    cd("C:\Users\johnbodhi\Documents\GitHub\gradient-decent\histGradient_V1\Data\Excel Data");
    
    % Test sequences not included in training data...
    
    % dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.
    
    dataSet = readmatrix( 'testRGB (1).csv' ); % Unsupervised test sequence. 

    % dataSet = readmatrix( 'testRGB (2).csv' ); % Unsupervised test sequence. 

    % dataSet = readmatrix( 'testRGB (3).csv' ); % Unsupervised test sequence. 

    % dataSet = readmatrix( 'testRGB (4).csv' ); % Unsupervised test sequence. 
    
    cd("C:\Users\johnbodhi\\Documents\GitHub\gradient-decent\histGradient_V1\MATLAB Code");
    
    dataSet = histogramization( dataSet, N, Observation );  
    
    [ D, E ] = classifier( dataSet, Observation );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ PREC(T,1) REC(T,1) ACC(T,1) F1(T,1) ] = fMeasure( D, E );
    
    AVE = [ mean(PREC(:,1)) mean(REC(:,1)) mean(ACC(:,1)) mean(F1(:,1)) ];
    
    T = T + 1;
end

toc;