clear all; close all; clc; 

tic;

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\Data\Image Data");

photoToArray(); % Pre-process all images.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\MATLAB Code"); 

% We can use K-Means clustering, and a RGB running average within a Gradient Decent / Ascent.

global classType classGroups imageLength

Classes = 2; % Classes per group...

classType = zeros( 1, Classes ); 

Groups = 4; % Groups per classification...

classGroups = zeros( 1, Groups );

Np = 40; Mp = 40; imageLength = Np*Mp; % Photo length, and width. 

% Number of images per class to classify.

N = zeros(size(numImages,2),1);

for i = 1:1:size(N,1)

    % N(i,1) = numImages(i); % Number of objects per class.
    N(i,1) = 5; 
end
totalN = sum(N); 

% We can generate an objective label vector to keep track of our errors
% with unsupervised data...

verObservation = verificationList( N, totalN );

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\Data\Excel Data");

dataSet = readmatrix( 'trainRGB.csv' );   % Supervised training data.

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\MATLAB Code");

pixelClassifierTraining( dataSet );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Verification / Test...

hh = 1; cc = 1; T = 0;       

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\Data\Excel Data");

% Test sequences not included in training data...

% dataSet = readmatrix( 'verificationRGB.csv' ); % Supervised test sequence.

dataSet = readmatrix( 'testRGB.csv' ); % Unsupervised test sequence.   

% We can randomize all data frames over all groups.

% dataSet = randomizeAll( dataSet, Np, Mp, N ); % Randomize all photos.

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

% Unsupervised test scoop.

ii = 1;
for i = ( 1 + T ):1:( N(cc,1)*imageLength + T ) % We can choose more than one frame per class.
    if ( dataSet( i, C ) == 0 )

        dataSet_( ii, 1:C ) = dataSet( i, 1:size(dataSet,2) ); ii = ii + 1;
    end
end    
T = T + N(cc,1)*imageLength; 

% We can randomize all frames within each scoop for N > 1...        

M = N(cc,1); cc = cc + 1;

%        dataSet_ = randomizeClass( dataSet_, Np, Mp, L, C, M ); % We can randomize all frames in each class.
% 
%        dataSet_ = randomizeGroup( dataSet_, Np, Mp, L, C, M ); % We can randomize all frames in each group.

clear dataSet

dataSet = dataSet_( :, 1:size(dataSet,2) ); 

testObservation = dataSet( :, size(dataSet,2) ); % Supervised observations.

clear dataSet_

cd("C:\Users\johnm\OneDrive\Documents\GitHub\gradient-decent\decisionTreeGradient_V2\MATLAB Code");    

[ D, E ] = imageClassification( dataSet, testObservation, verObservation );

if( hh <= size( Nl, 1 ) )

    [ PREC( hh ), REC( hh ), ACC( hh ), F1( hh ) ] = fMeasure( D, E ); 
end

hh = hh + 1;    

AVE = [ mean(PREC) mean(REC) mean(ACC) mean(F1) ];

toc;