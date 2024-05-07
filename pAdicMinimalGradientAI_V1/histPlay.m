clc; clear all; close all; tic

global frameLength classType classGroups BOOST RA B_...
    DATARANGE BINS Randomized Supervision Noise...
    BPI ii ij

frameLength = 1250; classType = [0 0]; classGroups = [0 0 0]; 

DATARANGE = 256; BINS = 16;

Randomized = 0; Supervision = 0; Noise = 0;

% Objective labels...

l = [ 1; 1; 1; 1; 2; 2; 2; 2; 3; 3; 3; 3; 4; 4; 4; 4 ];
  
N = size(l,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = readmatrix('trainRGB (1).csv'); N = floor ( A / frameLength );

A_ = histogramization(A,N,l); A_ = mean(A_(:,:,:),3);

% A_ = complexHistogramization(A,N,[]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_ = zeros(size(A_,1),size(A_,2),size(A_,3));

D_ = zeros(size(A_,1),size(A_,2),size(A_,3));

for k =1:1:size(A_,3)
   for i =1:1:size(A_,1)
    
       [C_(i,:,k), D_(i,:,k)] = sort(A_(i,:,k),2,'descend');
   end
end

for k = 1:1:size(A_,3)
    
   G(:,:,k)  = C_(:,:,k);
   
   Gj(:,:,k) = D_(:,:,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BOOST = 1;

RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)-2+BOOST);

[ RA ] = CNN( A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ RA ] = filterOptimization( RA, A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BOOST = 0;

RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)+BOOST);

[ RA ] = GMRES( A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Noise = 1;

if( Noise )
        
    AWGN = rand([size(A_,1),size(A_,2)]);
        
    A_ = A_ + AWGN;
 end

[ A ] = denoise( A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

decisionTree(); % Ensemble

toc
