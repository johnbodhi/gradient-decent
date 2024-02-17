clc; clear all; close all; tic

global frameLength classType classGroups RA B_...
    DATARANGE BINS Randomized Supervision...
    BPI ii ij N

frameLength = 1250; classType = [0 0]; classGroups = [0 0 0]; 

DATARANGE = 256; BINS = 16; Randomized = 0; Supervision = 0;

N = 16;

% We can read in supervised or unsupervised n-dimensional data.

A  = readmatrix('trainRGB (4).csv');

% We want to homogenize the data, and organize the samples 
% symmetrically or asymmetrically into a single array.

A_ = histogramization(A,N,[]);

C_ = zeros(size(A_,1),size(A_,2),size(A_,3));

D_ = zeros(size(A_,1),size(A_,2),size(A_,3));

for k =1:1:size(A_,3)
   for i =1:1:size(A_,1)
    
       [C_(i,:,k), D_(i,:,k)] = sort(A_(i,:,k),2,'descend');
   end
end

% We can look at the sorted turbulence of our hypersurfaces, 
% as well as the column indexes. This is a random forest.

for k = 1:1:size(A_,3)
    
   G(:,:,k)  = C_(:,:,k); 
   
   Gj(:,:,k) = D_(:,:,k);
end

% We need to boost our filter allocation to have globally 
% defined values in our learning rate.

BOOST = 0;

RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)+BOOST);

% We can pass the sample set into a convolutional version of 
% krylov subspace for a single filter update with a pre existing filter,
% or we can generate an initial condition for a supergroup of classes,
% and proceed with further sample sets with the same method or with
% gradient decent. This is error correction.

% [ RA ] = BiCGSTAB( RA, A_ );

% These are objective labels, however the data may or may not be labeled.

l = [ 1; 1; 1; 1; 2; 2; 2; 2; 3; 3; 3; 3; 4; 4; 4; 4 ]; % Optional labels...

% We need to define initial conditions of our filter without knowing 
% what we are sampling, with or without supervision. For this, we use
% a convolutional sub-gradient process to exhaust all combinations of
% the samples in the averages before classifying against the same data
% for which the filter is being composed. We also continually optimize 
% the filter by deleting unecessary elements belonging to the
% seperate categories, because this saves time and space in the future.

% Now, we need to back propagate on the same data, and in the same
% structure, to amend the filter manifold with a new category. 
% This allows us to extend our range of classifications are far as 
% necessary. 

% We can contiunally optimize the embedding as we introduce 
% new data, and extend our filter manifold to infinity. 

% We also introduce error smoothing in a p-adic tree structure 
% in the event the construction failes to converge properly by 
% failing to avoid repetition in the averages, or yield results
% that depend on overly exclusive averages.

ii = 1; ij = 2;

BPI = [ ii ij ];

BOOST = 1;

RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)+BOOST);

NN = sum((1:1:N));

B_ = size(classGroups,2); % Convergence criterion.

while( B_ )
    
    RA = cat(1,RA,zeros(size(classType,2),BINS,size(classGroups,2)));

    [ RA, TL_ ] = buildManifold( A_, l ); 
    
    for k = 1:1:size(classGroups,2)
        
        if( (find(TL_(:,1,k) == TL_(:,2,k))) &&...
                ( ((sum(TL_(:,1,k) + sum((TL_(:,2,k)) < NN ) ||...
                ( ( sum(TL_(:,1,k) + sum((TL_(:,2,k)) > NN) ) )
            
            % We need to continually inject new data 
            % into suspicious categories that causes row-wise
            % samples to be unused or used repetitively in a binary 
            % categorification.
            
            % [ RA, TL ] = decisionTree();
        else
            
            B_ = B_ - 1;        
        end
    end
end

[ RA ] = filterOptimization(RA, A_, l );

% We now have a preconditioned, or optimized, array that we can 
% send into Krylov Subspace. This may hasten verification, or
% classification of new data.

BOOST = 0;

RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)+BOOST);

RA = RA(2:end,2:end,2:end);

[ RA ] = BiCGSTAB( RA, A_ );

toc



