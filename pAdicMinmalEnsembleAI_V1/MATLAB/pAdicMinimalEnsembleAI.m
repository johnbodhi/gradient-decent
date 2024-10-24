clc; clear all; close all; tic

global frameLength classType classGroups BOOST RA W...
    DATARANGE BINS Randomized Supervision Noise...
    Deterministic Probabalistic...
    GD BiCGSTAB

Deterministic = 1; Probabalistic = 0;

classType = [0 0]; classGroups = [0 0 0]; 

frameLength = 10000; DATARANGE = 256; BINS = 16;

Supervision = 0; Randomized = 0; Noise = 0;

GD = 1; BiCGSTAB = 0; HMM = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = readmatrix(  ' ' ); 

N = floor( size(A,1) / frameLength );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( Deterministic )

    % CNN

    A_ = histogramization(A,N,[]); A_ = mean(A_(:,:,:),3);
    
    % A_ = complexHistogramization(A,N,[]);

elseif ( Probabalistic )

    % HMM / CNN

    A_ = tansitionHistogramization(A,N,[]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if( HMM )

    % S = HMM(A_); 

% if ( fHMM )
    
    % S = fHMM(A_);

% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hypersurfaces...

% C_ = zeros(size(A_,1),size(A_,2),size(A_,3));
% 
% D_ = zeros(size(A_,1),size(A_,2),size(A_,3));
% 
% for k =1:1:size(A_,3)
%    for i =1:1:size(A_,1)
% 
%        [C_(i,:,k), D_(i,:,k)] = sort(A_(i,:,k),2,'descend');
%    end
% end
% 
% for k = 1:1:size(A_,3)
% 
%    G(:,:,k)  = C_(:,:,k);
% 
%    Gj(:,:,k) = D_(:,:,k);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BOOST = 1;
% 
% RA = zeros(size(classType,2)+BOOST,BINS+BOOST,size(classGroups,2)-2+BOOST);
% 
% W = zeros(size(RA,1),size(RA,2),size(RA,3)); W(1,:,1) = 1e1; W(2,:,1) = 1e2;

% [ RA ] = CNN( A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [ RA ] = filterOptimization( RA, A_ );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We can put our system into a dyadic tree structure. This allows us to
% track the loss function over all data, and minimize cost until all 
% samples are classified. The size of the structure and the termination 
% leaf determines the value of N_l ( NxN, Leaf ).

% decisionTree();

toc