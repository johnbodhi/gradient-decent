clear all; close all; clc; tic

A = readmatrix("C:\Users\jmgar\OneDrive\Documents\GitHub\project-euler\Project_Euler_Problem_67\MaximumPathSumII\MaximumPathSumII\Matrices\triangle_tiny.csv");
% A = readmatrix("C:\Users\jmgar\OneDrive\Documents\GitHub\project-euler\Project_Euler_Problem_67\MaximumPathSumII\MaximumPathSumII\Matrices\triangle_small.csv"); 
% A = readmatrix("C:\Users\jmgar\OneDrive\Documents\GitHub\project-euler\Project_Euler_Problem_67\MaximumPathSumII\MaximumPathSumII\Matrices\triangle_large.csv"); 

N = size(A,1); M = size(A,2);

% N = 10; M = N; % Maximal for N odd...
% for j = 1:M
%     for i = 1:N
%          if( i < j || i == j)
% 
%             A( i, j ) = randi( [ 0, N ] );
%             % A( i, j ) = 0;
%          end
%      end
% end
% A = A';

% Flip matrix A to take advantage of symmetry.

AS = zeros(size(A,1),size(A,2));
for i = 1:1:size(A,1)

    AS(i,1:i) = flip(A(i,1:i),2);
end

% Transpose matrix A to march through natrually.

A  = A';
AS = AS';

% Generate 3-D storage of all branch subspaces (leaf) within decision constraint. 

D = zeros( size(A,1), size(A,2), size(A,1) ); kk = 1;

for k = size( D, 3 ):-1:1
    for j = 1:1:size( D, 2 )
        for i = 1:1:size( D, 1 )               
            
            if( k == size(D,3) )
                cc = 0;
                for j = 1:size(D,2)
                    for i = 1:size(D,1)
                        if( i == j )
                            D( k - cc, j, k ) = k - cc;
                            cc = cc + 1;
                        end
                    end
                end

            elseif( k < size(D,2) || k > 1 )

                hh = size(D,2) - k + 1; 
                dd = 1;
                for p = k:-1:1
                    D( p, dd:hh, k ) = p;
                    dd = dd + 1;
                    hh = hh + 1; 
                end                

            elseif( k == 1 )
                    
                for j = 1:size(D,2)
                    for i = 1:size(D,1)
                        if( i == size(D,1) )
                            D( size(D,1), j, k ) = k;
                        end
                    end
                end
                
            end
        end
    end
end

% Re-assign all true values in matrix A associated with the Trellis 
% indexes to matrix R.

D = flip(D,2);
for k = 1:size(D,3)
    for j = 1:size(D,2)
        for i = 1:size(D,1)
            if( D( i, j, k ) )

                RA( i, j, k ) = A( D( i, j, k ), j );

                RS( i, j, k ) = AS( D( i, j, k ), j );
            end            
        end
    end
end

% We want to march through the first half of the tree, and it's
% reflection to minimize the search for eigenvectors that define the 
% paths through the dyadic tree.

if( mod(N,2) ~= 0 )

    R = cat(3,RA(:,:,1:floor(N/2)+1),RS(:,:,1:floor(N/2)));

    R(:,:,floor(N/2)+2:N) = flip(R(:,:,floor(N/2)+2:N),3);
else

    R = cat(3,RA(:,:,1:floor(N/2)),RS(:,:,1:floor(N/2)));

    R(:,:,floor(N/2)+1:N) = flip(R(:,:,floor(N/2)+1:N),3);
end
D = flip(D,2); % R = flip(R,2);

R = circshift(R,1,1);

% Find all legal paths through the Trellis.

Z = diag( flip( pascal( size( A, 1 ) ), 2 ) ); % Number of paths through the triangle from each bottom index.

% We need to find all binary combinations below M-1. These are the all the
% paths specified by the opposite diagonal of Pascal's matrix.

qq = 0; pp = 0; kk = 1;
for i = 1:1:N

    while( pp < Z(i,1) )

        B = permn([0 1],M-1,kk);

        if( sum( B, 2 ) == qq )

            pp = pp + 1;

            V_(:,pp,i) = B;

        end
        kk = kk + 1;
    end
    qq = qq + 1; pp = 0; kk = 1;
end

% We need to utilize symmetry once again by re-using the 
% first half of the binaries.

% V_ = readmatrix("binaryDirectives");

if( mod(N,2) ~= 0 )

    C_ = flip(V_(:,:,1:floor(N/2)),3);

    V  = cat(3,V_(:,:,1:floor(N/2)+1),C_);    
else

    C_ = flip(V_(:,:,1:floor(N/2)),3);

    V  = cat(3,V_(:,:,1:floor(N/2)),C_);    
end

% This routine will march through every directive given by the binary
% streams, and add all corresponding elements in the tree.

SS = zeros(1,2);

% RF = flip(R,2);

UP = 0; DOWN = 1;

DIRECTION = [UP DOWN]; % [ Toward the vertex. Toward the edge. ]

for kk = 1:1:size(Z,1)
    
    if( DIRECTION(1,1) )
        
        for k = 1:1:size(R,1)
        
            if( R( k, 1, kk ) )
                
                ii_ = k;
                break;
            end
        end
        ii = ii_; jj = 1; % Bottom-to-top...
    
    elseif( DIRECTION(1,2) )
        
        % R = RF; 
        
        ii_ = 2;
        
        ii = ii_; jj = 1; % Top-to-bottom...
     
    end
    
    S = R(ii,jj,kk);

    for j = 1:1:Z(kk)
        
        for i = 1:1:size(V,1)
    
            if( ~V(i,j,kk) )
        
                if( jj <= N - 1 )
        
                    jj = jj + 1; 
                    
                    S = S + R(ii,jj,kk);  
                    
                end           
           
            elseif( V(i,j,kk) )
        
        
                if( ii > N - M + 1 && jj <= N - 1 )
        
                    if( DIRECTION(1,1) )
                        
                        ii = ii - 1; 
                        
                    elseif( DIRECTION(1,2) )
                        
                        ii = ii + 1;
                        
                    end
                    
                    jj = jj + 1;    
                    
                    S = S + R(ii,jj,kk); 
                   
                end
                
            end

        end

        SS(1,1) = S;
    
        if( SS( 1 ) < SS( 2 ) )
    
            SS( 1 ) = 0;

        elseif( SS( 1 ) > SS( 2 ) )
    
            SS( 2 ) = 0; 
            SS = circshift( SS, 1, 2 );
        end
        ii = ii_; jj = 1; 
        
        S = R(ii,jj,kk);
    end

end
H = max(SS);

toc;