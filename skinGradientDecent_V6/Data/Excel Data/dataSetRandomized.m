function [ X ] = dataSetRandomized(dataSet,L,C)

    if ( isfile( 'dataSetRandomized.csv' ) )

        X = readmatrix('dataSetRandomized.csv');

    else

        dataSetRandomized = zeros( L, C ); 

        i = zeros(L,1);

        for k = 1:L

            i( k ) = randi( [1, L ] );

            while ( find( i( 1:k-1 ) == i( k ) ) )
                i( k ) = randi( [ 1, L ] );      
            end

            dataSetRandomized( k, 1:C )  = dataSet( i( k ), 1:C );
        end        

            writematrix(dataSetRandomized,'dataSetRandomized.csv');

            X = dataSetRandomized;
    end    
end