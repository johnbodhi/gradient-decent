function [ RA ] = SVM( dataSet, N, Observation )

    global A RA Supervision Randomized Optimized Train

    A = histogramization( dataSet, N, Observation ); 

    if ( ( ~Supervision && Train && ~Randomized ) ||...
         ( Supervision && Train && Randomized ) )

        [ RA ] = filterCreation( RA, A );       
        
    elseif( ( ~Supervision && ~Train && ~Randomized ) ||...
            ( Supervision && ~Train && Randomized ) ) 
        
        [ RA ] = filterUpdate( dataSet );

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( A, Observation );
    end

end