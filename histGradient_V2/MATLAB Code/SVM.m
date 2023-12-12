function [ RA ] = SVM( dataSet, N, Observation )

    global Supervision Randomized Optimized Train

    Ar = histogramization( dataSet, N, Observation );

    Ac = complexHistogramization( dataSet, N, Observation );

    if ( ( ~Supervision && Train && ~Randomized ) ||...
         ( Supervision && Train && Randomized ) )

        [ RA ] = filterCreation( Ar );  

        % [ RA ] = filterUpdate( Ar );
        
    elseif( ( ~Supervision && ~Train && ~Randomized ) ||...
            ( Supervision && ~Train && Randomized ) ) 
        
        [ RA ] = filterUpdate( dataSet );

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( A, Observation );
    end

end
