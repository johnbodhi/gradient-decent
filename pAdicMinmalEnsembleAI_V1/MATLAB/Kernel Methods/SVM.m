function [ RA ] = SVM( dataSet )

    global Supervision Randomized Optimized Train

    A_ = histogramization( dataSet, [] ); 

    % A_ = complexHistogramization( dataSet, [] ); 

    % A_ = transitionHistogramization( dataSet, [] ); 

    if ( ( ~Supervision && Train && ~Randomized ) ||...
         ( Supervision && Train && Randomized ) )

        [ RA ] = filterCreation( dataSet );       
        
    elseif( ( ~Supervision && ~Train && ~Randomized ) ||...
            ( Supervision && ~Train && Randomized ) ) 
        
        [ RA ] = filterUpdate( dataSet );

    end

    if ( Optimized )        
    
        [ RA ] = filterOptimization( RA );
    end
    
end