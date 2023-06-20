function [ Z ] = imageDecision( Y ) 

    global V

    if ( V )

        Z = gradientDecentS( Y );       
    else 

        Z = gradientDecentUS( Y );    
    end

end