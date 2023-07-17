function [ Z ] = imageDecision( Y ) 

    % We can implement decision trees on RA for an unknown input signal.

    Z = gradientDecent( Y );

    % We can implement dyadic decision trees...

%     M = size(ZA,2);
% 
%     for i = 1:1:size(ZA,1)
% 
%         X = ZA(:,(i-1)*M/2+1:i*M/2,:);
% 
%         Z_(i,1) = gradientDecent( Y, X );
%     end
%     
%     Z = min(Z_);
    

end