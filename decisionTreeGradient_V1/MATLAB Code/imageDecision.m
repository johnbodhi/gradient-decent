function [ Z ] = imageDecision( Y ) 

    global RA numRA

    % We can implement decision trees on RA for an unknown input signal.

    Z = gradientDecent( Y );

    % We can implement dyadic decision trees over groups...

%     for i = 1:1:numRA
% 
%         X = RA(:,:,i);
% 
%         Z_(i,1) = gradientDecent( Y, X );
%     end
%     
%     Z = min(Z_);
end

