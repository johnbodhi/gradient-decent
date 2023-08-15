function plotDist(RA)

    % Instead of weigting, we can flip and dope our distributions...

    T = (1:1:256);

    figure( 'name', "The Histogram Gradient Method");

    subplot(3,1,1);
    plot(T(1:1:256),RA(7,:,1),T(1:1:256),RA(8,:,1));
    title("The ol' Dope 'n Flip: Step 1 (RA 7,:,1 V.S. RA 8,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    subplot(3,1,2);
    plot(T(1:1:256),RA(7,:,1),T(1:1:256),flip(RA(8,:,1),2));
    title("The ol' Dope 'n Flip: Step 2 (RA 7,:,1 V.S. RA 8,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    RA(7,220:256,1) = 0;    
    RA(8,241:256,1) = 0;

    subplot(3,1,3)
    plot(T(1:1:256),RA(7,:,1),T(1:1:256),flip(RA(8,:,1),2));
    title("The ol' Dope 'n Flip: Step 3 (RA 7,:,1 V.S. RA 8,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    figure( 'name', "The Histogram Gradient Method");

    subplot(3,1,1);
    plot(T(1:1:256),RA(5,:,1),T(1:1:256),RA(6,:,1));
    title("The ol' Dope 'n Flip: Step 1 (RA 5,:,1 V.S. RA 6,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    subplot(3,1,2);
    plot(T(1:1:256),RA(5,:,1),T(1:1:256),flip(RA(6,:,1),2));
    title("The ol' Dope 'n Flip: Step 2 (RA 5,:,1 V.S. RA 6,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    RA(5,1:7,1) = 0;    
    RA(6,1:7,1) = 0;

    subplot(3,1,3)
    plot(T(1:1:256),RA(5,:,1),T(1:1:256),flip(RA(6,:,1),2));
    title("The ol' Dope 'n Flip: Step 3 (RA 5,:,1 V.S. RA 6,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );

  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure( 'name', "The Histogram Gradient Method");

    subplot(3,1,1);
    plot(T(1:1:256),RA(3,:,1),T(1:1:256),RA(4,:,1));
    title("The ol' Dope 'n Flip: Step 1 (RA 3,:,1 V.S. RA 4,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    subplot(3,1,2);
    plot(T(1:1:256),RA(3,:,1),T(1:1:256),flip(RA(4,:,1),2));
    title("The ol' Dope 'n Flip: Step 2 (RA 3,:,1 V.S. RA 4,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    RA(3,1:50,1) = 0;    
    RA(4,1:80,1) = 0;

    RA(3,225:256,1) = 0;    
    RA(4,225:256,1) = 0;

    RA(3,50:124,1) = 0; 

    subplot(3,1,3)
    plot(T(1:1:256),RA(3,:,1),T(1:1:256),flip(RA(4,:,1),2));
    title("The ol' Dope 'n Flip: Step 3 (RA 3,:,1 V.S. RA 4,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure( 'name', "The Histogram Gradient Method");

    subplot(3,1,1);
    plot(T(1:1:256),RA(1,:,1),T(1:1:256),RA(2,:,1));
    title("The ol' Dope 'n Flip: Step 1 (RA 1,:,1 V.S. RA 2,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    subplot(3,1,2);
    plot(T(1:1:256),RA(1,:,1),T(1:1:256),flip(RA(2,:,1),2));
    title("The ol' Dope 'n Flip: Step 2 (RA 1,:,1 V.S. RA 2,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );
    
    RA(1,1:50,1) = 0;    
    RA(2,1:50,1) = 0;

    subplot(3,1,3)
    plot(T(1:1:256),RA(1,:,1),T(1:1:256),flip(RA(2,:,1),2));
    title("The ol' Dope 'n Flip: Step 3 (RA 1,:,1 V.S. RA 2,:,1)");
    xlabel("X_n");
    ylabel("H");
    % yline(0);
    legend( 'Target Signal 1', 'Target Signal 2' );

end

