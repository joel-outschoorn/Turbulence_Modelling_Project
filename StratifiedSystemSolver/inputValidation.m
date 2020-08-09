function  [flow,plotFigure,N,alpha,beta,Re,Fh,theta,Sc] = inputValidation(flow,plotFigure,N,alpha,beta,Re,Fh,theta,Sc)

    % ==================== INTRODUCTION TO FUNCTION =======================
    %
    % Function to read, validate and display the main parameters inputs
    %
    % INPUTS:
    % ( 1 ) Chosen shear flow ( pCf = 1; pPf = 2 ) [ flow ]  
    % ( 2 ) Request whether to produce plots ( Yes = 1; No = 0 ) [ plotFigure ]
    % ( 3 ) Number of Chebyshev grid points [ N ]                                                                                      
    % ( 4 ) Wavenumber ( x-direction ) [ alpha ]                                             	
    % ( 5 ) Wavenumber ( z-direction ) [ beta ] 
    % ( 6 ) Reynolds number [ Re ]  
    % ( 7 ) Froude number ( horizontal ) [ Fh ]            	
    % ( 8 ) Shear/stratification unalignment angle [ theta ] ( degrees )
    % 
    % OUTPUTS:
    % Fully validated inputs
    
    % -------------------------- BEGIN FUNCTION ---------------------------
    
    while ( flow ~= 1 && flow ~= 2 ) % Validating shear flow
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen shear flow is either 1 (pCf) or 2 (pPf)!')
        flow = input('Chosen shear flow (pCf or pPf): ');
    end
    
    while ( plotFigure ~= 0 && plotFigure ~= 1 ) % Validating request to plot figures
        disp('ERROR FROM INPUT')
        disp('Please ensure that plotFigure is either 0 (do not plot) or 1 (plot)')
        plotFigure = input('Would you like the eigenspectra/eigenfunction to be plotted: ');
    end
    
    while ( N ~= round(N) || N <= 2 ) % Validating grid points
        disp('ERROR FROM INPUT')
        disp('Please ensure the number of grid points is an appropriate value!')
        N  =  input('Number of grid points: ');
    end   
    
    while ( alpha < 0 ) % Validating wave number (x-dir)
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen wave number (x-dir) is an appropriate value!')
        alpha = input('Chosen wave number (x-dir): ');
    end  
    
    while ( beta < 0 ) % Validating wave number (z-dir)
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen wave number (z-dir) is an appropriate value!')
        beta = input('Chosen wave number (z-dir): ');
    end  
    
    while ( Re <= 0 ) % Validating Reynolds number
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen Reynolds number is an appropriate value!')
        Re = input('Chosen Reynold number: ');
    end
    
    while ( Fh <= 0 ) % Validating Froude number
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen Froude number is an appropriate value!')
        Fh = input('Chosen Froude number: ');
    end    
    
    while ( Sc <= 0 ) % Validating Schmidt number
        disp('ERROR FROM INPUT')
        disp('Please ensure the chosen Schmidt number is an appropriate value!')
        Sc = input('Chosen Schmidt number: ');
    end    
    
    % Theta requires no validation
    
    % Displaying value to user
    disp('User Inputs: ')
    
    switch flow
        case 1
            fprintf(1, '%s:          %s \n', 'flow', 'plane Couette flow')
        case 2
            fprintf(1, '%s:          %s \n', 'flow', 'plane Poiseuille flow')
    end    
    
    switch plotFigure
        case 0
            fprintf(1, '%s:    %s \n', 'plotFigure', 'Eigenspectra/Eigenfunctions will not be plotted')
        case 1
            fprintf(1, '%s:    %s \n', 'plotFigure', 'Eigenspectra/Eigenfunctions to be plotted')
    end   
    
    fprintf(1, '%s:             %d \n', 'N', N)
    fprintf(1, '%s:         %.3f \n', 'alpha', alpha)
    fprintf(1, '%s:          %.3f \n', 'beta', beta)
    fprintf(1, '%s:            %.3f \n', 'Re', Re)
    fprintf(1, '%s:            %.3f \n', 'Fh', Fh)
    fprintf(1, '%s:         %.3f %s \n', 'theta', theta, 'degrees')
    disp('  ')
    disp('Constant parameter:')
    fprintf(1, '%s:             %.3f \n', 'Sc', Sc)
    
    %-------------------------- END OF FUNCTION ---------------------------
    
end