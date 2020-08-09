function [omega_stratSort,eFunction,ymesh] = OS_SQ_Stratified_Solver(flow,plotFigure,N,alpha,beta,Re,Fh,theta)

    % ==================== INTRODUCTION TO FUNCTION =======================
    %
    % Linear Stability Analysis of a density stratified plane shear flow
    % with a linear stratification and shear inalignment by angle theta
    %
    % Shear flows:
    % pCf - plane Couette flow
    % pPf - plane Poiseuille flow
    %
    % Edited by: Joel Outschoorn
    % {Ref 1}: Outschoorn; Transition In Stratified Shear Flows: Linear 
    % Stability Analysis With A Shear And Stratification Inalignment
    %
    % {Ref 2}: Schmid and Henningson; Stability and Transition in Shear Flows
    %
    % ====================== Main user inputs =============================                
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
    % ( 1 ) Complex frequency sorted by the growth rate from least stable
    %       eigenmodes to most stable [ omega_stratSort ]
    % ( 2 ) Eigenfunctions of the most unstable eigenmode [ eFunction ]
    % ( 3 ) Domain grid points between -1 to 1 [ ymesh ]
    %
    % ================= Governing eigenvalue problem ======================
    %
    % Eigenvalue problem - Equations (5.23) and (5.24) in {Ref 1}
    %
    % L_strat*q_strat = c_strat*M_strat*q_strat
    %
    % L_strat = [ Los,    zeros,    Losb;
    %             zi*beta*DU,    Lsq,    zi*alpha*sind(theta);
    %             Lv,    (sind(theta)/(Fh^2*k2))*zi*alpha,    Lb ];
    %
    % M_strat = [(k2*I - D2), zeros, zeros;
    %                  zeros,   eye, zeros;
    %                  zeros, zeros,  eye ];
    %
    % q_strat = [v_tilde eta_tilde b_tilde]
    %
    % c_strat -> eigenvalue
    % c_strat = i*omega_strat -> omega_strat = -i*c_strat
    % omega_strat is the growth rate
    %
    % ================ Unstratified system (reference) ====================
    %
    % Unstratified System - Equations (3.33) - (3.35) in {Ref 2}
    %
    % Orr-Somm
    % M_os = (k2*I - D2)
    % L_os = Los
    %
    % Orr-Somm and Squire
    % M = [(k2*I - D2), zeros; zeros, eye]
    % L = [ Los, zeros; zi*beta*DU, Lsq]
    %
    % ======================= Velocity parameters =========================
    %
    % Velocity profile U
    % Velocity gradient DU
    % Second velocity derivative DDU
    % pCf - linear velocity profile between -1 to 1
    % pPf - parabolic velocity profile with no slip at each boundary, max
    % velocity 1 at the midpoint
    
    % -------------------------- BEGIN FUNCTION ---------------------------
    
    % Figure formatting
    set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');

    % Introduction
    disp('===============================================================')
    disp('           LINEAR STABLILTY ANALYSIS OF AN UNALIGNED           ')
    disp('          DENSITY STRATIFICATION TO A PLANE SHEAR FLOW         ')
    disp('===============================================================')

    %-------------------- Validating input parameters ---------------------

    disp('---------------------- INPUT PARAMETERS -----------------------')
    
    Sc = 700;               % Schmid number

    % Validation and displaying value
    [flow,plotFigure,N,alpha,beta,Re,Fh,theta,Sc] = ...
        inputValidation(flow,plotFigure,N,alpha,beta,Re,Fh,theta,Sc);

    % Shear flow velocity profile
    [U,DU,DDU,ymesh] = shearFlow(flow,N);
    
    % Constants
    zi = sqrt(-1);          % Imaginary number
    I  = eye(N-2);          % Identity matrix
    
    k2 = alpha^2 + beta^2;  % k constant

    disp('---------------------------------------------------------------')

    %------------------- Computing Chebyshev derivatives ------------------

    % Functions to compute Chebyshev matrices provided by Dr. Hwang
    [~, D1] = chebdif(N,1);         % Function to compute 1st derivative matrix
    D1      = D1(2:N-1, 2:N-1);     % First derivative ignoring boundaries
    [~, D2] = chebdif(N,2);         % Function to compute 2nd derivative matrix
    D2      = D2(2:N-1, 2:N-1, 2);  % Second derivative ignoring boundaries
    [~, D4] = cheb4c(N);            % Function to compute 4th derivative matrix

    %-------------------- Assembling eigenvalue system --------------------

    % Terms in matrix system, Equations (5.24)
    Los  = zi*alpha*U*(k2*I - D2) + zi*alpha*DDU + (1/Re)*(k2*k2*I - 2*k2*D2 + D4);
    Losb = (k2*I*cosd(theta) - zi*beta*D1*sind(theta))/(Fh^2);
    Lsq  = zi*alpha*U + (1/Re)*(k2*I - D2);
    Lv   = (sind(theta)/(Fh^2*k2))*zi*beta*D1 - (cosd(theta)/(Fh^2))*I;
    Lb   = (zi*alpha*U + (1/(Re*Sc))*(k2*I - D2))/(Fh^2);

    % Coefficients to L and M matrix
    ZERO = zeros(N-2);
    ONE   = eye(N-2);

    % M matrix
    M_strat = [(k2*I - D2), ZERO, ZERO;
                ZERO, ONE, ZERO;
                ZERO, ZERO, ONE ];

    % L matrix
    L_strat = [ Los, ZERO, Losb;
                zi*beta*DU, Lsq, (zi*alpha*sind(theta)*I)/(Fh^2);
                Lv, (sind(theta)/(Fh^2*k2))*zi*alpha*I, Lb ];

    %----------------- Solving and sorting eigenvalues --------------------

    [V,c_strat] = eig(L_strat,M_strat);  % System e/spectra and e/functions
    c_strat     = diag(c_strat);         % Returning eigenmodes into vector form                    
    omega_strat = -zi*c_strat;           % Stratified growth rate

    % Sorting e/values from least stable to most stable in terms of the
    % imaginary part
    omega_stratSort = eigenSort(omega_strat);  
    
    % Finding eigenfunction for most unstable eigenmode
    % if statement to find whether most unstable eigenmode has a conjugate
    % pair
    if abs(imag(omega_stratSort(1)) - imag(omega_stratSort(2))) < 1e-5
        
        % Position of most unstable eigenmode in unsorted eigenspectra
        position1 = find(imag(omega_strat) == imag(omega_stratSort(1)));
        % Position of most unstable conjugate pair
        position2 = find(imag(omega_strat) == imag(omega_stratSort(2)));
        
        % Eigenfunction of the most unstable eigenmode
        eFunction(:,1) = V(:,position1); 
        % Eigenfunction of the conjugate pair
        eFunction(:,2) = V(:,position2); 
        
        % Informing the user of the pair of unstable eigenmodes
        disp('------------------------ Information --------------------------')
        disp('        The most unstable eigenmode has a conjugate pair       ')
        disp('       Eigenfunctions of both eigenmodes will be plotted       ')
        disp('---------------------------------------------------------------')
    else
        
        % Position of most unstable eigenmode in unsorted eigenspectra
        position = find(imag(omega_strat) == imag(omega_stratSort(1)));
        
        % Eigenfunction of the most unstable eigenmode
        eFunction = V(:,position); 
        
    end
    
    % Plotting eigenspectra and eigenfunction should the user request it
    if plotFigure == 1
       
        % Eigenspectra plot
        figure
        plot(omega_stratSort,'kx','Markersize',10,'Linewidth',1.5)
        hold on
        % Highlighting the most unstable eigenvalues
        plot(omega_stratSort(1),'ro','Markersize',15,'Linewidth',1.5)
        % Conjugate pairs
        if abs(imag(omega_stratSort(1)) - imag(omega_stratSort(2))) < 1e-5
            plot(omega_stratSort(2),'ro','Markersize',15,'Linewidth',1.5)
        end
        % Highlighting x and y axis
        plot([-2 2],[0 0],'k-'); plot([0 0],[-2 2],'k-');
        grid minor; box on; set(gca,'Fontsize',20);
        ylim([-1 0.1]);
        xlabel('Real($\omega$)'); ylabel('Imag($\omega$)')
        title('Eigenspectra')
        
        % Figure title
        figTitle = {'Eigenfunction', 'Eigenfunction of conjugate pair'};
        for i = 1:size(eFunction,2)
            % Eigenfunction plot
            figure
            % v-perturbation plot
            subplot(1,3,1)
            plot(abs(eFunction(1:end/3,i)),ymesh,'k','Linewidth',1.5)
            hold on
            plot(real(eFunction(1:end/3,i)),ymesh,'b--','Linewidth',1.5)
            plot(imag(eFunction(1:end/3,i)),ymesh,'r--','Linewidth',1.5)
            legend('Abs($\omega$)','Real($\omega$)','Imag($\omega$)')
            xlabel('$\tilde{v}$')
            ylabel('$y$')
            grid minor; box on; set(gca,'Fontsize',25)
            title('$\tilde{v}$-Perturbation Profile')
            
            % eta-perturbation plot
            subplot(1,3,2)
            plot(abs(eFunction(end/3+1:2*end/3,i)),ymesh,'k','Linewidth',1.5)
            hold on
            plot(real(eFunction(end/3+1:2*end/3,i)),ymesh,'b--','Linewidth',1.5)
            plot(imag(eFunction(end/3+1:2*end/3,i)),ymesh,'r--','Linewidth',1.5)
            xlabel('$\tilde{\eta}$')
            grid minor; box on; set(gca,'Fontsize',25)
            title('$\tilde{\eta}$-Perturbation Profile')

            % b-perturbation plot
            subplot(1,3,3)
            plot(abs(eFunction(2*end/3+1:end,i)),ymesh,'k','Linewidth',1.5)
            hold on
            plot(real(eFunction(2*end/3+1:end,i)),ymesh,'b--','Linewidth',1.5)
            plot(imag(eFunction(2*end/3+1:end,i)),ymesh,'r--','Linewidth',1.5)
            xlabel('$\tilde{b}$')
            grid minor; box on; set(gca,'Fontsize',25)
            title('$\tilde{b}$-Perturbation Profile')
            sgtitle(figTitle(i))
        end
   
    end
    
    %-------------------------- END OF FUNCTION ---------------------------

end
