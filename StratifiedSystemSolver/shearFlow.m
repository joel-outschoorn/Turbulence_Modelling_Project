function [U,DU,DDU,ymesh] = shearFlow(flow,N)

    % ==================== INTRODUCTION TO FUNCTION =======================
    %
    % Function to determine velocity profile of chosen the shear flow
    % Plots of the velocity and stratification profile are also displayed
    %
    % INPUTS:
    % ( 1 ) Chosen shear flow [ flow ]
    % ( 2 ) Number of Chebyshev grid points [ N ]    
    % 
    % OUTPUTS:
    % ( 1 ) Velocity profile [ U ]
    % ( 2 ) Velocity gradient profile [ DU ]
    % ( 3 ) Second derivative of velocity profile [ DDU ]
    % ( 4 ) Domain grid points between -1 to 1 [ ymesh ]
   
    % -------------------------- BEGIN FUNCTION ---------------------------  
    
    % Grid points
    y = sin(pi*(N-1:-2:1-N)'/(2*(N-1))); % Compute Chebyshev points
    ymesh = y(2:N-1); % Ignoring boundary points 

    % Shear flow either pCf or pPf
    % User may add other shear flows by adding further cases
    switch flow 
        case 1 % pCf
        
            % pCf velocity profiles
            U    =      diag(ymesh); 
            DU   =     eye(N-2,N-2); 
            DDU  =   zeros(N-2,N-2); 
            
            % Displaying velocity and stratification profile to user
            disp('  ')
            disp('Velocity profile:               Density profile:')
            disp('(x,y) plane                     (X,Y) plane   ')
            
            v_pCf = pCf_plot; % Plot of velocity and stratified profile
            
        case 2 % pPf

            % pPf velocity profiles
            U    =   diag(1 - ymesh.^2); 
            DU   =       diag(-2*ymesh); 
            DDU  =      -2*eye(N-2,N-2);
            
            % Displaying velocity and stratification profile to user
            disp('  ')
            disp('Velocity profile:               Density profile:')
            disp('(x,y) plane                     (X,Y) plane   ')
            
            v_pPf = pPf_plot; % Plot of velocity and stratified profile
            
    end
    
    % Functions to plot of velocity and stratified profile
    % Profiles are hard coded - do not edit these functions
    
    % pCf plot
    function v_pCf = pCf_plot
        
            % Command window plot
            v_pCf = char(ones(12,55));
            v_pCf(:) = ' ';

            for i = 1:6
                v_pCf(i,11:22-(2*i)) = '-';
                v_pCf(i,22-((2*i)-1)) = '>';
            end
            
            for i = 7:11
                v_pCf(i,22-((2*i-2))) = '<';                 
                v_pCf(i,11-2*(i-7):10) = '-';
            end
            
            for i = 1:11
                v_pCf(i,33:32 + 2*i) = '-';
                v_pCf(i,32 + 2*i) = '>';
            end
            
            % Walls
            v_pCf(1,1:19) = '=';
            v_pCf(12,1:19) = '=';
            v_pCf(1,20:21) = ' ';
            
            v_pCf(1,33:54) = '=';
            v_pCf(12,33:54) = '=';
            
            % Final plot
            disp(v_pCf);
        
    end
    
    % pPf plot
    function v_pPf = pPf_plot
        
            % Command window plot
            v_pPf = char(ones(12,55));
            v_pPf(:) = ' ';       
            
            % Parabola profile
            v_pPf(2,9) = '>';
            v_pPf(3,13) = '>';
            v_pPf(4,17) = '>';
            v_pPf(5,19) = '>';
            v_pPf(6,20) = '>';
            
            for i = 2:6 
               for j = 1:20
                   if ( v_pPf(i,j) == '>')
                       break   
                   else
                        v_pPf(i,j) = '-';
                   end
               end
            end
            
            % Parabolic symmetry
            for i = 7:11
                v_pPf(i,:) = v_pPf(13-i,:);
            end   
            
            for i = 1:11
                v_pPf(i,33:32 + 2*i) = '-';
                v_pPf(i,32 + 2*i) = '>';
            end  
            
            % Walls
            v_pPf(1,1:20) = '=';
            v_pPf(12,1:20) = '=';          
            
            v_pPf(1,33:54) = '=';
            v_pPf(12,33:54) = '=';
            
            % Final plot
            disp(v_pPf);
        
    end
       
    %-------------------------- END OF FUNCTION ---------------------------
    
end