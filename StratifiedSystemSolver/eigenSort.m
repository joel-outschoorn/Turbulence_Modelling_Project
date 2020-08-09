function eValueSort = eigenSort(eValue)

    % ==================== INTRODUCTION TO FUNCTION =======================
    %
    % Function to sort eigenvalues in terms of decreasing imaginary part
    % Sorted in increasing stablity of eigenmodes
    %
    % INPUTS:
    % ( 1 ) Unsorted eigenvalues [ eValue ]
    % 
    % OUTPUTS:
    % ( 1 ) Fully sorted eigenvalues [ eValueSort ]
    
    % -------------------------- BEGIN FUNCTION ---------------------------    

    % sorting e/values in terms of imaginary part
    [~, sortCounter] = sort(imag(eValue));
    
    % sorted e/values in increasing imaginary part
    eValueSort = eValue(sortCounter);
    
    % sorted e/values in decreasing imaginary part (more negative)
    eValueSort = flip(eValueSort);
    
    %-------------------------- END OF FUNCTION ---------------------------    

end