
% ====================== INTRODUCTION TO PROGRAM ==========================
%
% Program to run different test cases of the Orr-Sommerfeld (OS) and Squire
% (SQ) Stratified Solver function: OS_SQ_Stratified_Solver.m
% Used to produce useful results and plots for analysis of the statiblity
% of an inaligned density stratification to a plane shear flow
%
% ====================== Main parameter inputs ============================                
%
% OS_SQ_Stratified_Solver.m function:
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

% ---------------------------- BEGIN PROGRAM ------------------------------

clear
clc
close all

% Example input parameters
pCf        = 1; 
pPf        = 2;
plotFigure = 1;
N          = 100;
alpha      = 1;
beta       = 0;
Re         = 1000;
Fh         = 1;
theta      = 30;

% Example simulation for pCf
[omega_pCf,eFunction_pCf,ymesh_pCf] = OS_SQ_Stratified_Solver(pCf,plotFigure,N,alpha,beta,Re,Fh,theta);

% Example simulation for pPf
[omega_pPf,eFunction_pPf,ymesh_pPf] = OS_SQ_Stratified_Solver(pPf,plotFigure,N,alpha,beta,Re,Fh,theta);

%---------------------------- END OF PROGRAM ------------------------------

