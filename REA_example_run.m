%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Envelope Analysis (REA)
%
% Created by Daniel Du, 05/20/2017
% Proteomics and Metabolimics Core Facility
% MD Anderson Cancer Center, Houston, TX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note

%This script performs response envelope analysis on a published data set. 
%The data set can be found in Griner PNAS 2014. Response envelope analysis 
%(REA) is a tool to quantitatively determine combination effects including
%synergy, additivity, and antagonism. 

%Borisy, A.A. et al, PNAS 2003. doi: 10.1073/pnas.1337088100
%Griner, L.A.M. et al, PNAS 2014. doi: 10.1073/pnas.1311846111

%--------------------------------------------------------------------------
clear
clc

%% 1 general parameters
trim = 0.1; %trim the x scale and y scale this amount beyond the data points for 2-D representation
ft = 9; %font size
lw = 3; %line width

%% 2 read data file
data = csvread('Borisy_PNAS_2003.csv'); 
%data = csvread('Griner_PNAS_2014.csv'); 

%the data are formatted as the following:
%column 1, concentrations of drug 1, in uM
%column 2, concentrations of drug 2, in uM
%column 3, corresponding survival rates in %
%if the data show affected ratio or drug response instead, then the 
%following step is needed
%data(:,3) = 100-data(:,3); 

%% 3 perform the analysis
f = REA_package(data,trim,ft,lw,'pentamidine','chlorpromazine',1); 
%f is defined as a vector that contains SI and AI