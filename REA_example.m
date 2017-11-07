%Note
%
%This script compares response envelope analysis with other exisiting 
%methods (alpha, beta, gamma, etc) on public datasets from screens. 
%Response envelope analysis (REA) is a tool to quantitatively determine 
%combination effects including synergy, additivity, and antagonism. 
%
%References:
%1. Borisy, A.A. et al, PNAS 2003. doi: 10.1073/pnas.1337088100
%   The pentamidine/chlorpromazine data are from this paper.
%
%2. Griner, L.A.M. et al, PNAS 2014. doi: 10.1073/pnas.1311846111
%   The ibrutinib/MK-2206 data are from this paper
%
%the data are formatted as the following:
%column 1, concentrations of drug 1, in uM
%column 2, concentrations of drug 2, in uM
%column 3, corresponding survival rates in percentage

%--------------------------------------------------------------------------
clear
clc

%% 1 general parameters
trim = 0.1; %trim the x scale and y scale this amount beyond the data points for 2-D representation
ft = 9; %font size
lw = 3; %line width

%% 2 read data file
%file = 'pentamidine_chlorpromazine';
file = 'ibrutinib_nintedanib'; 
data = csvread([file,'.csv']); %read data file

%% 3 perform the analysis
spl = strsplit(file,'_'); drug1 = spl{1}; drug2 = spl{2}; %split the file name to get the names of the two drugs
siai = REA_package(data,trim,ft,lw,drug1,drug2,0); %REA
[alpha,beta,gamma] = RSRA_package(data,-1,1); %RSRA