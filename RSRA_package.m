function [alpha,beta,gamma] = RSRA_package(data,ca,cb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Surface Regression Analysis
%
% Di Du, Ph.D.
% Department of Bioinfomatics and Computational Biology
% University of Texas MD Anderson Cancer Center, Houston, TX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes
%
%Response Surface Regression Analysis is a set of methods to quantitatively 
%determine the global combination effect between drugs. Three different
%parameters, alpha, beta, and gamma are calculated based on Loewe
%Additivity model, Bliss Independence model,and Gaddum's Non-interaction
%model, respectively. 
%
%Instuctions:
%
%Usually one takes s0 = 0 in order for the method to be non-parametric. See
%Cokol, M., et al, MSB, 2011. One can also use the s0 from package_envelope
%_s0.
%
%References:
%
%1. Greco, W.R. et al, Cancer Res. 1990. PMID: 2386940
%2. Cokol, M. et al, Mol. Syst. Biol. 2011. doi: 10.1038/msb.2011.71
%
%Updates:
%10/6/2017 added beta and gamma

%--------------------------------------------------------------------------
data = sortrows(data,[-2 -1]); %sort rows based on column 2 and then column 1
ndose1 = length(unique(data(:,1))); %number of doses for drug 1
ndose2 = length(unique(data(:,2))); %number of doses for drug 2

%single drug response curve for drug 1
[dose1,id1] = sort(data(1:ndose1,1)); %concentrations of drug 1 
dose1 = dose1(2:end); %truncated concentrations of drug 1 to exclude 0
base1 = dose1(1)^2/dose1(3); 
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv1 = data((data(:,2) == 0),3)/100; %unaffected fraction of drug 1
surv1 = surv1(id1);  %sort surv1
surv1 = surv1(2:end); %truncated unaffected fraction of drug 1 to exclude 0

%single drug response curve for drug 2
[dose2,id2] = sort(data(1:(size(data,1)/ndose2):end,2)); %concentrations of drug 2 
dose2 = dose2(2:end); %truncated concentrations of drug 2 to exclude 0
base2 = dose2(1)^2/dose2(3);
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
surv2 = data((data(:,1) == 0),3)/100; %unaffected fraction of drug 2
surv2 = surv2(id2);  %sort surv2
surv2 = surv2(2:end); %truncated unaffected fraction of drug 2 to exclude 0

%combination data for drug 1 and 2
surv12 = zeros(length(dose1),length(dose2)); %unaffected fraction in a matrix
surv12_ar = zeros(length(dose1)*length(dose2),3); %unaffected fraction in an array
k = 0; %index
for i = 1:length(dose1)
    for j = 1:length(dose2)
        k = k + 1;
        surv12(i,j) = data(intersect(find(data(:,1) == dose1(i)),...
            find(data(:,2) == dose2(j))),3)/100; %unaffected fraction in a matrix
        surv12_ar(k,:) = [dose1(i) dose2(j) surv12(i,j)]; %unaffected fraction in an array
    end 
end

%% 2 calculate Hill parameters for individual drugs
%first regression 
bot = min([surv1' surv2']); %initial guess of S0, assay background
c1 = package_envelope_hill([median(dose1) 2],[1;surv1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = c1(2); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 2],[1;surv2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = c2(2); %EC50 and Hill slope for drug 2

%second regression with adjusted S0, assay background
c = package_envelope_s0([lam1,h1,bot,lam2,h2],surv1,dose1,surv2,dose2); %regression for S0
bot = c(3); %new S0, assay background
c1 = package_envelope_hill([median(dose1) 2],[1;surv1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = (c1(2)); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 2],[1;surv2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = (c2(2)); %EC50 and Hill slope for drug 2

%% 3 calculate SAP (alpha)
alpha = package_alpha(ca,surv12_ar,lam1,h1,lam2,h2,bot);

%% 4 calculate beta and gamma from Bliss and Gaddum
[beta,gamma] = package_beta_gamma(cb,data,0);
