function [f,g] = package_beta_gamma(c0,data,s0)
% convert array data into matrix data
dose1 = sort(unique((data(:,1))));
dose2 = sort(unique((data(:,2))));
for i = 1:length(dose1)
    data(data(:,1) == dose1(i),1) = i;
end
for j = 1:length(dose2)
    data(data(:,2) == dose2(j),2) = j;
end
surv = zeros(length(dose1),length(dose2));
for k = 1:size(data,1)
    surv(data(k,1),data(k,2)) = data(k,3)/100;
end
fu = (surv-s0)/(1-s0); %uneffected ratio

% non-linear regression setting
options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxIter',1e4,'MaxFunEval',1e6,'Display','off');

% calculate beta based on Bliss independence
function S = sub2(beta)
S = 0;
for ii = 1:length(dose1) 
    for jj = 1:length(dose2)
        S = S + (fu(ii,jj) - beta*fu(ii,1)*fu(1,jj)).^2;
    end
end
end
f = fminsearch(@sub2,c0,options);  %beta

% calculate gamma based on Gaddum's non-interaction model
function S = sub3(gamma)
S = 0;
for ii = 1:length(dose1)
    for jj = 1:length(dose2)
        S = S + (fu(ii,jj) - gamma*min((fu(ii,1)),(fu(1,jj)))).^2; %If the unaffected ratio is considered, gamma < 1 gives synergy
    end
end
end
g = fminsearch(@sub3,c0,options);  %gamma

end