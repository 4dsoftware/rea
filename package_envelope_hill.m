function f = package_envelope_hill(c0,surv,dose,bot)
function S = sub2(c)
y = (1-bot)./((dose/c(1)).^(c(2)) + 1) + bot; 
S = ((norm(y-surv))^2);
end
options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxIter',1e4,'MaxFunEval',1e4,'Display','off');
f = fmincon(@sub2,c0,[],[],[],[],[0 0],[Inf 4],[],options); %Hill coefficients and EC50's cannot be negative
end