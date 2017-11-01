function f = package_alpha(c0,data,lam1,h1,lam2,h2,bot)

dose1 = data(:,1);
dose2 = data(:,2);
surv = data(:,3);

%c = [lam1 h1 lam2 h2 bot alp]
function S = sub2(alp)
y = lam2*dose1.*((1-surv).^(1/h2)).*(surv-bot).^(1/h1) + lam1*dose2.*((1-surv).^(1/h1)).*(surv-bot).^(1/h2) +...
    alp*dose1.*dose2.*((1-surv).*(surv-bot)).^(1/2/h1 + 1/2/h2) - ...
    lam1*lam2*(1-surv).^(1/h1+1/h2);
S = norm(y);
end

options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxIter',1e4,'MaxFunEval',1e6,'Display','off');
f = fminsearch(@sub2,c0,options);

end