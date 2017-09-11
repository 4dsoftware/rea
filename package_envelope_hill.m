function f = package_envelope_hill(c0,effu,dose,bot)

function S = sub2(c)
y = (1-bot)./((dose/c(1)).^abs(c(2)) + 1) + bot; %h1 and h2 cannot be negative
S = ((norm(y-effu))^2);
end

options = optimset('TolX',1e-8,'TolFun',1e-8,'MaxIter',5e3,'MaxFunEval',5e3);
f = fminsearch(@sub2,c0,options);

end