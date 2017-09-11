function f = package_envelope_s0(c0,effu1,dose1,effu2,dose2)

function S = sub3(c)
y1 = (1-abs(c(3)))./((dose1/c(1)).^abs(c(2)) + 1) + abs(c(3));
y2 = (1-abs(c(3)))./((dose2/c(4)).^abs(c(5)) + 1) + abs(c(3));
S = ((norm(y1-effu1))^2+(norm(y2-effu2))^2);
end

options = optimset('TolX',1e-9,'TolFun',1e-9,'MaxIter',2e4,'MaxFunEval',2e4);
f = fminsearch(@sub3,c0,options);

end