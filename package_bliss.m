function f = package_bliss(x1,x2,a,b,lam1,lam2,h1,h2)
f = b+(a-b)./(1.+(x1/lam1).^h1)./(1.+(x2/lam2).^h2);