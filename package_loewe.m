function [f1,f2] = package_loewe(x1,x2,a,b,lam1,lam2,h1,h2)
f1 = b+(a-b)./(1.+(x1/lam1+(x2/lam2).^(h2/h1)).^h1);
f2 = b+(a-b)./(1.+(x2/lam2+(x1/lam1).^(h1/h2)).^h2);