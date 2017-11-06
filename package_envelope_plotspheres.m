function package_envelope_plotspheres(x,y,z,a,alpha,asp1,asp2,colstr)
[X,Y,Z] = ellipsoid(x,y,z,asp1*a,asp2*a,a);
hs = surf(X,Y,Z);
set(hs,'edgecolor','none','facecolor',colstr,'facealpha',alpha);
