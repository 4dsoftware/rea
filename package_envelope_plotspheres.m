function package_envelope_plotspheres(x,y,z,a,alpha,asp1,asp2,colstr)
[X,Y,Z] = ellipsoid(x,y,z,asp1*a,asp2*a,a);
hs = surf(X,Y,Z);
%set(hs,'cdata',ones(length(X)));
%colormap(hot);
%set(hs, 'FaceAlpha',alpha, 'EdgeAlpha', 0);
set(hs,'edgecolor','none','facecolor',colstr,'facealpha',alpha);
