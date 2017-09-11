function void = package_envelope_plotter(lam1,lam2,a,b,h1,h2,range11,range12,range21,range22,grids,ft)

% Response Envelope
% Created by Daniel Du, 4/1/2016
% This program is used to study drug combination effect using response
% envelope. 

x1 = 10.^(linspace(range11,range12,grids));
x2 = 10.^(linspace(range21,range22,grids));
f1 = zeros(grids);
f2 = zeros(grids);
f3 = zeros(grids);
fupp = zeros(grids);
fbot = zeros(grids);

for i = 1:grids
    for j = 1:grids
        [f1(i,j),f2(i,j)] = package_loewe(x1(i),x2(j),a,b,lam1,lam2,h1,h2);
        f3(i,j) = package_bliss(x1(i),x2(j),a,b,lam1,lam2,h1,h2);
        fupp(i,j) = max([f1(i,j) f2(i,j) f3(i,j)]);
        fbot(i,j) = min([f1(i,j) f2(i,j) f3(i,j)]);
    end
end
%set 1
h1 = surf(log10(x1),log10(x2),fupp');
set(h1,'edgecolor','none','facecolor',0.7*[1 1 1],'facealpha',0.4);
hold on
h2 = surf(log10(x1),log10(x2),fbot');
set(h2,'edgecolor','none','facecolor',0.3*[1 1 1],'facealpha',0.4);
hold off

xlabel('Dose 1','fontsize',ft);
ylabel('Dose 2','fontsize',ft);
zlabel('S','fontsize',ft);
set(gca,'fontsize',ft);
view(130,30);

return
