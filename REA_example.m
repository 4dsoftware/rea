%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response Envelope Analysis (REA)
%
% Created by Daniel Du, Jan 2017
% Proteomics and Metabolimics Core Facility
% MD Anderson Cancer Center, Houston, TX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note

%This script performs response envelope analysis on a published data set. 
%The data set can be found in Griner PNAS 2014. Response envelope analysis 
%(REA) is a tool to quantitatively determine combination effects including
%synergy, additivity, and antagonism. 

%Griner, L.A.M. et al, PNAS 2014. doi: 10.1073/pnas.1311846111

%--------------------------------------------------------------------------

%% 0 parameter setting
clear
clc
mksz = 0.03; %marker size 
trim = 0.1; %trim the x scale and y scale this amount beyond the data points for 2-D representation
ft = 9; %font size
lw = 3; %line width


%% 1 load the data
%read data file
data = csvread('Griner_PNAS_2014_1.csv'); 
%load the data set into the environment, the data are formatted as the
%following:
%column 1, concentrations of drug 1
%column 2, concentrations of drug 2
%column 3, corresponding survival rates
%if the data show affected ratio or drug response instead, then the 
%following step is needed
%data(:,3) = 100-data(:,3); 

%single drug response curve for drug 1
dose1 = sort(unique(data(:,1))); %concentrations of drug 1 
dose1 = dose1(2:end); %truncated concentrations of drug 1 to exclude 0
base1 = dose1(1)^2/dose1(3); 
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
effu1 = data((data(:,2) == 0),3)/100; %unaffected fraction of drug 1
effu1 = effu1(2:end); %truncated unaffected fraction of drug 1 to exclude 0

%single drug response curve for drug 2
dose2 = sort(unique(data(:,2))); %concentrations of drug 2 
dose2 = dose2(2:end); %truncated concentrations of drug 2 to exclude 0
base2 = dose2(1)^2/dose2(3);
%on logarithmic scale, there can't be zeros. Here we set base1 to be zero,
%which is way below the lowest experiment concentration on this axis. 
effu2 = data((data(:,1) == 0),3)/100; %unaffected fraction of drug 2
effu2 = effu2(2:end); %truncated unaffected fraction of drug 2 to exclude 0

%combination data for drug 1 and 2
effu12 = zeros(length(dose1),length(dose2)); %unaffected fraction in a matrix
effu12_ar = zeros(length(dose1)*length(dose2),3); %unaffected fraction in an array
k = 0; %index
for i = 1:length(dose1)
    for j = 1:length(dose2)
        k = k + 1;
        effu12(i,j) = data(intersect(find(data(:,1) == dose1(i)),...
            find(data(:,2) == dose2(j))),3)/100; %unaffected fraction in a matrix
        effu12_ar(k,:) = [dose1(i) dose2(j) effu12(i,j)]; %unaffected fraction in an array
    end 
end

%% 2 calculate Hill parameters for individual drugs
%first regression 
bot = min([effu1' effu2']); %initial guess of S0, assay background
c1 = package_envelope_hill([median(dose1) 0.5],[1;effu1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = c1(2); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 0.5],[1;effu2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = c2(2); %EC50 and Hill slope for drug 2

%second regression with adjusted S0, assay background
c = package_envelope_s0([lam1,h1,bot,lam2,h2],effu1,dose1,effu2,dose2); %regression for S0
bot = c(3); %new S0, assay background
c1 = package_envelope_hill([median(dose1) 0.5],[1;effu1],[base1;dose1],bot); %regression for drug 1
lam1 = c1(1); h1 = c1(2); %EC50 and Hill slope for drug 1
c2 = package_envelope_hill([median(dose2) 0.5],[1;effu2],[base2;dose2],bot); %regression for drug 2
lam2 = c2(1); h2 = c2(2); %EC50 and Hill slope for drug 2

%% 3 visualization
%initialize figure object
figure('position',[400 50 1500 310]);

%---------------------------  PANEL 1  ------------------------------------
%This panel provides visualization of the data and the envelope in 3-D
subplot(141);

%envelope plotter
package_envelope_plotter(lam1,lam2,1,bot,h1,h2,log10(base1)-0.3,...
    log10(max(dose1))+0.3,log10(base2)-0.3,log10(max(dose2))+0.3,200,ft); 

%define 3-d markers
asp1 = (log10(max(dose1))-log10(base1))/1.2; %define aspect ratio for axis 1
asp2 = (log10(max(dose2))-log10(base2))/1.2; %define aspect ratio for axis 2
alpha = 0.7; %transparency

%single drug response for drug 1
for i = 1:length(dose1)
        hold on
        package_envelope_plotspheres(log10(dose1(i)),log10(base2),effu1(i),mksz,alpha,asp1,asp2,[0 0 0]);
end

%single drug response for drug 2
for j = 1:length(dose2)
        hold on
        package_envelope_plotspheres(log10(base1),log10(dose2(j)),effu2(j),mksz,alpha,asp1,asp2,[0 0 0]);
end

%combination data for drug 1 and 2
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2); %the effect calculated from generalized Loewe
    f3 = package_bliss(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2); %the effect calculated from Bliss
    if effu12_ar(k,3) < min([f1 f2 f3]) %red for synergy
        hold on
        package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),effu12_ar(k,3),mksz,alpha,asp1,asp2,[1 0 0]);
    elseif effu12_ar(k,3) > max([f1 f2 f3]) %blue for antagonism
        hold on
        package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),effu12_ar(k,3),mksz,alpha,asp1,asp2,[0 0 1]);
    else %white for additivity
        hold on 
        package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),effu12_ar(k,3),mksz,alpha,asp1,asp2,[1 1 1]);
    end
end

%notations
xlim([log10(base1)-mksz*asp1 log10(max(dose1))+mksz*asp1]); %x axis limit
ylim([log10(base2)-mksz*asp2 log10(max(dose2))+mksz*asp2]); %y axis limit
xlabel('ibrutinib (\muM)','fontsize',ft); %x axis label
ylabel('MK-2206 (\muM)','fontsize',ft); %y axis label
zlim([0 1.2]); %z axis limit
camlight right; %define camera light
lighting phong; %define light type
title('ibrutinib and MK-2206','fontsize',ft); %title of the figure

%replace logarithmic tick labels with algebraic ones
xlb = 2.^(-1:2:6);
ylb = 2.^(3:2:12);
set(gca,'fontsize',ft,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb))...
    ,'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
view(150,20); %camera angle

%---------------------------  PANEL 2  ------------------------------------
%This panel provides 2-D projection of the combination data shown in panel
%1, thus no comments are shown for similar scripts. 
subplot(142);

mksz = 12; %redefine marker size for this panel
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    if effu12_ar(k,3) < min([f1 f2 f3]) %red for synergy
        hold on
        plot(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),'rs','markersize',mksz,'markerfacecolor','r');
    elseif effu12_ar(k,3) > max([f1 f2 f3]) %blue for antagonism
        hold on
        plot(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),'bs','markersize',mksz,'markerfacecolor','b');
    end
end
box on
xlim([log10(min(dose1))-trim log10(max(dose1))+trim]);
ylim([log10(min(dose2))-trim log10(max(dose2))+trim]);
xlabel('ibrutinib (\muM)','fontsize',ft);
ylabel('MK-2206 (\muM)','fontsize',ft);
xlb = 2.^(-1:2:7);
ylb = 2.^(3:2:12);
set(gca,'fontsize',ft,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb))...
    ,'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
set(gca,'fontsize',ft);

%---------------------------  PANEL 3  ------------------------------------
%This panel identifies the largest islands of synergy and antagonism using 
%connected-component labeling. 
subplot(143);

%label the data points
ff = zeros(1,length(dose1)*length(dose2)); %synergy labels
gg = zeros(1,length(dose1)*length(dose2)); %antagonism labels
for k = 1:(length(dose1)*length(dose2))
    [f1,f2] = package_loewe(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    if effu12_ar(k,3) < min([f1 f2 f3])
        ff(k) = 1;
    end
    if effu12_ar(k,3) > max([f1 f2 f3])
        gg(k) = 1;
    end
end

%find the largest island of synergy
ff = reshape(ff,9,9); %reshape synergy label array to matrix
CC = bwconncomp(ff,4); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,id] = max(numPixels);
reg = CC.PixelIdxList{id}'; %data list corresponding to the largest island of synergy 

%find the largest island of antagonism
gg = reshape(gg,9,9); %reshape antagonism label array to matrix
CC = bwconncomp(gg,4); %connected-component labeling
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest2,id] = max(numPixels);
reg2 = CC.PixelIdxList{id}'; %data list corresponding to the largest island of antagonism 

%visualize the island
plot(log10(effu12_ar(reg,1)),log10(effu12_ar(reg,2)),'rs','markersize',mksz,'markerfacecolor','r');
hold on
plot(log10(effu12_ar(reg2,1)),log10(effu12_ar(reg2,2)),'bs','markersize',mksz,'markerfacecolor','b');
hold off
xlim([log10(min(dose1))-trim log10(max(dose1))+trim]);
ylim([log10(min(dose2))-trim log10(max(dose2))+trim]);
xlabel('ibrutinib (\muM)','fontsize',ft);
ylabel('MK-2206 (\muM)','fontsize',ft);
xlb = 2.^(-1:2:7);
ylb = 2.^(3:2:12);
set(gca,'fontsize',ft,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb))...
    ,'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
set(gca,'fontsize',ft);

%---------------------------  PANEL 4  ------------------------------------
%This panel visualizes the difference of the measured data from the 
%envelope for synergy and antagonism islands and calculates the synergy 
%index (SI) and antagonism index (AI). 
subplot(144);

%envelope plotter
mksz = 0.015; %redefine marker size for this panel
package_envelope_plotter(lam1,lam2,1,bot,h1,h2,log10(base1)-0.3,...
    log10(max(dose1))+0.3,log10(base2)-0.3,log10(max(dose2))+0.3,200,ft);
hold on

%difference of the measured data from the envelope for synergy
nn = 0;
dif = zeros(1,length(reg));
for k = reg
    nn = nn + 1;
    package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),effu12_ar(k,3),mksz,alpha,asp1,asp2,[1 0 0]);
    hold on
    [f1,f2] = package_loewe(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    flow = min([f1 f2 f3]);
    plot3([log10(effu12_ar(k,1)) log10(effu12_ar(k,1))],[log10(effu12_ar(k,2)) log10(effu12_ar(k,2))],...
        [effu12_ar(k,3) flow],'linewidth',lw,'color','r')
    hold on
    package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),flow,mksz,alpha,asp1,asp2,[1 0 0]);
    hold on
    dif(nn) = abs(effu12_ar(k,3) - flow);
end

%difference of the measured data from the envelope for antagonism
nn = 0;
dif2= zeros(1,length(reg2));
for k = reg2
    nn = nn + 1;
    package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),effu12_ar(k,3),mksz,alpha,asp1,asp2,[0 0 1]);
    hold on
    [f1,f2] = package_loewe(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    f3 = package_bliss(effu12_ar(k,1),effu12_ar(k,2),1,bot,lam1,lam2,h1,h2);
    flow = max([f1 f2 f3]);
    plot3([log10(effu12_ar(k,1)) log10(effu12_ar(k,1))],[log10(effu12_ar(k,2)) log10(effu12_ar(k,2))],...
        [effu12_ar(k,3) flow],'linewidth',lw,'color','b');
    hold on
    package_envelope_plotspheres(log10(effu12_ar(k,1)),log10(effu12_ar(k,2)),flow,mksz,alpha,asp1,asp2,[0 0 1]);
    hold on
    dif2(nn) = abs(effu12_ar(k,3) - flow);
end
hold off

si = sum(dif)/size(data,1); %SI
ai = sum(dif2)/size(data,1); %AI

xlim([log10(base1)-mksz*asp1 log10(max(dose1))+mksz*asp1]);
ylim([log10(base2)-mksz*asp2 log10(max(dose2))+mksz*asp2]);
xlabel('ibrutinib (\muM)','fontsize',ft);
ylabel('MK-2206 (\muM)','fontsize',ft);
zlim([0 1.2]);
camlight right; 
lighting phong;
title(['SI = ',num2str(si,2),', AI = ',num2str(ai,2)],'fontsize',ft);
xlb = 2.^(-1:2:6);
ylb = 2.^(3:2:12);
set(gca,'fontsize',ft,'xtick',(log10(xlb)),'xticklabel',num2cell((xlb))...
    ,'ytick',(log10(ylb)),'yticklabel',num2cell((ylb)));
set(gca,'fontsize',ft);
view(150,20);

