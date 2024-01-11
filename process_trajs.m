%% Overall Definitions

clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/2D-Particle-Tracking-and-Postprocessing/'));

Fs = 112;

fname = 'EXP11';

% pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
%     '/ISS-CASIS/Experimental Design/Earth tests/' name '/images'];
% pathout = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
%     '/ISS-CASIS/Experimental Design/Earth tests/' name '/post processing'];

pathin = '/Volumes/landau2/ISS/EXP11/images/';
pathout = '/Volumes/landau2/ISS/EXP11/postproc/'; mkdir(pathout)
cd(pathout);

load('trajf','trajf')
trajs_conc=trajf;
%% Nice colors
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
%% Compute pdfs

pdfX(1) = mkpdf5(trajs_conc,'Xf',256,10);
pdfX(2) = mkpdf5(trajs_conc,'Yf',256,10);

pdfV(1) = mkpdf5(trajs_conc,'Vx',256,10);
pdfV(2) = mkpdf5(trajs_conc,'Vy',256,10);

pdfA(1) = mkpdf5(trajs_conc,'Ax',256,20);
pdfA(2) = mkpdf5(trajs_conc,'Ay',256,20);


%% Plot Normalized PDFs

figure;
semilogy(pdfX(1).xpdfn,pdfX(1).pdfn,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on;
semilogy(pdfX(2).xpdfn,pdfX(2).pdfn,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);

%set(gca,FontSize=15)
legend('$X$','$Y$','$Z$','interpreter','latex',Location='northeast');
%title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(\frac{X-\langle X \rangle}{std(X)})$','interpreter','latex',FontWeight='bold')
xlabel('$\frac{X-\langle X \rangle}{std(X)}$','interpreter','latex',FontWeight='bold')
grid off
axis padded
ylim([5e-7 1])

folderout = 'pdfs/';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_x'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_x'],8,6,'fig')
%%%%%%%%%


figure;
semilogy(pdfV(1).xpdfn,pdfV(1).pdfn,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on;
semilogy(pdfV(2).xpdfn,pdfV(2).pdfn,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

legend('$V_x$','$V_y$','Gaussian','interpreter','latex',Location='northeast');
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('$\frac{V-\langle V \rangle}{std(V)}$','interpreter','latex',FontWeight='bold')
grid off
axis padded
ylim([5e-7 1])

folderout = 'pdfs/';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_v'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_v'],8,6,'fig')
%%%%%%%%%

figure;

semilogy(pdfA(1).xpdfn,pdfA(1).pdfn,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2); hold on
semilogy(pdfA(2).xpdfn,pdfA(2).pdfn,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);

xpdf=linspace(-5,5,1024);
plot(xpdf,normpdf(xpdf,0,1),Color=color1,LineWidth=3);

ylabel('$PDF(\frac{A-\langle A \rangle}{std(A)})$','interpreter','latex',FontWeight='bold')
legend('$A_x$','$A_y$','Gaussian','interpreter','latex',Location='northeast');
xlabel('$\frac{A-\langle A \rangle}{std(A)}$','interpreter','latex',FontWeight='bold')
grid off
%xlim([-5 5])
ylim([5e-7 1])


folderout = 'pdfs/';
mkdir(folderout)
savefig_FC([folderout 'PDF_a'],8,6,'pdf')
savefig_FC([folderout 'PDF_a'],8,6,'fig')

%% Table with moments of distribution
maketable(pdfA,pdfV,folderout)

%% Mean Square Separation
MSD(1) = structFunc_struct(trajs_conc,'Xf',2);
MSD(2) = structFunc_struct(trajs_conc,'Yf',2);

%% Plot
figure;
loglog(MSD(1).tau/Fs,MSD(1).mean,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(MSD(2).tau/Fs,MSD(2).mean,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);

xMSD = linspace(1,350,1000)/Fs;
loglog(xMSD,0.5e5*xMSD.^2,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$MSD^x$','$MSD^y$','Interpreter','latex', 'Location','southeast')
%title('$MSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$MSD(m^2)$','interpreter','latex',FontWeight='bold')
xlabel('$\tau(s)$','interpreter','latex',FontWeight='bold')
text(2e-3,3,'$\tau^2$','interpreter','latex',FontWeight='bold',FontSize=20)
grid on
axis padded


folderout = 'MSS/';
mkdir(folderout)
savefig_FC([folderout 'MSS'],8,6,'pdf')
savefig_FC([folderout 'MSS'],8,6,'fig')

%% Longitudinal S2

S2L(1)= structFunc_struct(trajs_conc,'Vx',2);
S2L(2)= structFunc_struct(trajs_conc,'Vy',2);

%% Plot S2L 
figure;
loglog(S2L(1).tau/Fs,S2L(1).mean,'d',MarkerSize=3,Color=color3(1,:),LineWidth=2);hold on
loglog(S2L(2).tau/Fs,S2L(2).mean,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);

xS2L = linspace(1,15,100)/Fs;
loglog(xS2L,2e8*xS2L.^2,'--',Color=color1,LineWidth=2)
xS2L = linspace(16,200,100)/Fs;
loglog(xS2L,9.5e5*xS2L.^1,'--',Color=color1,LineWidth=2)


legend('$S_2^L(x)$','$S_2^L(y)$','interpreter','latex',Location='southeast')
ylabel('$S_2^L$','interpreter','latex',FontWeight='bold')
xlabel('$\tau / sec$','interpreter','latex',FontWeight='bold')
text(8e-4,8e2,'$\tau^2$','interpreter','latex','FontWeight','bold','FontSize',20)
text(1e-2,1.5e4,'$\tau$','interpreter','latex','FontWeight','bold','FontSize',20)
grid on
axis padded

folderout = 'S2L/';
mkdir(folderout)
savefig_FC([folderout 'S2L'],8,6,'pdf')
savefig_FC([folderout 'S2L'],8,6,'fig')


%% Velocity and Acceleration Correlations

Ruu(1) = xcorr_struct(trajs_conc,'Vx',1);
Ruu(2) = xcorr_struct(trajs_conc,'Vy',1);

Raa(1) = xcorr_struct(trajs_conc,'Ax',1);
Raa(2) = xcorr_struct(trajs_conc,'Ay',1);

%% Fit correlation 

nCorrFitV = 270;
nCorrFitA = 20;

if2layersV = 2;
if2layersA = 99;

ifboundedV = [1 1 1]; 
ifboundedA = [1 1 1];

Ruufit(1) = correlationFit(Ruu(1),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(1));
Ruufit(2) = correlationFit(Ruu(2),Fs,1,nCorrFitV,'V',if2layersV,ifboundedV(2));
Raafit(1) = correlationFit(Raa(1),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(1));
Raafit(2) = correlationFit(Raa(2),Fs,1,nCorrFitA,'A',if2layersA,ifboundedA(2));
clear nCorrFitV nCorrFitA
clear if2layersV  if2layersA
clear ifboundedV ifboundedA

%% Plot Correlation function and fit

f1 = figure;
tiledlayout(2,1)
nexttile
plot(Ruu(1).tau/Fs,Ruu(1).mean/Ruu(1).mean(1),'d',Color=color3(1,:),MarkerSize=1);hold on
plot(Ruu(2).tau/Fs,Ruu(2).mean/Ruu(2).mean(1),'d',Color=color3(2,:),MarkerSize=1);
plot(Ruufit(1).x,Ruufit(1).yfit,'-',Color=color3(1,:));hold on
plot(Ruufit(2).x,Ruufit(2).yfit,'-',Color=color3(2,:))

legend('$x$','$y$',Location='eastoutside')
title('$R_{uu}$')
ylabel('$\frac{\langle u(t)u(t+\tau) \rangle} {\langle u^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

nexttile
plot(Raa(1).tau/Fs,Raa(1).mean/Raa(1).mean(1),'o',Color=color3(1,:),MarkerSize=1);hold on
plot(Raa(2).tau/Fs,Raa(2).mean/Raa(2).mean(1),'o',Color=color3(2,:),MarkerSize=1);
plot(Raafit(1).x,Raafit(1).yfit,'-',Color=color3(1,:))
plot(Raafit(2).x,Raafit(2).yfit,'-',Color=color3(2,:))

legend('$x$','$y$',Location='eastoutside')
title('$R_{aa}$')
ylabel('$\frac{\langle a(t)a(t+\tau) \rangle} {\langle a^2(t) \rangle}$');
xlabel('$\tau$/s')
grid on
axis tight

folderout = 'corr/';
mkdir(folderout)
savefig_FC([folderout 'corr_classic'],8,6,'pdf')
savefig_FC([folderout 'corr_classic'],8,6,'fig')

stop
%% Fractal dimensions

frac_dim = boxcount_no_moisy(trajs_conc, 10, 1)


%% Fourier transforms

%%%% ftt

load('traj_concf.mat')
d = vertcat(trajf.vf);
fps = 60;

L = size(d,1);
T = 1/fps;
Fs = fps;

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
dfft = fft(detrend(d),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

%%% find peaks of TFs
[peakValues_perp, indexesOfPeaks_perp] = findpeaks(2*abs(dfft(1:NFFT/2+1)), 'MINPEAKHEIGHT', max(2*abs(dfft(1:NFFT/2+1))*0.6));

% Plot single-sided amplitude spectrum.
figure(20);clf

plot(f,2*abs(dfft(1:NFFT/2+1)),'-.') , hold on, plot(f(indexesOfPeaks_perp),peakValues_perp,'r*')
title('Vel perp - Single-Sided Amplitude Spectrum of d(t)')
xlabel('Frequency (Hz)')
ylabel('|fft|'), grid on

%%% Pwelch
figure(21);clf

L = size(d,1); %length of signal
N = nextpow2(L)-1; % Next power of 2 from length of y
[p_vr,f_vr] = pwelch(d,hamming(2^(N)),2^(N-1),2^(N),fps);
semilogy(f_vr,p_vr), hold on
[fitresult, ~] = gaussian_on_pwelch(f_vr,p_vr, 1);


%% Pair dispersion simple 2D %this takes a while

figure(15);
load('traj_concf.mat')

Dconc = [];
Dall =[];
for oo = 1:numel(trajf)-1

    for ll = 1:numel(trajf)-oo

        if length(trajf(oo).t) > length(trajf(oo+ll).t)
            diffe = trajf(oo).t - [trajf(oo+ll).t ; ones(length(trajf(oo).t) - length(trajf(oo+ll).t),1)];
        else
            diffe = trajf(oo+ll).t - [trajf(oo).t ; ones(length(trajf(oo+ll).t) - length(trajf(oo).t),1)];
        end

        I = find(diffe==0);

        D = pdist2([trajf(oo).xf(I) trajf(oo).yf(I)].*pxtomm,[trajf(oo+ll).xf(I) trajf(oo+ll).yf(I)].*pxtomm,'euclidean');
        D = diag(D)-D(1,1);
        Dall{oo,ll} = D;

        Dconc = [D ; Dconc];

        plot(trajf(oo).t_sec(I),D,'.');hold on
        xlabel('t [s]'); ylabel('Sep [mm]')

    end

end

savefig('pairdisp')
