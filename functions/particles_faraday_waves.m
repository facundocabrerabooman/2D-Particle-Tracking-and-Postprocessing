%% detect particles on faraday waves
clear all, close all, clear

name = '200Hz';

%pathin = ['/Users/FC/GoogleDrive/.shortcut-targets-by-id/1Pig6Q_rbK2uqGJyDkfT1vrOXI-BjZY3x/DDT_projects/Faraday_waves/particle_videos_ground/200Hz/' name filesep];
pathin = cd;
vids = dir([pathin filesep '*.MP4']);
cd(pathin)


%%% cut frames so to only see tank
%corners = [570 1300 1300 570; 115 115 845 845]; % upl-upr-lwr-lwl
a=1500; b=50; d=870; c=650; % 
corners = [c a a c; b b d d]; %upl-upr-lwr-lwl
%%%
pxtomm = 0.1428;

for ll=1:numel(vids) % loop over each video


    video = VideoReader([vids(ll).folder filesep vids(ll).name]);
    numFrames = video.numFrames;
        video_track = VideoWriter([pathin filesep name num2str(ll)], 'MPEG-4'); %create video
        video_track.FrameRate = 20;
        vide_track.Quality = 100;
        open(video_track);
    idx=0;
    CC=[];
    R=[];

    for ii = 1:numFrames % for over frames
        clear fr_cc fr_r
        ii/numFrames

        idx = idx+1;

        frame = read(video,ii);
        frame = rgb2gray(frame);
   %     framecut_or = frame(corners(2,2):corners(2,3),corners(1,1):corners(1,2));
        framecut_or = frame;
        framecut_or = imadjust(framecut_or);
        %
        framecut=framecut_or;
        framecut(framecut<240)=0;
        framecut_bw=imbinarize(framecut);

        %BW = imregionalmax(framecut,8);

         % size_strel = 1;
   %dd = strel('diamond',size_strel); 
   %BW = imopen(BW,dd);

    BW = framecut_bw;
        % figure(1);clf
        % subplot(1,2,1)
        % imshow(framecut_bw)
        % subplot(1,2,2);
        % imshow(afterOpening)

        %[cc, r] = imfindcircles(framecut,[1 10]);
        PP=regionprops(BW, 'Centroid', 'Eccentricity', 'EquivDiameter');
        PP([PP.Eccentricity]>.96 ,:) = [];
        %PP([PP.EquivDiameter]>9,:) = [];
        PP([PP.EquivDiameter]>9,:) = [];
        PP([PP.EquivDiameter]<3,:) = [];
        r = [PP.EquivDiameter].*2;
        cc = cat(1,PP.Centroid);

        %%% debug plot
%                 f=figure('visible', 'off');
%                  clf;imshow(frame);hold on
%                  viscircles(cc, r,'EdgeColor','r');
%                  title(num2str(ii))
%                  pause(0.1)
%                  writeVideo(video_track,getframe(f));
        %%%


        %cc = cc.*pxtomm; % to mm
        %r = r.*pxtomm;

        fr_cc = [ones(size(cc,1),1).*idx cc zeros(size(cc,1),1)]; % zeros are 3rd coordinate that we do not have
        fr_r = [ones(size(cc,1),1).*idx r' zeros(size(cc,1),1)];
        CC = [CC ; fr_cc];
        R = [R ; fr_r];

    end
    close(video_track)

    save([pathin filesep 'CC_R_' name '_' num2str(ll) '.mat'])

end

%% Track particles
name = '200Hz';

folderin = cd;
folderout = folderin;

traj_conc = [];

for oo = 1
    fname = [name '_' num2str(oo)];
    load(['CC_R_' fname],'CC')

    maxdist=10;
    lmin=40;
    flag_pred=1;
    npriormax=1;
    porder=3;
    flag_conf=1;
    nframes = 99999;
    [traj,tracks]=track2d_faradaywaves(CC,folderin,folderout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf,[],nframes);

    traj_conc=[traj_conc traj];
end

save('traj_conc','traj_conc')
%% Plot trajs 2D

load('traj_conc.mat')
traj=traj_conc;
pxtomm = 0.1428;

xt=vertcat(traj.x).*pxtomm;
yt=vertcat(traj.y).*pxtomm;

color=zeros(length(xt),1);
c=1;
for k=1:length(traj)
    temp=traj(k).length;
    color(c:c+temp-1)=k;
    c=c+temp;
end

%figure(4);clf
figure
[~,uniq] = unique(color);
for ii=1:numel(uniq)-1
    plot(xt(uniq(ii):uniq(ii+1)-1),yt(uniq(ii):uniq(ii+1)-1),'.-');hold on
end
axis equal
box on; grid on
xlabel('mm')
ylabel('mm')
%ylim([-5 95])
%xlim([-5 95])

savefig('traj_conc')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Traj Analysis

load('traj_conc.mat')

%%% compute velocity through a gaussian filter

w = 50; % filter intensity
fps = 500;
[~, trajf]=compute_vel_acc_traj(traj_conc,fps,w);

save('traj_concf','trajf','pxtomm')


%% Plot filtered traj

load('traj_concf.mat')
pxtomm = 0.1428;

xt=vertcat(trajf.xf).*pxtomm;
yt=vertcat(trajf.yf).*pxtomm;
color=zeros(length(xt),1);
c=1;
for k=1:length(trajf)
    temp=trajf(k).lengthf;
    color(c:c+temp-1)=k;
    c=c+temp;
end


%figure(4);clf
figure
[~,uniq] = unique(color);
for ii=1:numel(uniq)-1
    plot(xt(uniq(ii):uniq(ii+1)-1),yt(uniq(ii):uniq(ii+1)-1),'.-');hold on
end
axis equal
box on; grid on
xlabel('Pos. (mm)')
ylabel('Pos. (mm)')
%ylim([-5 95])
%xlim([-5 95])

fname = 'traj_concf';
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')


%% Histogram of position
load('traj_concf.mat')

figure(12);clf
subplot(1,2,1)
histogram(vertcat(trajf.xf).*pxtomm,50,'FaceColor','r');hold on
subplot(1,2,2)
histogram(vertcat(trajf.yf).*pxtomm,50,'FaceColor','g');

title('Hist. Positions [mm]'); %ylim([0 3000])

%xline(mean(vertcat(trajf.xf)),'r',mean(vertcat(trajf.xf)))
%xline(mean(vertcat(trajf.yf)),'g',mean(vertcat(trajf.yf)))

fname = 'hist_pos';
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')

%% Histogram vel, same fig
figure;clf
vx = vertcat(trajf.uf);
vy = vertcat(trajf.vf);

histogram((vx - mean(vx))./std(vx),30);hold on
histogram((vy - mean(vy))./std(vy),30);hold on

legend({'$V_x$','$V_y$'},'Interpreter', 'latex')
xlabel('$\frac{V_i - \langle V_i \rangle}{\sigma_{V_i}}$','Interpreter', 'latex')
grid;

fname = 'hist_vel_1fig';
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')

%% Histogram of velocity  + gaussian
figure(13);clf

subplot(1,2,1)
histfit(vertcat(trajf.uf).*pxtomm,30,'normal')
pdu = fitdist(vertcat(trajf.uf).*pxtomm,'Normal');
set(gca,'FontSize',17);
xlabel('Vel (mm/s)')
ylabel('Counts')
grid; title('$V_x$','Interpreter', 'latex')

subplot(1,2,2)
histfit(vertcat(trajf.vf).*pxtomm,30,'normal')
pdv = fitdist(vertcat(trajf.vf).*pxtomm,'Normal');
grid; title('$V_y$','Interpreter', 'latex')

xticks([-50:25:50])

disp('Mean velx [mm/s]')
mean(vertcat(trajf.uf))
disp('Mean vely')
mean(vertcat(trajf.vf))

fname = 'hist_vel';
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')


%% Histogram of velocity + gaussian -- turbulent way

% pdfVx = mkpdf6(trajf,'uf',65);
% pdfVy = mkpdf6(trajf,'vf',65);
% pdfAx = mkpdf6(trajf,'af',65);
% pdfAy = mkpdf6(trajf,'bf',65);
% 
% %%%
% figure;hold on;
%   semilogy(pdfVx.xpdfn,pdfAx.pdfn,'LineWidth',3)
%   semilogy(pdfVy.xpdfn,pdfAy.pdfn,'LineWidth',3)
%   semilogy(linspace(-5,5,128),gauss(linspace(-5,5,128),[1/sqrt(2*pi),0,1,0]),'LineWidth',3,'LineStyle','-.','Color','k');
%   legend('x','y','z','gauss','Location','South');
%   grid;
%   xlabel('v','Interpreter','latex');
%   ylabel('PDF','Interpreter','latex');
%   set(gca,'FontSize',24);
%   set(gca,'Xscale','linear','Yscale','log');
%   title('Velocity PDF');
%   
%   
%   stop
%   fname = 'pdf_vel';
%   savefig_FC(fname,8,6,'fig')
%   savefig_FC(fname,8,6,'pdf')


%% Histogram of acceleration -- FIT BY GAUSSSIAN (dataa should not be gaussian)
figure(14);clf

subplot(1,2,1)
histfit(vertcat(trajf.af).*pxtomm,30,'normal')
pdu = fitdist(vertcat(trajf.af).*pxtomm,'Normal');
set(gca,'FontSize',17);
xlabel('Vel (mm/s)')
ylabel('Counts')
grid; title('$A_x$','Interpreter', 'latex')

subplot(1,2,2)
histfit(vertcat(trajf.bf).*pxtomm,30,'normal')
pdv = fitdist(vertcat(trajf.bf).*pxtomm,'Normal');
grid; title('$A_y$','Interpreter', 'latex')

%xticks([-50:25:50])
%xlim([-1000,1000])
disp('Mean accx [mm/s]')
mean(vertcat(trajf.af))
disp('Mean accy')
mean(vertcat(trajf.bf))

fname = 'hist_acc';
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')
%% Pair dispersion simple 2D

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

        plot(trajf(oo).t_sec(I),D,'.-');hold on
        xlabel('t [s]'); ylabel('Sep [mm]')

    end

end

savefig('pairdisp')

%% compute lag stats

[Rvx,D2x,Nptsvx,Ntrackvx]=lagstats_tracks(trajf,'uf',2);
[Rvy,D2y,Nptsvy,Ntrackvy]=lagstats_tracks(trajf,'vf',2);

[Rax,~,Nptsax,Ntrackax]=lagstats_tracks(trajf,'af',2);
[Ray,~,Nptsay,Ntrackay]=lagstats_tracks(trajf,'bf',2);

%% correlations

%fps=58.9;
fps = 500;
t=(0:length(Rvx)-1)'/fps;

figure;
plot(t,Rvx./Rvx(1),'r-s'); hold on;
plot(t,Rvy./Rvy(1),'b-s'); hold on;
title('Rv')
fname = 'Rv';
xlim([0 1])
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')


figure;
plot(t,Rax./Rax(1),'r-s'); hold on;
plot(t,Ray./Ray(1),'b-s'); hold on;
title('Ra')
fname = 'Ra';
xlim([0 0.4])
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')

%% structure function
fps = 58.9;

t=(0:length(D2x)-1)'/fps;
D2x=D2x*(1e-3)^2; %from mm to m
D2y=D2y*(1e-3)^2; %from mm to m

% figure;
% loglog(t,D2x,'r-s');grid on
% set(gca,'fontsize',17) 
%   fname = 'D2Lx';
%   xlabel('$\tau$','Interpreter','latex')
%   ylabel('$D^2_L(\tau) = \langle [v_x(t+\tau) - v_x(t)]^2 \rangle$','Interpreter','latex')
% savefig_FC(fname,8,6,'fig')
% savefig_FC(fname,8,6,'pdf')
% 
% figure;
% loglog(t,D2y,'b-s');grid on
%   xlabel('$\tau$','Interpreter','latex')
%   ylabel('$D^2_L(\tau) = \langle [v_y(t+\tau) - v_y(t)]^2 \rangle$','Interpreter','latex')
%   fname = 'D2Ly';
% savefig_FC(fname,8,6,'fig')
% savefig_FC(fname,8,6,'pdf')


figure;
loglog(t,D2x,'r-s');grid on; hold on
loglog(t,D2y,'b-s');grid on
set(gca,'fontsize',17) 
  fname = 'D2L';
  xlabel('$\tau$','Interpreter','latex')
  ylabel('$D^2_L(\tau) = \langle [v_x(t+\tau) - v_x(t)]^2 \rangle$','Interpreter','latex')
savefig_FC(fname,8,6,'fig')
savefig_FC(fname,8,6,'pdf')
%% Fractal dimensions
load('traj_concf.mat')

frac_dim = boxcount_no_moisy(vertcat(trajf(:).xf), 10, 1)

%[n,r] = boxcount(ccc,'slope')

%%


















































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


%% Autocorrelations

load('traj_concf.mat')
d = vertcat(trajf.vf);

Rx = xcorr(vertcat(trajf.uf),'normalized');
figure;
plot(Rx)
title('autocorrelation uf')

Ry = xcorr(vertcat(trajf.vf),'normalized');
figure
plot(Ry)
title('autocorrelation vf')

Rposx = xcorr(vertcat(trajf.xf),'normalized');
figure
plot(Rposx)
title('autocorrelation xf')

Rposy = xcorr(vertcat(trajf.yf),'normalized');
figure
plot(Rposy)
title('autocorrelation yf')

stop
%%%%%%%%
Rx = xcorr(trajf(1).uf,trajf(2).uf);
figure;
plot(Rx)
title('crosssed correlation')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reshape for pair_dispersion2D.m -- did not work, OUTPUT file has sigle values when it should have vectors
% run pair_dispersion2D.m using this output

load('traj_concf.mat')

clearvars -except trajf
idx=0;
trajf_formatted=[];
for ii = 1:numel(trajf)

    idx=idx+1;

    trajf_formatted(idx).xf = trajf(ii).xf;
    trajf_formatted(idx).yf = trajf(ii).yf;
    trajf_formatted(idx).zf = trajf(ii).zf;
    trajf_formatted(idx).uf = trajf(ii).uf;
    trajf_formatted(idx).vf = trajf(ii).vf;
    trajf_formatted(idx).wf = trajf(ii).wf;
    trajf_formatted(idx).tf = trajf(ii).t_sec;


    trajf_formatted(idx).xv = trajf(ii).xf;
    trajf_formatted(idx).yv = trajf(ii).yf;
    trajf_formatted(idx).zv = trajf(ii).zf;
    trajf_formatted(idx).vx = trajf(ii).uf;
    trajf_formatted(idx).vy = trajf(ii).vf;
    trajf_formatted(idx).vz = trajf(ii).wf;
    trajf_formatted(idx).tf = trajf(ii).t_sec;
    trajf_formatted(idx).tv = trajf(ii).t_sec;
    trajf_formatted(idx).t = trajf(ii).t;
    trajf_formatted(idx).index = ii;
    trajf_formatted(idx).lengthf = numel(trajf(ii));
end

%stop
save('trajectories_to_pair_disp','trajf_formatted')




