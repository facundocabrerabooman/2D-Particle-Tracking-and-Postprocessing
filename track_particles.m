%% Detect particles 
%%
clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/FC/Documents/GitHub/2D-Particle-Tracking-and-Postprocessing/'));

%%

clear all, close all, clear

name = 'Test 10';

pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/images 2'];
pathout = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/post processing 2'];
files = dir([pathin filesep '*.tiff']);
numFrames = numel(files)-20;

cd(pathout)


%%% cut frames so to only see tank
corners = [5 570 298 575]; %[top bottom left right]
%%%
pxtomm = 1;


video_track = VideoWriter([pathout filesep name],'MPEG-4'); %create video
video_track.FrameRate = 5;
vide_track.Quality = 20;
open(video_track);
idx=0;
CC=[];
R=[];


for ii = 1:1:11870 % for over frames
    clear fr_cc fr_r
    ii

    idx = idx+1;

    frame = imread([pathin filesep files(ii).name]);

    %%%%%
    framecut_or = frame(corners(1):corners(2), corners(3):corners(4));
    % imshow(framecut_or);

    frame = imadjust(framecut_or);

    % frame(frame<11000)=0;       %CHECK the intensity of the image
    % frame_bw=imbinarize(frame);
    % imshow(frame_bw);

    %if background is darker than particles
    frame_bw = ~frame;
    % imshow(frame_bw);


    %BW = imregionalmax(frame_bw,8);

    %[cc, r] = imfindcircles(framecut,[1 10]);
    PP=regionprops(frame_bw, 'Centroid', 'Eccentricity', 'EquivDiameter');


    % PP([PP.Eccentricity]>.95 ,:) = [];
    % PP([PP.EquivDiameter]<2,:) = [];
    PP([PP.EquivDiameter]<3,:) = [];      % CHECK take things out that are too big

    r = [PP.EquivDiameter].*2;
    cc = cat(1,PP.Centroid);

    %%% debug plot

    if mod(ii,50)==0
        f=figure('visible', 'off');
        clf;imshow(frame_bw);hold on
        viscircles(cc, r,'EdgeColor','r');
        title(num2str(ii))
        pause(0.1)
        writeVideo(video_track,getframe(gcf))
    end



    %cc = cc.*pxtomm; % to mm
    %r = r.*pxtomm;

    fr_cc = [ones(size(cc,1),1).*idx cc zeros(size(cc,1),1)]; % zeros are 3rd coordinate that we do not have
    fr_r = [ones(size(cc,1),1).*idx r' zeros(size(cc,1),1)];
    CC = [CC ; fr_cc];
    R = [R ; fr_r];

end
%close(video_track)

save([pathout filesep 'CC_R_' name '.mat'])


%% Track particles

folderin = cd;
folderout = folderin;

traj_conc = [];
name = 'Test ';

for oo = 10
    fname = [name num2str(oo)];
    load(['CC_R_' fname],'CC')

    maxdist=10;     %if particles are moving faster, then you need to increase this
    lmin=10;        %use this if particles disappear. if a particle is detected for more than 10 frames, consider it an track it


    flag_pred=1;
    npriormax=1;
    porder=3;
    flag_conf=1;
    nframes = 99999;
    [traj,tracks]=track2d(CC,folderin,folderout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf,[],nframes);

    traj_conc=[traj_conc traj];
end

save('traj_conc','traj_conc')
%% Plot trajs 2D

load('traj_conc.mat')
traj=traj_conc;
pxtomm = 1;

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
% ylim([-5 95])
% xlim([-5 95])

savefig('traj_conc')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Traj Analysis

load('traj_conc.mat')

%%% compute velocity through a gaussian filter

w = 3; % filter intensity
fps = 112;
[~, trajf]=compute_vel_acc_traj(traj_conc,fps,w);

save('traj_concf','trajf','pxtomm')


%% Plot filtered traj

load('traj_concf.mat')

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
xlabel('mm')
ylabel('mm')
% ylim([-5 95])
% xlim([-5 95])

savefig('traj_concf')


%% Histogram of position
load('traj_concf.mat')

figure(12);clf
subplot(1,2,1)
histogram(vertcat(trajf.xf).*pxtomm,50,'FaceColor','r');hold on
title('x Hist. Positions [mm]')
subplot(1,2,2)
histogram(vertcat(trajf.yf).*pxtomm,50,'FaceColor','g');

title('y Hist. Positions [mm]'); %ylim([0 3000])

%xline(mean(vertcat(trajf.xf)),'r',mean(vertcat(trajf.xf)))
%xline(mean(vertcat(trajf.yf)),'g',mean(vertcat(trajf.yf)))

savefig('hist_pos')


%% Histogram of velocity -- ADD A GAUSSIAN FIT ALSO NORMALIZE THIS
figure(13);clf

subplot(1,2,1)
histogram(vertcat(trajf.uf).*pxtomm,50,'FaceColor','r');hold on
subplot(1,2,2)
histogram(vertcat(trajf.vf).*pxtomm,50,'FaceColor','g');

title('Hist. Vel. [px/s]'); %xlim([-60 60]);% ylim([0 3000]);

%xline(mean(vertcat(trajf.uf)),'r',mean(vertcat(trajf.uf)))
%xline(mean(vertcat(trajf.vf)),'g',mean(vertcat(trajf.vf)))

disp('Mean velx [mm/s]')
mean(vertcat(trajf.uf))
subplot(1,2,1)
text(30,2000,['meanu = ' num2str(mean(vertcat(trajf.uf)))],'FontSize',12)
disp('Mean vely')
mean(vertcat(trajf.vf))
subplot(1,2,2)
text(30,2000,['meanv = ' num2str(mean(vertcat(trajf.vf)))],'FontSize',12)


savefig('hist_vel')
%% Velocity converge check
tot=numel(trajf);
numTraj=0;

for aa = 1:4
    
    numTraj = [tot,round(tot/2),round(tot/4),round(tot/8)];
    newTraj=trajf(1:numTraj(aa));
    meanu(aa)=mean(vertcat(newTraj.uf));
    meanv(aa)=mean(vertcat(newTraj.vf));

end 
    
figure()
plot(numTraj,meanu,'blue')
hold on
plot(numTraj,meanv,'red')

%% Histogram of acceleration -- FIT BY GAUSSSIAN
figure(14);clf

subplot(1,2,1)
histogram(vertcat(trajf.af),'FaceColor','r');hold on
legend({'af','bf'});title('Hist. Acc.'); xlim([-1e4 1e4]); %ylim([0 8000]);

subplot(1,2,2)
histogram(vertcat(trajf.bf),'FaceColor','g');
legend({'af','bf'});title('Hist. Acc.'); xlim([-1e4 1e4]); %ylim([0 8000]);


disp('Mean accx [mm/s2]')
mean(vertcat(trajf.af))
subplot(1,2,1)
text(30e2,1000,['mean = ' num2str(mean(vertcat(trajf.af)))],'FontSize',12)
disp('Mean accy')
mean(vertcat(trajf.bf))
subplot(1,2,2)
text(30e2,1000,['mean = ' num2str(mean(vertcat(trajf.bf)))],'FontSize',12)

savefig('hist_acc')

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

%% compute lag stats

[Rvx,D2x,Nptsvx,Ntrackvx]=lagstats_tracks(trajf,'uf',2);
[Rvy,D2y,Nptsvy,Ntrackvy]=lagstats_tracks(trajf,'vf',2);

[Rax,~,Nptsax,Ntrackax]=lagstats_tracks(trajf,'af',2);
[Ray,~,Nptsay,Ntrackay]=lagstats_tracks(trajf,'bf',2);

%% correlations RUN LAG STATS IF YOU WANT THIS STUFF OR IF YOU CHANGED THINGS AROUND

fps = 112;
t=(0:length(Rvx)-1)'/fps;

figure;
plot(t,Rvx./Rvx(1),'r-s'); hold on;
plot(t,Rvy./Rvy(1),'b-s'); hold on;
title('Rv')
fname = 'Rv';
xlim([0 0.1])
savefig('Rv')

figure;
plot(t,Rax./Rax(1),'r-s'); hold on;
plot(t,Ray./Ray(1),'b-s'); hold on;
title('Ra')
fname = 'Ra';
xlim([0 0.1])
savefig('Ra')

%% structure function lagrangian (this changes scale each time I rn without clearing??)
fps = 112;

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
xlim([0.003 0.6])
savefig('d2L')

%% Fractal dimensions
load('traj_concf.mat')

frac_dim = boxcount_no_moisy(coordmat, 10, 1)

%[n,r] = boxcount(ccc,'slope')


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
    trajf_formatted(idx).lengthf = trajf(ii).lengthf;
end

%stop
save('trajectories_to_pair_disp','trajf_formatted')




