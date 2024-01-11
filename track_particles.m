%% Overall Definitions

clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/2D-Particle-Tracking-and-Postprocessing/'));


pxtomm = 0.091;
fps = 150; 

fname = 'EXP11_a';

% pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
%     '/ISS-CASIS/Experimental Design/Earth tests/' name '/images'];
% pathout = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
%     '/ISS-CASIS/Experimental Design/Earth tests/' name '/post processing'];

pathin = '/Volumes/landau2/ISS/EXP11/images/';
pathout = '/Volumes/landau2/ISS/EXP11/postproc/'; mkdir(pathout)
cd(pathout);

%% Detect particles -- choose intensity_cut

%%% cut frames so to only see tank
corners = [126 590 255 715]; %[top bottom left right]
intensity_cut = 4.5e4;
%%%

files = dir([pathin filesep '*.tiff']);
numFrames = numel(files)-1;

video_track = VideoWriter([pathout filesep fname],'MPEG-4'); %create video
video_track.FrameRate = 5;
vide_track.Quality = 20;
open(video_track);
idx=0;
CC=[];
R=[];


for ii = 1:numFrames
    clear fr_cc fr_r
    ii

    idx = idx+1;

    frame = imread([pathin filesep files(ii).name]);
    %frame_subs = imsubstract(bkg,frame);
    %%%%%
    frame_cut = frame(corners(1):corners(2), corners(3):corners(4));
    % imshow(framecut_or);

    frame_cut = imadjust(frame_cut);

    frame_cut(frame_cut<intensity_cut)=0;       %CHECK the intensity of the image
    frame_bw=imbinarize(frame_cut);

    PP=regionprops(frame_bw, 'Centroid', 'Eccentricity', 'MajorAxisLength');

    PP([PP.MajorAxisLength]>7,:) = [];      % Take out things that are too big

    r = [PP.MajorAxisLength]./2;
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



    cc = cc.*pxtomm; % to mm
    r = r.*pxtomm;

    fr_cc = [ones(size(cc,1),1).*idx cc zeros(size(cc,1),1)]; % zeros are 3rd coordinate that we do not have
    fr_r = [ones(size(cc,1),1).*idx r' zeros(size(cc,1),1)];
    CC = [CC ; fr_cc];
    R = [R ; fr_r];

end

close(video_track)

save([pathout filesep 'CC_R_' fname '.mat'])


%% Track particles

pathin = cd;
pathout = pathin;

traj_conc = [];

load(['CC_R_' fname],'CC')

maxdist=10;     %if particles are moving faster, then you need to increase this
lmin=10;        %use this if particles disappear. if a particle is detected for more than 10 frames, consider it and track it
flag_pred=1;
npriormax=1;
porder=3;
flag_conf=1;
nframes = 99999;
[traj,tracks]=track2d(CC,pathin,pathout,fname,maxdist,lmin,flag_pred,npriormax,porder,flag_conf,[],nframes);

traj_conc=[traj_conc traj];

save('traj_conc','traj_conc')
%% Plot Trajs

%load('traj_conc.mat')

traj=traj_conc;

xt=vertcat(traj.x).*pxtomm;
yt=vertcat(traj.y).*pxtomm;

color=zeros(length(xt),1);
c=1;
for k=1:10:length(traj)
    temp=traj(k).length;
    color(c:c+temp-1)=k;
    c=c+temp;
end

figure
[~,uniq] = unique(color);
for ii=1:10:numel(uniq)-1
    plot(xt(uniq(ii):uniq(ii+1)-1),yt(uniq(ii):uniq(ii+1)-1),'.');hold on
end

axis equal
box on; grid on
xlabel('mm')
ylabel('mm')

folderout = 'rawtrajs/';
mkdir(folderout)
savefig_fcb([folderout 'raw_trajs'],8,6,'pdf')
savefig_fcb([folderout 'raw_trajs'],8,6,'fig')

%% Filter

%load('traj_conc.mat')
traj=traj_conc;

%%% Find proper filter width
if pi==pi % set pi==pi if you want this to run. No need to do it every time.

    [s(1), m(1), w]=findFilterWidth_PTV(traj,'x');
    [s(2), m(2), w]=findFilterWidth_PTV(traj,'y');


    mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
    color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
    color1 = '#476d76';
    
    figure;
    yyaxis left
    loglog(w,s(1).vx,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
    loglog(w,s(2).vx,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
    hold off
    
    yyaxis right
    loglog(w,s(1).ax,'^-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on;
    loglog(w,s(2).ax,'^-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
    
    plot([10 10],ylim,'--',Color=color1,LineWidth=2)
    
    yyaxis left
    legend('$V_x$','$V_y$','$A_x$','$A_y$','interpreter','latex',Location='southwest');
    title('$std.(w)$','interpreter','latex',FontWeight='bold')
    xlabel('$filter\ width\ w$','interpreter','latex',FontWeight='bold')
    
    yyaxis left
    ylabel('$\sigma_{v}$','interpreter','latex',FontWeight='bold')
    yyaxis right
    ylabel('$\sigma_{a}$','interpreter','latex',FontWeight='bold')
    
    grid on
    axis padded
    
    folderout = [pathout filesep 'filter_check' filesep];
    mkdir(folderout)
    savefig_fcb([folderout 'filter_check'],8,6,'pdf')
    savefig_fcb([folderout 'filter_check' ],8,6,'fig')
    save([pathout filesep 'filter_check' filesep 'output_filtering.mat'],'s','m','w')
end

%%% Compute velocity and acceleration through a gaussian filter

w = 10; % filter intensity

[~, trajf]=compute_vel_acc_traj(traj_conc,fps,w);

save('trajf','trajf')


%% Plot filtered traj

%load('trajf.mat')

xt=vertcat(trajf.xf);
yt=vertcat(trajf.yf);
color=zeros(length(xt),1);
c=1;
for k=1:10:length(trajf)
    temp=trajf(k).lengthf;
    color(c:c+temp-1)=k;
    c=c+temp;
end


figure
[~,uniq] = unique(color);
for ii=1:10:numel(uniq)-1
    plot(xt(uniq(ii):uniq(ii+1)-1),yt(uniq(ii):uniq(ii+1)-1),'.');hold on
end
axis equal
box on; grid on
xlabel('mm')
ylabel('mm')


folderout = 'filteredtrajs/';
mkdir(folderout)
savefig_fcb([folderout 'filteredtrajs'],8,6,'pdf')
savefig_fcb([folderout 'filteredtrajs'],8,6,'fig')

%% Concatenate trajectories

if 1==pi

    trajs_conc = [];
    
    load('trajs_EXP11_a.mat')
    trajs_conc = [trajs_conc trajsf];
    clear trajsf
    1
    load('trajs_EXP11_b.mat')
    trajs_conc = [trajs_conc trajsf];
    clear trajsf
    2
    load('trajs_EXP11_c.mat')
    trajs_conc = [trajs_conc trajsf];
    clear trajsf

end

save('traj_conc','trajs_conc')
