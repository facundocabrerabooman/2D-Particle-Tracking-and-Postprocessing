%% analysis faraday waves
close all, clear all, clear

name = '200Hz_drop';
angle_tube = atand(0.018);% (atan of slope)
corners = [588 1264 1264 588; 34 34 710 710]; % upl-upr-lwr-lwl

pathin = ['/Users/FC/Research/Drop_tower/Faraday_waves/Small_camera_canon/' name filesep];
video = VideoReader([pathin filesep name '.mp4']);
numFrames = video.numFrames;

video_track = VideoWriter([pathin filesep name]); %create video
video_track.FrameRate = 2;
vide_track.Quality = 100;
open(video_track);
idx=0;
for ii = 1:2:numFrames %

    ii/numFrames

 idx = idx+1;
 frame = read(video,ii);
 frame = imrotate(frame,angle_tube);
 frame = rgb2gray(frame);
 framecut = frame(corners(2,2):corners(2,3),corners(1,1):corners(1,2));

 %figure;imshow(framecut)
 %figure;imshow(BW)
 
 framecut = imadjust(framecut);
 framecut_f = imgaussfilt(framecut,2);

 BW = edge(framecut,'canny');
 %figure;imshow(BW)

 
 %%%%%% track structures (failed)

%[cc, r] = imfindcircles(BW,[2 10]);
 %PP=regionprops(BW, 'Centroid', 'Eccentricity', 'EquivDiameter');
% PP([PP.Eccentricity]>.45 ,:) = [];
% %PP([PP.EquivDiameter]>8,:) = [];
% %PP([PP.EquivDiameter]<3,:) = [];
% r = [PP.EquivDiameter]./2;
% cc = cat(1,PP.Centroid);


% clf;imshow(framecut);hold on
% viscircles(cc, r,'EdgeColor','b');
%%%%%%


%%%%%%%%% figs for video

%f28=figure(28);
%  imshow(BW)

f28=figure('visible', 'off');
pcolor(flipud(framecut))
shading interp
colormap(jet)
colorbar

% figure
% [X,Y] = meshgrid(0:size(framecut,2)-1,0:size(framecut,1)-1);        
% meshc(X, Y, framecut)                            
% grid on
% colormap(jet)         
% view(0,90)

%%%%%%%%%

writeVideo(video_track,getframe(f28));

%pause(0.5)
end
 close(video_track)
