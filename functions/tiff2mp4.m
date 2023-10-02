%% detect particles on faraday waves
clear all, close all, clear

name = 'Test 7';

pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Pinning edge/' name '/images'];
pathout = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Pinning edge/' name '/post processing'];
files = dir([pathin filesep '*.tiff']);
numFrames = numel(files);

cd(pathout)

video_track = VideoWriter([pathout filesep name],'MPEG-4'); %create video
video_track.FrameRate = 60;
vide_track.Quality = 100;

open(video_track);


for ii = 1:1:1963 
    ii
    
    frame = imread([pathin filesep files(ii).name]);

            f=figure('visible', 'off');
            clf;imshow(frame);
            writeVideo(video_track,getframe(gcf))

end 

close(video_track)

