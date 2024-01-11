%% compute velocity and acceleration from traj. structure
%it will add that information to the original structure, also filtered
%positions
function [traj,trajf]=compute_vel_acc_traj(traj,fps,w)

counter=0;
trajf = [];

for ii=1:numel(traj)
    clear d
    
    d(:,1) = traj(ii).x;
    d(:,2) = traj(ii).y;
    d(:,3) = traj(ii).z;
  
    if size(d,1)>w*3

        counter = counter+1;

        w1 = round(w/3);
        w1_vel = round(w/3);
        w1_acc = round(w/3);
        
        kpos = posfiltcoef(w1,w);
        trajf(counter).Xf = conv(d(:,1),kpos,'valid');
        trajf(counter).Yf = conv(d(:,2),kpos,'valid');
        trajf(counter).Zf = conv(d(:,3),kpos,'valid');
        
        kvel = velfiltcoef(w1_vel,w);
        trajf(counter).Vx = conv(d(:,1),kvel,'valid').*fps;
        trajf(counter).Vy = conv(d(:,2),kvel,'valid').*fps;
        trajf(counter).Vz = conv(d(:,3),kvel,'valid').*fps;
        
        kacc = accfiltcoef(w1_acc,w);
        trajf(counter).Ax = conv(d(:,1),kacc,'valid').*fps^2;
        trajf(counter).Ay = conv(d(:,2),kacc,'valid').*fps^2;
        trajf(counter).Az = conv(d(:,3),kacc,'valid').*fps^2;

        trajf(counter).lengthf = numel(trajf(counter).Ax);
        trajf(counter).t_sec = (1:1:numel(trajf(counter).Ax))'/fps;
        trajf(counter).t = (1:1:numel(trajf(counter).Ax))';
        
        trajf(counter).t_sec_abs = (traj(ii).frames(2:end-1))./fps;
        trajf(counter).frames = traj(ii).frames(2:end-1);
 end
  
  end