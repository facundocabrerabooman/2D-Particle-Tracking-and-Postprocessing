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
        %%% pre-butterw filter
%         [b,a] = butter(4,1/160,'low'); %i let pass 2300fps/160~14hz
%         d(:,4) = filtfilt(b, a, d(:,4));
%         d(:,3) = filtfilt(b, a, d(:,3));
%         d(:,2) = filtfilt(b, a, d(:,2));

        
        w1 = round(w/3);
        w1_vel = round(w/3);
        w1_acc = round(w/3);
        
        kpos = posfiltcoef(w1,w);
        trajf(counter).xf = conv(d(:,1),kpos,'valid');
        trajf(counter).yf = conv(d(:,2),kpos,'valid');
        trajf(counter).zf = conv(d(:,3),kpos,'valid');
        
        kvel = velfiltcoef(w1_vel,w);
        trajf(counter).uf = conv(d(:,1),kvel,'valid').*fps;
        trajf(counter).vf = conv(d(:,2),kvel,'valid').*fps;
        trajf(counter).wf = conv(d(:,3),kvel,'valid').*fps;
        
        kacc = accfiltcoef(w1_acc,w);
        trajf(counter).af = conv(d(:,1),kacc,'valid').*fps^2;
        trajf(counter).bf = conv(d(:,2),kacc,'valid').*fps^2;
        trajf(counter).cf = conv(d(:,3),kacc,'valid').*fps^2;

        trajf(counter).lengthf = numel(trajf(counter).af);
        trajf(counter).t_sec = (1:1:numel(trajf(counter).af))'/fps;
        trajf(counter).t = (1:1:numel(trajf(counter).af))';
        
        trajf(counter).t_sec_abs = (traj(ii).frames(2:end-1))./fps;
        trajf(counter).frames = traj(ii).frames(2:end-1);
        disp('careful with times, not absolute all start at 1')
        %figure(31), plot((1:1:numel(UVWf(:,1)))/fps,UVWf), title('translational velocity [mm/s]'),hold on
        
        %figure(32), plot((1:1:numel(ABCf(:,1)))/fps,ABCf), title('translational accel [mm/s2]'),hold on
  end
  
  end