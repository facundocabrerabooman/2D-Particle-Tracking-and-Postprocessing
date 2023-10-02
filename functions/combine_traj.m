%% load the files
clear all, close all, clear

name = 'Test 10';

%files
fol1 = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/post processing'];
fol2 = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/post processing 2'];
fol3 = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/post processing 3'];

%output
pathout = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu' ...
    '/ISS-CASIS/Experimental Design/Container design/' name '/total'];

%% get the trajectories

cd(fol1)
load('traj_concf.mat')
trajf1 = trajf;
ptm1 = pxtomm;
clear trajf pxtomm

cd(fol2)
load('traj_concf.mat')
trajf2 = trajf;
ptm2 = pxtomm;
clear trajf pxtomm

cd(fol3)
load('traj_concf.mat')
trajf3 = trajf;
ptm3 = pxtomm;
clear trajf pxtomm

cd(pathout)

%% Combine the trajectories
%convert structure to tables
trajf1t = struct2table(trajf1);
trajf2t = struct2table(trajf2);
trajf3t = struct2table(trajf3);

%combine tables vertically
merge_traj = vertcat(trajf1t,trajf2t,trajf3t);

%convert back to stucture
tot_traj = table2struct(merge_traj);

save('tot_traj')
clear trajf1 trajf2 trajf3 trajf1t trajf2t trajf3t
%% Plot trajectories
pxtomm = 1;

xt=vertcat(tot_traj.xf).*pxtomm;
yt=vertcat(tot_traj.yf).*pxtomm;
color=zeros(length(xt),1);
c=1;
for k=1:length(tot_traj)
    temp=tot_traj(k).lengthf;
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

savefig('traj_concf_tot')
%% Histogram of position PXTOMM NEEDS TO BE SAME FOR ALL??

figure()
subplot(1,2,1)
histogram(vertcat(tot_traj.xf).*ptm1,50,'FaceColor','r');hold on
title ('x histogram positions')
subplot(1,2,2)
histogram(vertcat(tot_traj.yf).*ptm1,50,'FaceColor','g')
title ('y histogram positions')

savefig('hist_totpos')

%% Histogram of velocity PXTOMM NEEDS TO BE SAME FOR ALL??

figure()
subplot(1,2,1)
histogram(vertcat(tot_traj.uf).*ptm1,50,'FaceColor','r');hold on
title ('x histogram velocity')
subplot(1,2,2)
histogram(vertcat(tot_traj.vf).*ptm1,50,'FaceColor','g')
title ('y histogram velocity')

disp('Mean velx')
mean(vertcat(tot_traj.uf))
subplot(1,2,1)
text(250,5000,['meanu = ' num2str(mean(vertcat(tot_traj.uf)))],'FontSize', 12)
disp('Mean vely')
mean(vertcat(tot_traj.vf))
subplot(1,2,2)
text(250,5000,['meanv = ' num2str(mean(vertcat(tot_traj.vf)))],'FontSize', 12)

savefig('hist_totvel')

%% Velocity converge check
tot=numel(tot_traj);
numTraj=0;

for aa = 1:4
    
    numTraj = [tot,round(tot/2),round(tot/4),round(tot/8)];
    newTraj=tot_traj(1:numTraj(aa));

    meanu(aa)=mean(vertcat(newTraj.uf));
    meanv(aa)=mean(vertcat(newTraj.vf));

end 

    
figure()
plot(numTraj,meanu,'blue')
hold on
plot(numTraj,meanv,'red')
legend('Mean u','Mean v')

%% compute lag stats

[Rvx,D2x,Nptsvx,Ntrackvx]=lagstats_tracks(tot_traj,'uf',2);
[Rvy,D2y,Nptsvy,Ntrackvy]=lagstats_tracks(tot_traj,'vf',2);

[Rax,~,Nptsax,Ntrackax]=lagstats_tracks(tot_traj,'af',2);
[Ray,~,Nptsay,Ntrackay]=lagstats_tracks(tot_traj,'bf',2);

%% correlations RUN LAG STATS IF YOU WANT THIS STUFF OR IF YOU CHANGED THINGS AROUND

fps = 112;
t=(0:length(Rvx)-1)'/fps;

figure;
plot(t,Rvx./Rvx(1),'r-s'); hold on;
plot(t,Rvy./Rvy(1),'b-s'); hold on;
title('Rv')
fname = 'Rv';
xlim([0 0.1])
% savefig_FC(fname,8,6,'fig')
% savefig_FC(fname,8,6,'pdf')
savefig('Rv_tot')


figure;
plot(t,Rax./Rax(1),'r-s'); hold on;
plot(t,Ray./Ray(1),'b-s'); hold on;
title('Ra')
fname = 'Ra';
xlim([0 0.1])
% savefig_FC(fname,8,6,'fig')
% savefig_FC(fname,8,6,'pdf')

savefig('Ra_tot')

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
savefig('d2L_tot')


