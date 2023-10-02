
% for 26 Hz
for oo = [1 4 7]

pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu/' ...
    'ISS-CASIS/Experimental Design/Container design/Test ' oo '/post processing'];
files = dir([pathin filesep '*.mat']);

um26(oo) = mean(vertcat(trajf.uf));
vm26(oo) = mean(vertcat(trajf.vf));

end 

% for 39 Hz
for oo = [2 5 8]

pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu/' ...
    'ISS-CASIS/Experimental Design/Container design/Test ' oo '/post processing/traj_concf.mat'];

um39(oo) = mean(vertcat(trajf.uf));
vm39(oo) = mean(vertcat(trajf.vf));

end 

% for 52 Hz
for oo = [3 6 9]

pathin = ['/Users/nataliefrank/Library/CloudStorage/GoogleDrive-natal8@pdx.edu/.shortcut-targets-by-id/15YNzxVZQM1XF1MixsqjektO5yZIwhJTu/' ...
    'ISS-CASIS/Experimental Design/Container design/Test ' oo '/post processing/traj_concf.mat'];

um52(oo) = mean(vertcat(trajf.uf));
vm52(oo) = mean(vertcat(trajf.vf));

end 
%% 
td = [3.1; 5.3; 8.1; 10.2; 14.8];
velvec = [31.59; 3.98; 1.71; 7.12; ];

figure()
plot(td, velvec, '--or','LineWidth',2,'MarkerSize',6)
hold on
% plot(td, um39, '--og','LineWidth',2,'MarkerSize',6)
% plot(td, um52, '--ob','LineWidth',2,'MarkerSize',6)
title ('Mean velocity across sizes')
%legend('26 Hz', '39 Hz', '52 Hz')

figure()
plot(td, vm26, '--om','LineWidth',2,'MarkerSize',6)
hold on
plot(td, vm39, '--ok','LineWidth',2,'MarkerSize',6)
plot(td, vm52, '--oc','LineWidth',2,'MarkerSize',6)
title ('V mean velocity across sizes')
legend('26 Hz', '39 Hz', '52 Hz')