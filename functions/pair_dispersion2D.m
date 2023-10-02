%% Pair dispersion processing 2D for faraday waves

load('trajectories_to_pair_disp.mat')

%%
delta=pairdisp_proc(trajf_formatted);
%%
%%% Pair dispersion statistics for separation
bin_type='log'; 
dmin=1;
dmax=342;
Nbin=10;
delta_d=pairdisp_stats_d(delta,bin_type,dmin,dmax,Nbin);

stop
save('OUTPUT_pairdispersion_d254-342_nbin20_logtype.mat')
%% Time evolution

load('OUTPUT_pairdispersion_d254-342_nbin20_logtype.mat')
fps=60;
pxtomm = 0.1428;

d0=vertcat(delta_d.d0);
figure(1);clf;hold on
for k=1:length(delta_d)
    y=delta_d(k).d2; %mean square separation
    %y=delta_d(k).Npts; %number of pairs
    x=(0:length(y)-1)'/fps; %times
    %plot(y,'x')
    %plot(x,y,'x')
    plot(y./x.^2,'x') %compensated
    %plot(x./d0(k)^(2/3),y./d0(k)^2,'x') %normalized
end
stop
set(gca,'XScale','log','YScale','log');grid on
d0_m=d0*pxtomm; %from px to mm
xlabel('time [frames]'), ylabel('S2 = deltaD/t^2')%, title('Distances in mts (legend)')

%mkdir('pairdisp')
%savefig(['pairdisp' filesep 'S2_vs_time'])
%% Time evolution compensated by D0^2/3
fps=2900;
d0=vertcat(delta_d.d0);
d0_m=d0*1e-3; %from mm to m
figure(2);%clf;
hold on
for k=1:length(delta_d)
    y=delta_d(k).d2; %mean square separation
    %y=delta_d(k).Npts; %number of pairs
    x=(0:length(y)-1)'/fps; %time
    %plot(y,'x')
    %plot(x,y,'x')
    %plot(y./x.^2,'x') %compensated
    plot(x./d0(k)^(2/3),y./d0(k)^2,'x') %normalized
end
set(gca,'XScale','log','YScale','log');grid on
%title('Initial distances (0.5mm - 12 cm)')%legend(num2str(d0_m),'Location','northeast'), 
xlabel('time/D0^{2/3}'), ylabel('DeltaD / (t^2 * D0^{2/3})')%, title('Distances in mts (legend)')
xlim([5e-6 0.1]), ylim([1e-7 100])

savefig(['pairdisp' filesep 'DeltaD_time_compensated_D23'])
%% Structure function
a=2;b=80;
S2=zeros(length(delta_d),1);
for k=1:18%length(S2)
    y=delta_d(k).d2;
    if ~isnan(y)
    x=(0:length(y)-1)'/fps;
    y_comp=y./x.^2;
    %S2(k)=median(y_comp(a:b))*1e-6; %from px to mts 
    S2(k)=median(y_comp(a:b))*1e-3; %from mm to mts 
    end
end

d0_m=d0*1e-3; %from mm to m
figure(3);clf,loglog(d0_m,S2,'x');grid on,xlabel('d0 [mts]'), ylabel('$S2 = \frac{44}{9}\epsilon^{2/3} D0^{2/3}$','Interpreter', 'latex', 'Fontsize', 24)
hold on, plot(d0_m,100.*d0_m.^(2/3),'r'), hold off

savefig(['pairdisp' filesep 'S2_vs_d0'])
%% Epsilon
epsilon=(((9/44)*S2).^1.5).*d0_m.^-1;
figure(4);clf,loglog(d0_m,epsilon,'x');grid on, xlabel('d0 [mts]'), ylabel('$\epsilon = (\frac{S2*9}{44})^{3/2}\frac{1}{D0}$','Interpreter', 'latex', 'Fontsize', 24)

savefig(['pairdisp' filesep 'epsilon'])