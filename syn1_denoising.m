clc
clear
%%  --------------------- loading data ------------------------------
addpath(genpath('./subroutines'));
load('data_noise_new7.mat');
data.org = data_noise;
load original_synthetic_data.mat;
pure_noise = data_noise - data_e;
data.org = pure_noise + data_e;
% load('syntNoisy3_z.mat');
% data.org = syntNoisy3_z;

%%  ------------------ papameters input ----------------------------
sample = 6000;  % Sampling rate of seismic record   
% ------------- 
%  1、If  the band range of the effective signal is known in advance 
%      the band range of the effective signal can be set
%  2、If  the band range of the effective signal is unknown in advance,
%   the band range of the effective signal can be set of 0~ 1/sample/2
% 
opt.f_s = 30;   % The cut-off frequency1       
opt.f_e = sample/2;  %  The cut-off frequency2 
opt.dt = 1/sample;  % sampling time (sec)
%   opt.nrs is the starting sample when calculating the ROV curve.  
opt.nrs = 20;
%   opt.bwconn is the connectivity form (4 or 8) of the pixel.
opt.bwconn = 8;

%%  --------------- Start the denoising program ---------------------------
tic
dn =  para_free_denoising5(data, opt);
toc

[noise_free_tf, ~] = wsst(data_e, sample);
% pure_noise = data_noise - data_e;
%%  ------------------ map the denoising result -----------------------------
t = 1:length(data_noise);
tt = t/sample;
t_max = 2400/sample;

%%
font1 = 15;
line1 = 1;
tf_max = max(max(abs(dn.org_tf)));
clim2 = 0.3;

figure
subplot 221
plot(tt, data_e, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 222
pcolor(tt, dn.org_f, abs(noise_free_tf));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;

subplot 223
plot(tt, data_noise, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 224
pcolor(tt, dn.org_f, (abs(dn.org_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;
%%
for i =1
figure
subplot 321
plot(tt, data_noise./max(abs(data_noise)), 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 322
pcolor(tt, dn.org_f, abs(dn.org_tf));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;

subplot 323
plot(tt, dn.gcv_dw./max(abs(dn.gcv_dw)), 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 324
pcolor(tt, dn.org_f, (abs(dn.gcv_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 325
% plot(tt, dn.final_dw./max(abs(dn.final_dw)), 'LineWidth', line1);
plot(tt, dn.final_dw, 'LineWidth', line1);
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 326
pcolor(tt, dn.org_f, (abs(dn.final_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;
end
%%   match  original waveform
for i = 1
figure
subplot 221
noisy1 = plot(tt, data_noise, 'color', [0.75 0.75 0.75],'LineWidth', line1 + 1);
hold on;
org1 = plot(tt, data_e./max(abs(data_e)), 'b','LineWidth', line1);
% hold on;
% plot(tt, data_noise./max(abs(data_noise)), 'color', [0.75 0.75 0.75],'LineWidth', line1);
hold on;
den1 = plot(tt, dn.final_dw, 'r','LineWidth', line1);
xlim([0 t_max]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out', 'xticklabel', []);
legend([org1, noisy1, den1],'Original data','Noisy data','Denoised data','location','northwest', 'FontWeight','bold');
grid minor;
grid on;

axes('Position',[0.18,0.60,0.3,0.3]); % Generate subfigure
x = t(0.16*sample: 0.21*sample);
y1 = data_e(x);
y2 = dn.final_dw(x);
plot(tt(x),y1,'b','linewidth', 2);
hold on;
plot(tt(x),y2,'r','linewidth', line1);
set(gca,'FontWeight','bold');
grid minor;
grid on;

subplot 222
plot(tt, data_noise - dn.final_dw,'LineWidth', line1+ 1)
xlim([0 t_max]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out', 'xticklabel', []);
legend('Extracted noise', 'FontWeight','bold');

grid minor;
grid on;

subplot 223
pcolor(tt, dn.org_f, (abs(dn.final_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 224
pcolor(tt, dn.org_f, abs(dn.org_tf) - abs(dn.final_tf));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:200:1000], 'TickDir', 'out');
clim([0 tf_max*0.2]);
grid minor;
grid on;
end
%%  
%  soft and hard thresholding
%  hard-thresholding
[org_wsst_Tx, f, hard_Tx, hard_dw] = dn_Thresh1(data_noise, 'hard', sample);
[~, ~, soft_Tx, soft_dw] = dn_Thresh1(data_noise, 'soft', sample);
%  bandpass filtering
[bp_dw] = Two_D_filter_bp(data_noise, 1/sample, 200, 450, '1');
bp_wsst = wsst(bp_dw, sample);

% figure
for i = 1
figure
subplot 621
plot(tt, data_e, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
% title('Original signal', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 622
pcolor(tt, dn.org_f, abs(noise_free_tf));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out', 'XTicklabel', []);
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;


subplot 623
plot(tt, data_noise, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
% title('Noisy signal', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 624
pcolor(tt, dn.org_f, (abs(dn.org_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out', 'XTicklabel', []);
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 625
plot(tt, bp_dw, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
% title('Bandpass filtering (250-600 Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 626
pcolor(tt, dn.org_f, abs(bp_wsst));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out', 'XTicklabel', []);
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;


subplot 627
plot(tt, hard_dw, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
% title('Hard-thresholding ', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 628
pcolor(tt, dn.org_f, abs(hard_Tx));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out', 'XTicklabel', []);
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 629
plot(tt, soft_dw, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
% title('Soft-thresholding', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot(6, 2, 10)
pcolor(tt, dn.org_f, abs(soft_Tx));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out', 'XTicklabel', []);
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;

subplot(6, 2, 11)
plot(tt, dn.final_dw, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
% title('The proposed', 'FontWeight','bold');
set(gca, 'fontsize', font1,'ytick', [-1:0.5:1], 'TickDir', 'out')
grid minor;
grid on;

subplot (6, 2, 12)
pcolor(tt, dn.org_f, (abs(dn.final_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([0 1000]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Freq (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [0:500:1000], 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;

%   zoomed details
zoom1 = 0.16;
zoom2 = 0.22;
figure
subplot 611
plot(tt, data_e, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 612
plot(tt, data_noise, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 613
plot(tt, bp_dw, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 614
plot(tt, hard_dw, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
% xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 615
plot(tt, soft_dw, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', [-1:0.5:1], 'TickDir', 'out', 'XTicklabel', [])
grid minor;
grid on;

subplot 616
plot(tt, dn.final_dw, 'LineWidth', line1)
xlim([zoom1 zoom2]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Amp', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;
end
%%  SNR and rms comparison
all_data = [data_noise, bp_dw, hard_dw, soft_dw, dn.final_dw];
arrivals = 0.166*sample;
snr_win = 260;
signa1_w = all_data(arrivals : arrivals + snr_win, :);
pure_noise_w = all_data(arrivals - snr_win + 1: arrivals, :);
snr = std(signa1_w)./std(pure_noise_w);
rms_comparison = rms(all_data - data_e);
for i = 1:5
    % [cc, lags] = xcorr(all_data(arrivals : arrivals + snr_win, i), data_e(arrivals : arrivals + snr_win), 'normalized');
    [cc, lags] = xcorr(all_data(:, i), data_e, 'normalized');
    % CC(i) = abs(cc(lags == 0));
    CC(i) = max(abs(cc));
end
%% wchoherence
% [WCOH , WCS, F,COI,WTX,WTY] = wcoherence(data_e, dn.final_dw, sample);
% 
% figure
% imagesc(abs(WTX))
% figure
% imagesc(abs(WTY))
% 
% wcoherence(data_e, dn.final_dw, 6000,'PhaseDisplayThreshold',0.7)
% xwt = cross(org_wsst_Tx, conj(dn.final_tf));
% xwt = org_wsst_Tx.*conj(dn.final_tf);
% xwt(data_e, dn.final_dw);
