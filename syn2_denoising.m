clc
clear
%%  --------------------- loading data ------------------------------
% load('data_noise_new7.mat');
% data.org = data_noise;
% load original_synthetic_data.mat;
% pure_noise = data_noise - data_e;
% data.org = pure_noise + data_e;
addpath(genpath('./subroutines'));
load('ex_synth+Noise.mat');
data.org = syntNoisy3_z;
load('ex_synth.mat');
data_e = synt_80_z(1:length(syntNoisy3_z));
data_noise = syntNoisy3_z;
figure
plot(data_e)

figure
plot(data_noise)
%%  ------------------ papameters input ----------------------------
sample = 10;  % Sampling rate of seismic record
% ------------- 
%  1、If  the band range of the effective signal is known in advance
%      the band range of the effective signal can be set
%  2、If  the band range of the effective signal is unknown in advance,
%   the band range of the effective signal can be set of 0~ 1/sample/2
%
opt.f_s = 0;   % The cut-off frequency1
opt.f_e = sample/2;  %  The cut-off frequency2
opt.dt = 1/sample;  % sampling time (sec)
%   opt.nrs is the starting sample when calculating the ROV curve.
opt.nrs = 8;
%   opt.bwconn is the connectivity form (4 or 8) of the pixel.
opt.bwconn = 4;

%%  --------------- Start the denoising program ---------------------------
tic
dn =  para_free_denoising0(data, opt);
toc

[noise_free_tf, ~] = wsst(data_e, sample);
% pure_noise = data_noise - data_e;
%%  ------------------ map the denoising result -----------------------------
t = 1:length(data_noise);
tt = t/sample;
t_max = length(data_noise)/sample;

%%
font1 = 15;
line1 = 1;
tf_max = max(max(abs(dn.org_tf)));
clim2 = 0.5;

f1 = 0;
f2 = sample/2;
f_tick = [f1:1:f2];

for i =1
figure
subplot 221
plot(tt, data_e(t), 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

subplot 222
pcolor(tt, dn.org_f, abs(noise_free_tf(:, t)));
shading interp;
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
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
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;
end
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
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
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
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 325
plot(tt, dn.final_dw./max(abs(dn.final_dw)), 'LineWidth', line1)
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
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;
end
%%
for i = 1
   figure
subplot 321
plot(tt, data_e(t), 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

axes('Position',[0.35656819610308 0.853122362869192 0.16 0.11]); % 生成子图 最简单的方式
x = t(26*sample: 34*sample);
y1 = data_e(x);
% y2 = dn.final_dw(x);
plot(tt(x),y1,'linewidth', 1);
set(gca,'FontWeight','bold');
ylim([-0.4 0.4])
grid minor;
grid on;

subplot 322
pcolor(tt, dn.org_f, abs(noise_free_tf(:, t)));
shading interp;
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 max(max(abs(noise_free_tf)))*clim2]);
grid minor;
grid on;


subplot 323
plot(tt, data_noise, 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

axes('Position',[0.355939660590823 0.553544303797465 0.16 0.11]); % 生成子图 最简单的方式
x = t(26*sample: 34*sample);
y1 = data_noise(x);
% y2 = dn.final_dw(x);
plot(tt(x),y1,'linewidth', 1);
set(gca,'FontWeight','bold');
ylim([-0.4 0.4])
grid minor;
grid on;

subplot 324
pcolor(tt, dn.org_f, (abs(dn.org_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on;

subplot 325
plot(tt, dn.final_dw./max(abs(dn.final_dw)), 'LineWidth', line1)
xlim([0 t_max]);
ylim([-1 1]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Normalized amplitude', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'TickDir', 'out')
grid minor;
grid on;

axes('Position',[0.35656819610308 0.253966244725737 0.16 0.11]); % 生成子图 最简单的方式
x = t(26*sample: 34*sample);
% y1 = data_e(x);
y1 = dn.final_dw(x);
plot(tt(x), y1,'linewidth', 1);
set(gca,'FontWeight','bold');
ylim([-0.4 0.4])
grid minor;
grid on;

subplot 326
pcolor(tt, dn.org_f, (abs(dn.final_tf)));
% imagesc(tt, (dn.org_f), (abs(dn.org_tf)));
shading interp;
ylim([f1 f2]);
xlim([0 t_max]);
xlabel('Time (s)', 'FontWeight','bold');
ylabel('Frequency (Hz)', 'FontWeight','bold');
set(gca, 'fontsize', font1, 'ytick', f_tick, 'TickDir', 'out');
clim([0 tf_max*clim2]);
grid minor;
grid on; 
end
