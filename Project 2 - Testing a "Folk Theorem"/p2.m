%% ECE310 - DSP Project 2
%% Tamar Bacalu, Mark Koszykowski, Sam Shersher

clc;
clear;
close all;

% Load and play original data for reference
data = load('projIA.mat');
soundsc(data.speech, data.fs);

pause(5);
%% a

samps = 100;

% Get the impulse and frequency response along with their domain vectors
[Htime,T] = impz(data.b,data.a,samps);
[Hfreq,W] = freqz(data.b,data.a,samps);

% Get the group delay of the filter
gd = grpdelay(data.b,data.a,samps);

% Plot impulse response
figure;
subplot(2,2,1);
stem(T, Htime);
xlabel('Samples (n)');
ylabel('Amplitude');
title('Impulse Response of All Pass Filter (N=1)');

% Plot magnitude response with respect to frequency
subplot(2,2,2);
plot(W, abs(Hfreq));
xlabel('Frequency (\omega)');
ylabel('Magnitude');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
ylim([0 2]);
title('Mag. of Freq Response of All Pass Filter (N=1)');

% Plot group delay
subplot(2,2,3);
plot(W,gd);
xlabel('Frequency (\omega)');
ylabel('Group Delay');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
title('Group Delay of All Pass Filter (N=1)');
%% b

% Plot the pole-zero plot
subplot(2,2,4);
zplane(data.b, data.a);
title('Pole-Zero Plot of All Pass Filter (N=1)');
%% c

% Run the sound file through the simple all-pass filter and listen to it
res = filter(data.b,data.a,data.speech);
soundsc(res, data.fs);

pause(5);

% No noticable audio distortion
%% d

% Convert the transfer function into second order sections
[z,p,k] = tf2zpk(data.b, data.a);
sos_1 = zp2sos(z,p,k);

% Calculate the filter coefficients to N=50
sos_50 = [];
for i = 1:50
    sos_50 = [sos_50;sos_1];
end
samps = 5000;





%DF1
Hdf1 = dfilt.df1(data.b,data.a);
Hdf1 = dfilt.cascade(repmat(Hdf1,1,50));

% Get the impulse and frequency response along with their domain vectors
[Htime,T] = impz(Hdf1,samps);
[Hfreq,W] = freqz(Hdf1,samps);

% Get the group delay of the filter
gd = grpdelay(Hdf1,samps);

% Plot impulse response
figure;
subplot(2,2,1);
stem(T,Htime);
xlabel('Samples (n)');
ylabel('Amplitude');
title('Impulse Response of DF1 All Pass Filter (N=50)');

% Plot magnitude response with respect to frequency
subplot(2,2,2);
plot(W,abs(Hfreq));
xlabel('Frequency (\omega)');
ylabel('Magnitude');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
ylim([0 2]);
title('Mag. of Freq Response of DF1 All Pass Filter (N=50)');

% Plot group delay
subplot(2,2,3);
plot(W,gd);
xlabel('Frequency (\omega)');
ylabel('Group Delay');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
title('Group Delay of DF1 All Pass Filter (N=50)');

% Plot the pole-zero plot
subplot(2,2,4);
[z,p] = Hdf1.zpk;
zplane(z,p);
title('Pole-Zero Plot of DF1 All Pass Filter (N=50)');





% DF1 SOS
Hdf1sos = dfilt.df1sos(sos_50);

% Get the impulse and frequency response along with their domain vectors
[Htime,T] = impz(Hdf1sos,samps);
[Hfreq,W] = freqz(Hdf1sos,samps);

% Get the group delay of the filter
gd = grpdelay(Hdf1sos,samps);

% Plot impulse response
figure;
subplot(2,2,1);
stem(T,Htime);
xlabel('Samples (n)');
ylabel('Amplitude');
title('Impulse Response of DF1 SOS All Pass Filter (N=50)');

% Plot magnitude response with respect to frequency
subplot(2,2,2);
plot(W,abs(Hfreq));
xlabel('Frequency (\omega)');
ylabel('Magnitude');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
ylim([0 2]);
title('Mag. of Freq Response of DF1 SOS All Pass Filter (N=50)');

% Plot group delay
subplot(2,2,3);
plot(W,gd);
xlabel('Frequency (\omega)');
ylabel('Group Delay');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
title('Group Delay of DF1 SOS All Pass Filter (N=50)');

% Plot the pole-zero plot
subplot(2,2,4);
[z,p] = Hdf1sos.zpk;
zplane(z,p);
title('Pole-Zero Plot of DF1 SOS All Pass Filter (N=50)');





% DF2 SOS
Hdf2sos = dfilt.df2sos(sos_50);

% Get the impulse and frequency response along with their domain vectors
[Htime,T] = impz(Hdf2sos,samps);
[Hfreq,W] = freqz(Hdf2sos,samps);

% Get the group delay of the filter
gd = grpdelay(Hdf2sos,samps);

% Plot impulse response
figure;
subplot(2,2,1);
stem(T,Htime);
xlabel('Samples (n)');
ylabel('Amplitude');
title('Impulse Response of DF2 SOS All Pass Filter (N=50)');

% Plot magnitude response with respect to frequency
subplot(2,2,2);
plot(W,abs(Hfreq));
xlabel('Frequency (\omega)');
ylabel('Magnitude');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
ylim([0 2]);
title('Mag. of Freq Response of DF2 SOS All Pass Filter (N=50)');

% Plot group delay
subplot(2,2,3);
plot(W,gd);
xlabel('Frequency (\omega)');
ylabel('Group Delay');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
title('Group Delay of DF2 SOS All Pass Filter (N=50)');

% Plot the pole-zero plot
subplot(2,2,4);
[z,p] = Hdf2sos.zpk;
zplane(z,p);
title('Pole-Zero Plot of DF2 SOS All Pass Filter (N=50)');





% DF2 Transposed SOS
Hdf2tsos = dfilt.df2tsos(sos_50);

% Get the impulse and frequency response along with their domain vectors
[Htime,T] = impz(Hdf2tsos,samps);
[Hfreq,W] = freqz(Hdf2tsos,samps);

% Get the group delay of the filter
gd = grpdelay(Hdf2tsos,samps);

% Plot impulse response
figure;
subplot(2,2,1);
stem(T,Htime);
xlabel('Samples (n)');
ylabel('Amplitude');
title('Impulse Response of DF2 Transposed SOS All Pass Filter (N=50)');

% Plot magnitude response with respect to frequency
subplot(2,2,2);
plot(W,abs(Hfreq));
xlabel('Frequency (\omega)');
ylabel('Magnitude');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
ylim([0 2]);
title('Mag. of Freq Response of DF2 Transposed SOS All Pass Filter (N=50)');

% Plot group delay
subplot(2,2,3);
plot(W,gd);
xlabel('Frequency (\omega)');
ylabel('Group Delay');
xlim([0 pi]);
xticks([0 pi/4 pi/2 3*pi/4 pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
title('Group Delay of DF2 Transposed SOS All Pass Filter (N=50)');

% Plot the pole-zero plot
subplot(2,2,4);
[z,p] = Hdf2tsos.zpk;
zplane(z,p);
title('Pole-Zero Plot of DF2 Transposed SOS All Pass Filter (N=50)');
%% e

% Run the sound file through the simple DF1 all-pass filter and listen to it
res = filter(Hdf1,data.speech);
soundsc(res, data.fs);

pause(5);

% Run the sound file through the simple DF1 SOS all-pass filter and listen to it
res = filter(Hdf1sos,data.speech);
soundsc(res, data.fs);

pause(5);

% Run the sound file through the simple DF2 SOS all-pass filter and listen to it
res = filter(Hdf2sos,data.speech);
soundsc(res, data.fs);

pause(5);

% Run the sound file through the simple DF2 Transposed SOS all-pass filter and listen to it
res = filter(Hdf2tsos,data.speech);
soundsc(res, data.fs);

% The audio is much more robotic and sounds kinda of like a slinky
% The reason the audio sounds so distorted is a result of the group delay
% which is significantly larger for N=50