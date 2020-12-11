%% ECE310 - DSP Project 4
%% Tamar Bacalu, Mark Koszykowski
clc;
clear;
close all;

load s1;
load s5;
load vowels;

%% 1

mu = 4e9;
fs = 5e6;
T = 2e-4;

t = 0:1/fs:T;
x = cos(2 * pi * mu * t.^2);

figure;
spectrogram(x, triang(256), 255, 256, fs, 'yaxis');
title("Spectrogram of Linear FM Chirp, \mu = 4e9");

%% 2

f1 = mu * t;
phi = 2 * pi * mu * t.^2;
f2 = (1/(2*pi)) * (diff(phi)./diff(t));

figure;
plot(t, f1, t(2:end), f2);
title('Instantaneous Frequency');
xlabel('t (s)');
ylabel('f (Hz)');
xlim([0 max(t)]);
legend('f_1(t)', 'f_2(t)', 'Location', 'northwest');
% f2 corresponds to the slope of the ridge

%% 3

mu = 1e10;

x2 = cos(2 * pi * mu * t.^2);

figure;
spectrogram(x2, triang(256), 255, 256, fs, 'yaxis');
title("Spectrogram of Linear FM Chirp, \mu = 1e10");
% Increasing the chirp rate causes the slope to increase dramatically
% and the frequency exceeds the Nyquist frequency

%% 4

fs = 8e3;

figure;
subplot(2, 1, 1);
spectrogram(s1, triang(1024), 1023, 1024, fs, 'yaxis');
title("Spectrogram of S1.mat Fine Frequency, Poor Termporal Resolution");
subplot(2, 1, 2);
spectrogram(s5, triang(1024), 1023, 1024, fs, 'yaxis');
title("Spectrogram of S5.mat Fine Frequency, Poor Termporal Resolution");

%% 5

figure;
subplot(2, 1, 1);
spectrogram(s1, triang(128), 127, 128, fs, 'yaxis');
title("Spectrogram of S1.mat Fine Temporal, Poor Frequency Resolution");
subplot(2, 1, 2);
spectrogram(s5, triang(128), 127, 128, fs, 'yaxis');
title("Spectrogram of S5.mat Fine Temporal, Poor Frequency Resolution");

%% 6

v = padarray(vowels, 2^nextpow2(length(vowels)) - length(vowels), 0, 'post');
s = spectrogram(v, rectwin(256), 128, 1024, fs, 'yaxis');

output = invstft(s, 1024);

figure;
subplot(3, 1, 1);
plot(vowels);
title('Vowels.mat');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);

subplot(3, 1, 2);
plot(output);
title('Inverse STFT Vowels.mat');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);

subplot(3, 1, 3);
plot(output(1:length(vowels)) - transpose(vowels));
title('Difference Between Original & Inverse STFT Vowels.mat');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);
ylim([-2 2]);

%% 7

sDown = downsample(transpose(s), 2);
sInv = invstft(transpose(sDown), 1024);

figure;
subplot(2, 1, 1);
plot(vowels);
title('Vowels.mat');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);

subplot(2, 1, 2);
plot(sInv);
title('Downsampled Vowels.mat');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)/2]);

soundsc(vowels, fs);
pause(2);
soundsc(sInv, fs);

%% Function for 6 & 7

function [output] = invstft(stft, n)
    stft = [stft(1:end-1, :); flipud(stft)];
    [~, len] = size(stft);
    output = zeros(1, 128*(len + 1));
    inv = ifft(stft, n, 'symmetric');
    for i = 1:len
        ind = 128*(i - 1) + 1;
        output(ind:ind+255) = output(ind:ind+255) + real(transpose(inv(1:256, i)));
    end
    output(129:length(output) - 128) = output(129:length(output) - 128) ./ 2;
end