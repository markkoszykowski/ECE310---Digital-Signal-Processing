%% ECE310 - DSP Project 5
%% Tamar Bacalu, Mark Koszykowski, Sam Shersher
clc;
clear;
close all;

load pj2data;

%% A.1

y1 = y(1:32);
n = linspace(0, 62, 63);

cyy2 = xcorr(y1, y1, 'biased');
cyy3 = conv(y1, flip(y1));

figure;
subplot(2, 1, 1);
plot(n, cyy2);
title("Autocorrelation of y_1[n] using xcorr()");
xlabel("n");
ylabel("Amplitude");
xlim([min(n) max(n)]);

subplot(2, 1, 2);
plot(n, cyy3);
title("Autocorrelation of y_1[n] using conv()");
xlabel("n");
ylabel("Amplitde");
xlim([min(n) max(n)]);

% If xcorr is set to 'unbiased', it calculates the unbiased estimate of the
% cross-correlation

%% A.2

% Fourier transform of autocorrelation function is positive, real function
% because it is the power density spectrum

phiy1y1 = cyy3/32;

cyy2dft = fft(cyy2, 64);
k = linspace(0, 63, 64);

figure;
subplot(2, 1, 1);
plot(k, abs(cyy2dft));
title("Absolute Value of the DFT of c_{y1y1}[m]");
xlabel("k");
ylabel("|C_{y1y1}(e^{j\omega_k})|");
xlim([min(k) max(k)]);

subplot(2, 1, 2);
plot(k, angle(cyy2dft));
title("Phase of the DFT of c_{y1y1}[m]");
xlabel("k");
ylabel("\angle C_{y1y1}(e^{j\omega_k})");
xlim([min(k) max(k)]);

% The 64 point DFT of phi_y1y1 is the same as the 64 point DFT of c_y1y1

phiy1y1dft = fft(phiy1y1, 64);

figure;
subplot(2, 1, 1);
plot(k, abs(phiy1y1dft));
title("Absolute Value of the DFT of \Phi^d_{y1y1}[m]");
xlabel("k");
ylabel("|\Phi^d_{y1y1}(e^{j\omega_k})|");
xlim([min(k) max(k)]);

subplot(2, 1, 2);
plot(k, angle(phiy1y1dft));
title("Phase of the DFT of \Phi^d_{y1y1}[m]");
xlabel("k");
ylabel("\angle \Phi^d_{y1y1}(e^{j\omega_k})");
xlim([min(k) max(k)]);

%% A.3

Y1 = abs(fft(y1, 64));
Y1 = Y1.^2;

Y2 = abs(fft(y(1:64), 64));
Y2 = Y2.^2;

figure;
subplot(3, 1, 1);
plot(k, abs(phiy1y1dft));
title("Absolute Value of 64-point DFT of \phi^d_{y1y1}[m]");
xlabel("k");
ylabel("|\Phi^d_{y1y1}(e^{j\omega_k})|");
xlim([min(k) max(k)]);

subplot(3, 1, 2);
plot(k, Y1);
title("Magnitude Squared of 64-point DFT of y_1[n]");
xlabel("k");
ylabel("|Y_1(e^{j\omega_k})|^2");
xlim([min(k) max(k)]);

subplot(3, 1, 3);
plot(k, Y2);
title("Magnitude Squared of 64-point DFT of first 64 points of y[n]");
xlabel("k");
ylabel("|Y_2(e^{j\omega_k})|^2");
xlim([min(k) max(k)]);

% A.3a and A.3b are the same except for a scaling factor. 

%% B.1

I_64 = Y1/64;
freq_des = downsample(Hejw2, floor(512/64));

figure;
plot(k, I_64, k, freq_des);
legend('64-point periodogram','Desired freq. response');
title('64-Point Periodogram and Desired Freq. Response'); 
xlabel('k');
ylabel('Magnitude');
xlim([min(k) max(k)]);

err = sum(abs((I_64-freq_des)).^2)/64;
disp( "The estimation error with only first 32-points of data is: " + err );

%% B.2

I_1024 = abs(fft(y, 1024)).^2/1024;
k2=0:1023;

figure;
plot(k2, I_1024, k*1024/64, freq_des);
legend('1024-point periodogram','Desired freq. response');
title('1024-Point Periodogram and Desired Freq. Response'); 
xlabel('k');
ylabel('Magnitude');
xlim([min(k2) max(k2)]);

err = sum(abs((downsample(I_1024,floor(1024/64))-freq_des)).^2)/64;
disp( "The estimation error with all 512-points of data is: " + err );

%% B.3

windows = reshape(y, [32 16]);
pdgrams = abs(fft(windows,64)).^2/64;
pdgram_avg = sum(pdgrams,2)/16;

figure;
plot(k, pdgram_avg, k, freq_des);
title('Periodogram Average and Desired Frequency');
legend('Periodogram Average','Desired Frequency');
xlabel('k');
ylabel('Magnitude');
xlim([min(k) max(k)]);

err = sum(abs(pdgram_avg'-freq_des).^2)/64;
disp( "The estimation error with periodogram averaging is: " + err );

% First periodogram that is less than the desired frequency response

%% B.4

phiyy = xcorr(y, y, 'unbiased');
phiyy = phiyy(512-15:512+15);
phiyydft = abs(fft(phiyy, 64));

figure;
plot(k, phiyydft, k, freq_des);
title('Blackman-Turkey Method and Desired Frequency');
legend('Blackman-Turkey Method','Desired Frequency');
xlabel('k');
ylabel('Magnitude');
xlim([min(k) max(k)]);

err = sum(abs(phiyydft-freq_des).^2)/64;
disp( "The estimation error with the Blackman-Turkey method is: " + err );

%% B.5

% Table in Command Window

% Blackman-Turkey method produced the lowest error

% The Blackman-Turkey method does so well without averaging because instead
% of attempting to do averaging, it smooths the output by using a window

win = triang(31);
phiyy = phiyy.*win';
phiyydft = abs(fft(phiyy, 64));

figure;
plot(k, phiyydft, k, freq_des);
title('Blackman-Turkey Method w/ Triangular Window and Desired Frequency');
legend('Blackman-Turkey Method','Desired Frequency');
xlabel('k');
ylabel('Magnitude');
xlim([min(k) max(k)]);

err = sum(abs(phiyydft-freq_des).^2)/64;
disp( "The estimation error with the Blackman-Turkey method and triangular window is: " + err );

% The triangular window does a better job for smoothing the signal as 
% opposed to a rectangular window
