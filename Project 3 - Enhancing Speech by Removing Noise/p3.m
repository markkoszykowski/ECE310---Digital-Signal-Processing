%% ECE310 - DSP Project 3
%% Tamar Bacalu, Mark Koszykowski, Sam Shersher
clc;
clear;
close all;

load projIB;
soundsc(noisy,fs);
pause(5);

% Find cutoff frequencies in Hz and rad
passband = 2500;
stopband = 4000;
wp = passband*2/fs; 
ws = stopband*2/fs;

% Set variables for the stop/pass band ripples
iirrp = 3;
firrp = 1.5;
iirG = 10^(40/20);
firG = 10^(38.5/20);
rs = 100;

%% Butterworth
close all;

% Calcuate order
[n, Wn] = buttord(wp,ws,iirrp,rs);
fprintf("Order of Butterworth filter is %d \n", n)

% Convert to SOS
[z, p, k] = butter(n,Wn);
[sos, g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

%In direct form implementation each of the poles and zeros are 1 multiply
fprintf("Num of multiplies: %d per sample\n",2*n + 1);

gain = dfilt.scalar(iirG);
% Apply amplification
butt = cascade(gain,hd);
% Apply filter
filtered = filter(butt,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(butt,'Butterworth')

%% Cheby1
close all;

% Calculate order
[n, Wn] = cheb1ord(wp,ws,iirrp,rs);
fprintf("Order of Cheby1 filter is %d \n", n)

% Convert to SOS
[z, p, k] = cheby1(n,iirrp,wp);
[sos, g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

fprintf("Num of multiplies: %d per sample\n",2*n + 1);

gain = dfilt.scalar(iirG);
% Apply Amplification
cheby1 = cascade(gain,hd);
% Apply filter
filtered = filter(cheby1,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(cheby1,'Cheby1')

%% Cheby2
close all;

% Calculate order
[n, Wn] = cheb2ord(wp,ws,iirrp,rs);
fprintf("Order of Cheby2 filter is %d \n", n)

% Convert to SOS
[z, p, k] = cheby2(n,rs,ws);
[sos, g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

fprintf("Num of multiplies: %d per sample\n",2*n+1);

gain = dfilt.scalar(iirG);
% Apply amplification
cheby2 = cascade(gain,hd);
% Apply filter
filtered = filter(cheby2,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(cheby2,'Cheby2')

%% Elliptic
close all;

% Calculate order
[n, Wn] = ellipord(wp,ws,iirrp,rs);
fprintf("Order of Elliptic filter is %d \n", n)

% Convert to SOS
[z, p, k] = ellip(n,iirrp,rs,wp);
[sos, g] = zp2sos(z,p,k);
hd = dfilt.df2sos(sos,g);

fprintf("Num of multiplies: %d per sample\n",2*n + 1);

gain = dfilt.scalar(iirG);
% Apply amplification
ellip = cascade(gain, hd);
% Apply filter
filtered = filter(ellip,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(ellip,'Elliptic');

%% Parks McLellan
close all;

% Calculate order
[n, fo, ao, w] = firpmord([passband stopband],[1 0], [(10^(firrp/20)-1)/(10^(firrp/20)+1)  10^(-rs/20)], fs);
fprintf("Order of Parks McLellan filter is %d \n", n)

% Convert to DF1
b = firpm(n,fo,ao,w);
hd = dfilt.df1(b);

%# of multiplies in FIR without exploiting symmetry is just n+1
fprintf("Num multiplies is %d per unit sample \n",n+1)

gain = dfilt.scalar(firG);
% Apply amplification
parks = cascade(gain, hd);
% Apply filter
filtered = filter(parks,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(parks,'Parks McLellan');

%% Kaiser
close all;

% Calculate order
[n, Wn, beta, ftype] = kaiserord([passband stopband], [1 0], [(10^(firrp/20)-1)/(10^(firrp/20)+1)  10^(-rs/20)], fs);
fprintf("Order of Kaiser filter is %d \n", n)

% Convert to DF1
b = fir1(n,Wn,ftype,kaiser(n+1,beta));
hd = dfilt.df1(b);

fprintf("Num multiplies is %d per unit sample \n",n+1)

gain = dfilt.scalar(firG);
% Apply amplification
kais = cascade(gain, hd);
% Apply filter
filtered = filter(kais,noisy);
soundsc(filtered,fs);
pause(5);

% Get graphs
filtinfo(kais,'Kaiser');

%% All the filters seem to do a reasonable job of blocking out the noise and making the track sound coherent

function filtinfo(filt,name)
    figure;
    subplot(3,1,1);
    [H,W] = freqz(filt);
    plot(W,20*log10(abs(H)));
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency (rad/sample)');
    ylabel('Magnitude (dB)');
    title(['Frequency Response of ',name]);
    
    subplot(3,1,2);
    plot(W,abs(H));
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency (rad/sample)');
    ylabel('Magnitude');
    title(['Frequency Response of ',name]);

    subplot(3,1,3);
    [g, w] = grpdelay(filt);
    plot(w,g);
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency (rad/sample)');
    ylabel('Group Delay (Samples)');
    title(['Group Delay of ',name]);

    zplane(filt);
    title(['Pole-Zero Plot of ',name]);

    figure;
    imp=zeros(1,100);
    imp(1)=1;
    impresponse = filter(filt,imp);
    samps = 0:99;
    stem(samps, impresponse);
    xlabel('Samples (n)');
    ylabel('Amplitude');
    title(['Impulse Response of ',name]);
end