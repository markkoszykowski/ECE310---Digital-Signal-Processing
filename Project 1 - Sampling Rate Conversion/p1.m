clear all;
close all;
clc;
%% Tamar Bacalu, Mark Koszykowski, Sam Shersher
%% ECE310 - Project 1

% Get data from audio file
fs2 = 24000;
[data, fs1] = audioread('Wagner.wav');
%soundsc(data, fs1)

pause(10);

% Test the simple Sample Rate Conversion function
test1 = srconvertSingleStage([1 zeros(1, 3000)]);
verify(test1);
% Use the simple Sample Rate Conversion function on audio data
out1 = srconvertSingleStage(data);
soundsc(out1, fs2);

pause(10);

% Test the PolyPhase Sample Rate Conversion function
test2 = srconvertPolyPhase([1; zeros(3000, 1)]);
verify(test2);
% Use the PolyPhase Sample Rate Conversion function on audio data
out2 = srconvertPolyPhase(data);
soundsc(out2, fs2);

function [out] = srconvertSingleStage(in) 
    upRate = 320;
    downRate = 147;
    
    filt = designfilt('lowpassfir', 'PassbandFrequency', (1/320), 'PassbandRipple', .03, 'StopbandFrequency', (1.2/320), 'StopbandAttenuation', 80, 'DesignMethod', 'kaiserwin');
    
    tic
        up = upsample(in, upRate);
        postFilter = filter(filt, up);
        out = downRate * downsample(postFilter, downRate);
    toc
    
    nMults = length(filt.Coefficients)*upRate;
    nAdds = nMults - 1;
    fprintf("Number of Multiplies Performed Per Input Sample: %d\n", nMults)
    fprintf("Number of Adds Performed Per Input Sample: %d\n", nAdds)
end

function [out] = srconvertPolyPhase(in) 
    % 320 = 2 * 2 * 2 * 2 * 2 * 2 * 5
    upRate = 320;
    % 147 = 3 * 7 * 7
    downRate = 147;
    
    nMults = 0; %Count the num of multiplies
    nAdds = 0; %Count the num of adds
    
    filt2 = designfilt('lowpassfir', 'PassbandFrequency', (1/2), 'PassbandRipple', .03, 'StopbandFrequency', (1.2/2), 'StopbandAttenuation', 90, 'DesignMethod', 'equiripple');
    filt5 = designfilt('lowpassfir', 'PassbandFrequency', (1/5), 'PassbandRipple', .03, 'StopbandFrequency', (1.2/5), 'StopbandAttenuation', 90, 'DesignMethod', 'equiripple');
    
    E2 = poly1(filt2.Coefficients, 2);
    E5 = poly1(filt5.Coefficients, 5);
    
    len2 = length(filt2.Coefficients); %Find the lengths of the
    len5 = length(filt5.Coefficients);
    lenSig = length(in);

    
    tic
        out = in;
    
        for val = 1:6
            temp1 = [];
            for i = 1:2
                temp2 = conv(E2(i,:), out);
                temp2 = [zeros(i-1, 1); upsample(temp2, 2)];
                temp1 = [temp1; zeros(length(temp2)-length(temp1), 1)] + temp2;
            end
            out = temp1;
            
            nMults = nMults + len2;
            nAdds = nAdds + 2*(len2/2 - 1);
        end
    
        temp1 = [];
        for i = 1:5
            temp2 = conv(E5(i,:), out);
            temp2 = [zeros(i-1, 1); upsample(temp2, 5)];
            temp1 = [temp1; zeros(length(temp2)-length(temp1), 1)] + temp2;
        end
        nMults = nMults + len5;
        nAdds = nAdds + 5*(len5/5 - 1);
        out = temp1;
    
        out = 3 * downsample(out, 3);
        for i = 1:2
            out = 7 * downsample(out, 7);
        end
    toc
    
    fprintf("Number of Multiplies Performed Per Input Sample: %d\n", nMults)
    fprintf("Number of Adds Performed Per Input Sample: %d\n", nAdds)
end