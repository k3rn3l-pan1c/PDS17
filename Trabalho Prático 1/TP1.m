%% Trabalho Prático 1
clear; clc; close all;

% load filters coefficients
load filter_coeffs

% Allow changes in coeefficients
change_coeffs = 0;

%% Input signal
% Sampling Frequency
Fs = 192 * 10^3;

% Load signal
x.time = audioread('reportagem.wav');

% Number of samples and sample vector
N = length(x.time);
n = 0:N-1;

% time vector
t = n / Fs;

% Time domain
figure(1)
plot(t, x.time)
xlabel('Time (s)')
ylabel('Amplitude')
title('Time Domain | x[n]')

% Frequency domain
f = (n - N/2) * Fs / N;
x.freq = fft(x.time);
figure(2)
stem(f/1e3, fftshift(abs(x.freq)), '.')
title('Frequency Domain | X(z)');
xlabel('Frequency (KHz)')
ylabel('Amplitude')

% Spectogram
figure(3)
spectrogram(x.time, 2048, 1024, 2048, Fs, 'yaxis')
colorbar
title('Spectrogram')

%% Coefficients for lowpass filter
% Response Type:  Lowpass
% Design Method:  IIR - Butterworth
% Fs:             192000
% Fpass:          10000
% Fstop:          24000
% Apass:          1 dB
% Astop:          80 dB
% Match exactly:  stopband
%
% Order: 11
if change_coeffs
    % load filter design tool
    fdatool
    
    % save coefficients
    iir_lowpass.Num = Num;
    iir_lowpass.Den = Den;
    save('filter_coeffs', 'iir_lowpass');
end;
%% Filter audible signal
% Sampling frequency for audible signal
Fs_low = Fs / 4;

% Filtering using low pass coeffs
y1.lowpass_time = filter(iir_lowpass.Num,iir_lowpass.Den, x.time);

% Time Domain
figure(4)
plot(t, y1.lowpass_time)
title('Input signal after Lowpass')
xlabel('Time (s)')
ylabel('Amplitude')

% Downsample
y1.time_filter_down = y1.lowpass_time(1:Fs/Fs_low:end);

% Play
sound(y1.time_filter_down, Fs_low)

%% Coefficients for bandpass filter
% Response Type:  bandpass
% Design Method:  IIR - Butterworth
% Fs:             192000
% Fstop1:         40000
% Fpass1:         56000
% Fpass2:         72000
% Fstop2:         96000
% Astop1:         60 dB
% Apass:          1 dB
% Astop2:         60 dB
% Match exactly:  stopband
%
% Order: 18
if change_coeffs
    % load filter design tool
    fdatool
    
    % save coefficients
    iir_bandpass.Num = Num;
    iir_bandpass.Den = Den;
    save('filter_coeffs', 'iir_bandpass', '-append');
end;

%% Filter high frequency signal
y2.time = filter(iir_bandpass.Num, iir_bandpass.Den, x.time);

% Time Domain
figure(5)
plot(t, y2.time)
title('High frequency signal after bandpass')
xlabel('Time (s)')
ylabel('Amplitude')


y2.freq = fft(y2.time);

% Frequency Domain
figure(6)
plot(f/1e3, abs(fftshift(y2.freq)))
title('High frequency signal after bandpass')
xlabel('Frequency (KHz)')
ylabel('Amplitude')


%% Demodulate high frequency signal
% Mininum frequency of the high frequency signal
f_low = 56e3;

% Multiplication with cosine
y2.demodulated = y2.time .* cos(2*pi*f_low*t');

% Time Domains
figure(7)
plot(t, y2.demodulated)
xlabel('time (KHz)')
title('Demodulated high frequency signal')
xlabel('Time (s)');
ylabel('Amplitude');

y2.demodulated_freq = fft(y2.demodulated);

% Frequency Domain
figure(8)
plot(f/1e3, abs(fftshift(y2.demodulated_freq)))
title('Demodulated High frequency signal ');
xlabel('Frequency (KHz)')
ylabel('Amplitude');

%% Coefficients for lowpass filter (to filter high frequency demodulated signal)
% Response Type:  Lowpass
% Design Method:  IIR - Butterworth
% Fs:             192000
% Fpass:          20000
% Fstop:          48000
% Apass:          1 dB
% Astop:          80 dB
% Match exactly:  stopband
%
% Order: 10
if change_coeffs
    % load filter design tool
    fdatool
    
    % save coefficients
    iir_lowpass_2.Num = Num;
    iir_lowpass_2.Den = Den;

    save('filter_coeffs', 'iir_lowpass_2', '-append');
end;

%% Filter demodulated high frequency signal
y2.base_band = filter(iir_lowpass_2.Num, iir_lowpass_2.Den, y2.demodulated);

% Time domain
figure(9)
plot(t, y2.base_band )
title('High frequency signal in base band')
xlabel('Time (s)')
ylabel('Amplitude')

y2.base_band_freq = fft(y2.base_band);

% Frequency domain
figure(10)
plot(f/1e3, abs(fftshift(y2.base_band_freq)))
title('High frequency signal in base band')
xlabel('Frequency (KHz)')
ylabel('Amplitude');
sound(y2.base_band, Fs_low);

%% Downsample to 48 KHz (Audible signal)
y2.base_band_downsampled = y2.base_band(1:Fs/Fs_low:end);

% Time domain
figure(11)
plot(f(1:Fs/Fs_low:end), y2.base_band_downsampled )
title('Downsampled High frequency signal in base band')
xlabel('Time (s)')
ylabel('Amplitude')


y2.base_band_downsampled_freq = fft(y2.base_band_downsampled);

% New frequency vector for lower frequency
N_low = N/(Fs/Fs_low);
n_low = 0: N_low;
f_low = (n_low - N_low/2) * Fs_low / N_low;

% Time domain (downsampled)
figure(12)
plot(f_low/1e3, abs(fftshift(y2.base_band_downsampled_freq)))
title('Downsampled High frequency signal in base band')
xlabel('Frequency (KHz)')
ylabel('Amplitude')

% Play high frequency sound in audible band
sound(y2.base_band_downsampled, Fs_low);

%% Listen to signals simultaneasly
y.time = y1.time_filter_down + y2.base_band_downsampled;
sound(y.time, Fs_low);
