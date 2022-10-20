% Input Settings
tau = 120;
fs = 960;
filtOrder = 8;
fc = 16.5 * (1:6);
freqNotch = fc;
bw = 1;

% Parameters
bwn = 2*bw/fs;
nFreq = length(fc);
fn = 2*fc/fs;
Q = fn/bwn;

% Generate Signal
t = (0:1/fs:tau-1/fs)';
nSamples = length(t);
x = sin(2*pi*fc.*t);
n = 0.1*randn(nSamples,1);
xn = sum(x + n,2);

% Filter
h = notchFilterDesign(fn,bwn,filtOrder);
xnf = filter(h,xn);

% Process Spectrum and Spectrogram
[Xn,f] = periodogram(xn,hanning(length(xn)),2^17,fs,'psd','onesided');
[Xnf,~] = periodogram(xnf,hanning(length(xnf)),2^17,fs,'psd','onesided');

% Plot (PSD)
figure, plot(f,10*log10(Xn))
hold on, plot(f,10*log10(Xnf),'g')
axis([f(1) f(end) -100 20])
xlabel('Frequency [Hz]')
ylabel('PSD [dB re 1V^2/Hz]')
title('Power Spectral Density (Notch Filtering)')
legend('Original','Filtered')
fname = 'PSD (Comparison)';
savefig(fname)
print(fname,'-dpng','-r250')

% Plot (PSD, Detail)
figure, plot(f,10*log10(Xn))
hold on, plot(f,10*log10(Xnf),'g')
axis([fc(1)-2 fc(1)+2 -100 20])
xlabel('Frequency [Hz]')
ylabel('PSD [dB re 1V^2/Hz]')
title('Power Spectral Density (Notch Filtering, Detail)')
legend('Original','Filtered')
fname = 'PSD (Comparison, Detail)';
savefig(fname)
print(fname,'-dpng','-r250')

% Plot (Spectrogram)
figure, spectrogram(xn,hanning(1024),512,1024,960);
caxis([-90 -10])
title('Spectrogram \rm(Original)')
fname = 'Spectrogram (Original)';
savefig(fname)
print(fname,'-dpng','-r250')

figure, spectrogram(xnf,hanning(1024),512,1024,960);
caxis([-90 -10])
title('Spectrogram \rm(Notch-Filtered)')
fname = 'Spectrogram (Notch-Filtered)';
savefig(fname)
print(fname,'-dpng','-r250')
