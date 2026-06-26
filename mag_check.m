%% mag_check.m  -- simplified, ONE figure per plot (no subplots)
%  Loads an injected/magnetometer pair, downconverts EACH file by its own
%  detected ~25 kHz tone (handles the +/- sign flip and small offsets),
%  and shows the few core plots. Requires Signal Processing Toolbox.
%
%  IMPORTANT: the cross-correlation / dead-time plot is only meaningful if
%  the two files were captured SIMULTANEOUSLY on a shared clock. Your current
%  files were recorded on different days, so that number is not a real dead
%  time -- the spectrum / inst-freq / BER plots are still valid.
clear; clc; close all;

% ----------------------- set your pair + params -----------------------
injFile = 'SDRuno_20260619_123725_25000HZ_INJ_SIN.wav';   % reference
magFile = 'SDRuno_20260623_150805_25000HZ_MAGSIN.wav';    % magnetometer
fc          = 25e3;     % nominal tone / centre frequency [Hz]
Rb          = 400;      % bit rate (FSK/MSK only)
pattern     = '0110100001101001';   % known 16-bit pattern
isModulated = false;    % false = plain sine, true = FSK/MSK
% ----------------------------------------------------------------------

[xi, fs ] = loadIQ(injFile);
[xm, fsm] = loadIQ(magFile);
if fs ~= fsm, xm = resample(xm, fs, fsm); end

%% Figure 1: raw spectra (see where the tone is, its sign, interference)
figure('Color','w','Name','Spectra');
[Pi,f] = pwelch(xi,[],[],[],fs,'centered');
[Pm,~] = pwelch(xm,[],[],[],fs,'centered');
plot(f/1e3,10*log10(Pi),'b'); hold on; plot(f/1e3,10*log10(Pm),'r');
grid on; xlabel('kHz'); ylabel('dB/Hz');
legend('injected','magnetometer','Location','best');
title('Raw spectra');

%% Per-file signed tone detection, then downconvert + tight low-pass
f0i = detectTone(xi, fs, fc, 2500);
f0m = detectTone(xm, fs, fc, 2500);
fprintf('Detected tone:  injected %+.1f Hz   magnetometer %+.1f Hz\n', f0i, f0m);
if sign(f0i) ~= sign(f0m)
    fprintf('** I/Q sign is opposite between the two files (spectral inversion).\n');
end
bw = max(3*Rb, 2500);
yi = baseband(xi, fs, f0i, bw);
ym = baseband(xm, fs, f0m, bw);
gi = ifreq(yi, fs);
gm = ifreq(ym, fs);

%% Figure 2: instantaneous-frequency overlay (modulation / timing view)
figure('Color','w','Name','Inst-freq overlay');
ns = min([numel(gi) numel(gm) round(0.06*fs)]);   % first ~60 ms
t  = (0:ns-1)/fs*1e3;
plot(t, gi(1:ns)/1e3,'b'); hold on; plot(t, gm(1:ns)/1e3,'r');
grid on; xlabel('ms'); ylabel('inst. freq [kHz]');
legend('injected','magnetometer','Location','best');
title('Instantaneous frequency (each downconverted to its own tone)');

%% Figure 3: cross-correlation (dead time) -- ONLY valid if simultaneous
L = min(numel(yi), numel(ym)); a = yi(1:L); b = ym(1:L);
[c1,lags] = xcorr(b, a);
[c2,~   ] = xcorr(conj(b), a);     % also try conjugate (handles IQ inversion)
if max(abs(c2)) > max(abs(c1)), c = c2; conjUsed = 1; else, c = c1; conjUsed = 0; end
win = abs(lags) <= round(5e-3*fs); cc = abs(c); cc(~win) = -inf;
[~,im] = max(cc); dt = lags(im)/fs;
figure('Color','w','Name','Cross-correlation');
plot(lags/fs*1e3, abs(c),'b'); hold on; plot(dt*1e3, abs(c(im)),'ro','MarkerFaceColor','r');
grid on; xlabel('lag [ms]'); ylabel('|cross-correlation|');
title(sprintf('Peak at %.3f ms (conj=%d) - valid ONLY for simultaneous capture', dt*1e3, conjUsed));

%% Figure 4: BER vs known pattern (modulated files only; no sync needed)
if isModulated
    pat = double(pattern(:) == '1');
    sps = fs/Rb;
    tau = bestTiming(gm, sps);
    bits = decideBits(gm, sps, tau*sps);
    [ber, off, inv] = berPattern(bits, pat);
    fprintf('Magnetometer BER vs pattern = %.3f  (offset %d, inv %d, %d bits)\n', ...
             ber, off, inv, numel(bits));
    base = circshift(pat, -off); if inv, base = 1 - base; end
    ref = repmat(base, ceil(numel(bits)/numel(pat))+1, 1); ref = ref(1:numel(bits));
    nb = min(64, numel(bits));
    figure('Color','w','Name','BER vs pattern');
    stem(0:nb-1, ref(1:nb), 'b', 'filled'); hold on;
    stem(0:nb-1, double(bits(1:nb))+0.05, 'r', 'filled');
    err = find(ref(1:nb) ~= bits(1:nb));
    plot(err-1, 1.1*ones(size(err)), 'kx', 'LineWidth', 1.5);
    ylim([-0.2 1.4]); grid on; xlabel('bit #');
    legend('sent (pattern)','magnetometer','errors','Location','best');
    title(sprintf('BER = %.3f   (0.5 means it is not locking)', ber));
end

%% ============================ helpers =================================
function [x, fs] = loadIQ(p)
    if exist(p,'file') ~= 2, error('File not found: %s', p); end
    [y, fs] = audioread(p);
    if size(y,2) >= 2, x = complex(y(:,1), y(:,2)); else, x = hilbert(y(:,1)); end
    x = x(:);
end

function f0 = detectTone(x, fs, fc, halfband)
    N = min(numel(x), 2^19);
    w = 0.5 - 0.5*cos(2*pi*(0:N-1).'/(N-1));
    X = fftshift(fft(x(1:N).*w));
    f = (-N/2:N/2-1).'*(fs/N);
    band = (abs(f) >= fc-halfband) & (abs(f) <= fc+halfband);
    Xb = abs(X); Xb(~band) = 0;
    [~,k] = max(Xb); f0 = f(k);
end

function y = baseband(x, fs, f0, bw)
    t = (0:numel(x)-1).'/fs;
    y = x .* exp(-1j*2*pi*f0*t);
    nord = min(400, 2*floor((numel(x)-1)/6)); nord = max(nord,8);
    b = fir1(nord, min(bw,0.45*fs)/(fs/2));
    y = filtfilt(b,1,real(y)) + 1j*filtfilt(b,1,imag(y));
end

function g = ifreq(x, fs)
    ph = unwrap(angle(x(:)));
    g  = [ph(2)-ph(1); diff(ph)] * fs/(2*pi);
end

function tau = bestTiming(g, sps)
    taus = 0:0.02:0.98; q = zeros(size(taus));
    for i = 1:numel(taus)
        s = sampleBits(g, sps, taus(i)*sps);
        if isempty(s), q(i)=0; else, q(i)=mean(abs(s-mean(s))); end
    end
    [~,j] = max(q); tau = taus(j);
end

function s = sampleBits(g, sps, d)
    Nb = max(0, floor((numel(g)-d)/sps)-1);
    idx = round(((0:Nb-1)+0.5)*sps + d) + 1;
    idx = idx(idx>=1 & idx<=numel(g)); s = g(idx);
end

function bits = decideBits(g, sps, d)
    s = sampleBits(g, sps, d);
    if isempty(s), bits = false(0,1); return; end
    ss = sort(s); lo = ss(max(1,round(0.1*numel(ss)))); hi = ss(round(0.9*numel(ss)));
    bits = s > 0.5*(lo+hi);
end

function [ber, off, inv] = berPattern(bits, pat)
    bits = logical(bits(:)); pat = logical(pat(:)); Lp = numel(pat); Lr = numel(bits);
    best = 1; off = 0; inv = false;
    if Lr < Lp, ber = NaN; return; end
    for p = 0:Lp-1
        ref = repmat(circshift(pat,-p), ceil(Lr/Lp)+1, 1); ref = ref(1:Lr);
        for fl = [false true]
            r = ref; if fl, r = ~r; end
            e = mean(r ~= bits);
            if e < best, best = e; off = p; inv = fl; end
        end
    end
    ber = best;
end
