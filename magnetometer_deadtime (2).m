%% ========================================================================
%  magnetometer_deadtime.m   (v2)
%  Estimate magnetometer latency ("dead time"), the FSK deviation, and BER
%  by comparing an INJECTED reference IQ recording to the MAGNETOMETER's IQ
%  output, for CW (sine), FSK, and MSK at a 25 kHz operating point.
%
%  What changed in v2:
%    - Bit rate = 400 bit/s.
%    - FSK deviation is UNKNOWN -> estimated from the recording.
%    - A KNOWN 16-bit FSK pattern is used as ground truth for BER.
%    - Dead time from complex-IQ cross-correlation (fine) + an instantaneous-
%      frequency cross-correlation as a cross-check, both windowed to avoid
%      the pattern-period ambiguity.
%
%  Requires: Signal Processing Toolbox (xcorr, fir1, filtfilt, pwelch,
%            spectrogram, resample, hann).
%
%  >>> Measurement assumption <<<
%  Reference and magnetometer captures must share ONE time base (simultaneous
%  capture / common trigger / common clock) for the dead-time number to be
%  the true injected->measured latency.
%
%  Quick start: leave SELFTEST = true and Run (known 60-sample delay + a
%  synthetic 700 Hz deviation to confirm everything), then set SELFTEST=false.
%% ========================================================================
clear; clc; close all;

%% -------------------------------- SETTINGS -----------------------------
SELFTEST = true;        % true: synthetic data (no files). false: load WAVs.

% --- WAV files (I in channel 1, Q in channel 2). Edit paths. ---
F.sine.ref = 'inj_sine.wav';   F.sine.rx = 'mag_sine.wav';
F.fsk.ref  = 'inj_fsk.wav';    F.fsk.rx  = 'mag_fsk.wav';
F.msk.ref  = 'inj_msk.wav';    F.msk.rx  = 'mag_msk.wav';

% --- System / modem parameters ---
P.fc      = 25e3;       % operating / tuned centre frequency [Hz]
P.Rb      = 400;        % bit rate [bit/s]
P.pattern = '0110100001101001';   % known 16-bit FSK test pattern (repeats)
P.fdev    = [];         % FSK deviation [Hz]. [] = auto-estimate from the data.

% --- IQ domain: 'auto' | 'passband' (contains ~25 kHz) | 'baseband' (DC) ---
P.domain = 'auto';

% --- Dead-time search & timing sweep ---
P.maxDeadtime_s = 5e-3; % search dead time within +/- this (< half the 40 ms
                        % pattern period, so the right xcorr peak is picked)
P.fineBits = 2.0;       % BER-vs-timing-offset sweep: +/- this many bit periods
P.fineStep = 0.02;      % step [bit periods]
% ------------------------------------------------------------------------

P.patBits = double(P.pattern(:) == '1');   % 0/1 column vector
if numel(P.patBits) ~= 16
    warning('Pattern has %d bits (expected 16).', numel(P.patBits));
end

if SELFTEST, S = makeSelfTest(P); else, S = []; end

% ---- estimate FSK deviation up front (from the cleanest reference) ----
if isempty(P.fdev)
    if SELFTEST || exist(F.fsk.ref,'file') == 2
        [rf,~,fsF] = getPair('fsk', F, S, SELFTEST);
        P.fdev = estimateFSKdev(rf, fsF, P);
        fprintf('Estimated FSK deviation = %.0f Hz  (mod. index h = %.2f)\n', ...
                 P.fdev, 2*P.fdev/P.Rb);
    else
        P.fdev = 2*P.Rb;
        fprintf('No FSK reference found; assuming deviation = %.0f Hz.\n', P.fdev);
    end
end

mods = {'sine','fsk','msk'};
R = struct();
for i = 1:numel(mods)
    m = mods{i};
    if ~SELFTEST && ~(exist(F.(m).ref,'file')==2 && exist(F.(m).rx,'file')==2)
        fprintf('Skipping %s (files not found).\n', m); continue;
    end
    [ref, rx, fs] = getPair(m, F, S, SELFTEST);
    R.(m) = analyze(m, ref, rx, fs, P);
end
printSummary(R, P, S, SELFTEST);

%% =========================== LOCAL FUNCTIONS ============================

function out = analyze(m, ref, rx, fs, P)
    L = min(numel(ref), numel(rx));
    ref = ref(1:L); rx = rx(1:L);
    out.fs = fs; out.type = m;

    dom = resolveDomain(P.domain, ref, fs, P.fc, P.Rb);
    out.domain = dom;
    bw = max(4*P.fdev, 6*P.Rb);
    if strcmp(dom,'passband')
        refb = lpf(mixdown(ref, fs, P.fc), fs, bw);
        rxb  = lpf(mixdown(rx,  fs, P.fc), fs, bw);
    else
        refb = lpf(ref, fs, bw);
        rxb  = lpf(rx,  fs, bw);
    end

    if strcmp(m,'sine')
        out = analyzeCW(out, ref, rx, refb, rxb, fs, P);
    else
        out = analyzeData(out, m, refb, rxb, fs, P);
    end
end

% ---------------------------- CW (sine) ---------------------------------
function out = analyzeCW(out, ref, rx, refb, rxb, fs, P)
    N = numel(ref);
    f = (0:N-1).'/N*fs; f(f>=fs/2) = f(f>=fs/2) - fs;
    Xr = fft(ref); Xx = fft(rx);
    [~,k] = max(abs(Xr)); ftone = abs(f(k));
    if ftone < P.Rb
        dphi = angle(mean(rxb)*conj(mean(refb))); fref = P.fc;
    else
        dphi = angle(Xx(k)*conj(Xr(k)));          fref = ftone;
    end
    tau = -dphi/(2*pi*fref); Tp = 1/fref;
    tau = mod(tau, Tp); if tau > Tp/2, tau = tau - Tp; end
    out.tone_Hz = fref; out.phaseDelay_s = tau; out.ambiguity_s = Tp;

    figure('Name','CW (sine) dead-time','Color','w');
    subplot(2,2,1);
        ns = min(N, round(5*fs/P.fc));
        plot((0:ns-1)/fs*1e6, real(ref(1:ns))); hold on;
        plot((0:ns-1)/fs*1e6, real(rx(1:ns)));
        xlabel('\mus'); ylabel('I'); legend('injected','magnetometer');
        title('Time domain (zoom)'); grid on;
    subplot(2,2,2);
        [c,lags] = xcorr(rx, ref); plot(lags/fs*1e6, abs(c)); hold on;
        xlim([-1 1]*5e6/P.fc);
        xlabel('lag [\mus]'); ylabel('|xcorr|');
        title('CW cross-correlation (ambiguous each period)'); grid on;
    subplot(2,2,3);
        [Pr,fp] = pwelch(ref,[],[],[],fs,'centered');
        [Px,~ ] = pwelch(rx ,[],[],[],fs,'centered');
        plot(fp/1e3,10*log10(Pr)); hold on; plot(fp/1e3,10*log10(Px));
        xlabel('kHz'); ylabel('dB/Hz'); legend('injected','magnetometer');
        title('Spectrum'); grid on;
    subplot(2,2,4);
        plot(real(ref(1:ns)), real(rx(1:ns)), '.'); axis equal; grid on;
        xlabel('I injected'); ylabel('I magnetometer');
        title(sprintf('Phase delay = %.3f \\mus  (\\pm%.1f \\mus ambig.)', ...
                       tau*1e6, Tp/2*1e6));
end

% --------------------------- FSK / MSK ----------------------------------
function out = analyzeData(out, m, refb, rxb, fs, P)
    sps = fs/P.Rb;
    fr  = ifreq(refb, fs);
    frx = ifreq(rxb,  fs);

    % ---- dead time: complex-IQ cross-correlation (windowed) ----
    [c,lags] = xcorr(rxb, refb);
    cc = abs(c); maxLag = round(P.maxDeadtime_s*fs);
    msk = abs(lags) <= maxLag; cc(~msk) = -inf;
    [~,im] = max(cc); lagIQ = parabolic(cc, im, lags);
    out.deadtime_s    = lagIQ/fs;
    out.deadtime_samp = lagIQ;
    out.deadtime_bits = lagIQ/sps;
    out.xcorr = abs(c); out.xcorrLags_s = lags/fs;

    % ---- cross-check: instantaneous-frequency cross-correlation ----
    [ci,~] = xcorr(frx - mean(frx), fr - mean(fr));
    ci = abs(ci); ci(~msk) = -inf;
    [~,im2] = max(ci); out.deadtime_if_s = parabolic(ci, im2, lags)/fs;

    % ---- demodulate measured signal & BER vs known pattern ----
    tauRx  = bestTiming(frx, sps);
    rxbits = decideBits(frx, sps, tauRx*sps);
    [berP, neP, nP, offP, invP] = berPattern(rxbits, P.patBits);
    out.BER = berP; out.bitErrors = neP; out.bitsCompared = nP;
    out.patOffset = offP; out.patInv = invP;

    % cross-check BER vs the injected reference's own demod
    txbits = decideBits(fr, sps, bestTiming(fr,sps)*sps);
    out.BER_ref = berRef(txbits, rxbits, 30);

    % measured deviation (for the histogram / report)
    sC = decideSamples(frx, sps, tauRx*sps); med = median(sC);
    hi = median(sC(sC>=med)); lo = median(sC(sC<med));
    out.dev_meas = 0.5*(hi-lo); out.devLines = [lo hi];

    % ---- BER vs sampling-timing offset (eye / timing margin) ----
    base = circshift(P.patBits(:), -offP); if invP, base = ~base; end
    offs = -P.fineBits:P.fineStep:P.fineBits; berf = zeros(size(offs));
    for i = 1:numel(offs)
        rb = decideBits(frx, sps, (tauRx+offs(i))*sps);
        berf(i) = berPatFixed(rb, base);
    end
    out.fineOff_s = offs/P.Rb; out.fineBER = berf;
    out.timingMargin_s = marginWidth(offs/P.Rb, berf, 1e-2);

    % ---- aligned baseband for plots ----
    rxb_al = delaySig(rxb, -round(lagIQ)); fr2 = ifreq(rxb_al, fs);

    % =========================== figure ================================
    figure('Name',sprintf('%s dead-time, deviation & BER',upper(m)),'Color','w');
    ns = min(numel(refb), round(8*sps));

    subplot(3,3,1);
        plot((0:ns-1)/fs*1e3, real(refb(1:ns))); hold on;
        plot((0:ns-1)/fs*1e3, real(rxb_al(1:ns)));
        xlabel('ms'); ylabel('I (baseband)'); legend('inj','mag (aligned)');
        title('Baseband I (zoom)'); grid on;
    subplot(3,3,2);
        plot((0:ns-1)/fs*1e3, fr(1:ns)/1e3); hold on;
        plot((0:ns-1)/fs*1e3, fr2(1:ns)/1e3);
        xlabel('ms'); ylabel('kHz'); legend('inj','mag');
        title('Instantaneous frequency'); grid on;
    subplot(3,3,3);
        plot(out.xcorrLags_s*1e6, out.xcorr); hold on;
        vline(out.deadtime_s*1e6,'r');
        xlim(out.deadtime_s*1e6 + [-1 1]*1e6*min(P.maxDeadtime_s,2e-3));
        xlabel('lag [\mus]'); ylabel('|xcorr|');
        title(sprintf('Dead time = %.1f \\mus', out.deadtime_s*1e6)); grid on;
    subplot(3,3,4);
        [Pr,fp] = pwelch(refb,[],[],[],fs,'centered');
        [Px,~ ] = pwelch(rxb ,[],[],[],fs,'centered');
        plot(fp/1e3,10*log10(Pr)); hold on; plot(fp/1e3,10*log10(Px));
        xlabel('kHz'); ylabel('dB/Hz'); legend('inj','mag');
        title('Baseband spectrum'); grid on;
    subplot(3,3,5);
        win = max(8, round(2*sps)); win = min(win, floor(numel(rxb)/4));
        spectrogram(rxb, hann(win), round(0.8*win), 1024, fs, 'centered','yaxis');
        title('Magnetometer spectrogram');
    subplot(3,3,6);
        eyeDiagram(frx/1e3, sps, 2);
        xlabel('symbol periods'); ylabel('freq [kHz]');
        title('Eye (mag inst-freq)'); grid on;
    subplot(3,3,7);
        histogram(sC/1e3, 50); hold on;
        vline(lo/1e3,'r'); vline(hi/1e3,'r');
        xlabel('kHz'); ylabel('count');
        title(sprintf('Tone freqs at bit centres | dev \\approx %.0f Hz', out.dev_meas));
    subplot(3,3,8);
        Lp = numel(P.patBits);
        ref0 = repmat(base, ceil(numel(rxbits)/Lp)+1, 1); ref0 = ref0(1:numel(rxbits));
        nb = min(64, numel(rxbits));
        stem(0:nb-1, ref0(1:nb), 'filled'); hold on;
        stem(0:nb-1, double(rxbits(1:nb))+0.04, 'filled');
        err = find(ref0(1:nb) ~= rxbits(1:nb));
        plot(err-1, 1.1*ones(size(err)), 'rx','LineWidth',1.5);
        ylim([-0.2 1.4]); xlabel('bit #'); legend('sent (pattern)','mag','errors');
        title(sprintf('Bits vs pattern  (BER = %.2g)', out.BER)); grid on;
    subplot(3,3,9);
        semilogy(out.fineOff_s*1e3, max(out.fineBER,1e-4)); hold on; vline(0,'k:');
        xlabel('timing offset [ms]'); ylabel('BER');
        title('BER vs sampling timing'); grid on;

    % standalone overlay so you can SEE how far the magnetometer lags
    timingOverlay(fr, frx, fs, sps, lagIQ, m, out.deadtime_s);
end

% ----------------------------- helpers ----------------------------------
function [ref, rx, fs] = getPair(m, F, S, SELFTEST)
    if SELFTEST, ref = S.(m).ref; rx = S.(m).rx; fs = S.fs; return; end
    [ref, fs1] = loadIQ(F.(m).ref);
    [rx,  fs2] = loadIQ(F.(m).rx);
    if fs1 ~= fs2, rx = resample(rx, fs1, fs2); end
    fs = fs1;
end

function [x, fs] = loadIQ(p)
    if exist(p,'file') ~= 2
        error('File not found: %s  (set SELFTEST = true to test without files).', p);
    end
    [y, fs] = audioread(p);
    if size(y,2) >= 2, x = complex(y(:,1), y(:,2));
    else,             x = hilbert(y(:,1)); end
    x = x(:);
end

function S = makeSelfTest(P)
    S.fs = 192e3; fs = S.fs;
    T = 0.4; N = round(T*fs); t = (0:N-1).'/fs;
    trueDelay = 60; snr = 20; rng(7);
    cw = exp(1j*2*pi*P.fc*t);
    S.sine.ref = cw;                              S.sine.rx = addAWGN(delaySig(cw,trueDelay),snr);
    devTest = 700;                                % synthetic FSK deviation to recover
    S.fsk.ref = cpfskPattern(P,fs,N,t,devTest,P.patBits);
    S.fsk.rx  = addAWGN(delaySig(S.fsk.ref ,trueDelay),snr);
    S.msk.ref = cpfskPattern(P,fs,N,t,P.Rb/4,P.patBits);
    S.msk.rx  = addAWGN(delaySig(S.msk.ref ,trueDelay),snr);
    S.trueDelay = trueDelay; S.devTest = devTest;
end

function x = cpfskPattern(P, fs, N, t, dev, pat)
    sps = fs/P.Rb; pat = pat(:); nb = ceil(N/sps);
    bits = pat(mod(0:nb-1, numel(pat)) + 1);
    sym  = 2*bits - 1;
    idx  = min(floor((0:N-1).'/sps)+1, nb);
    fb   = dev * sym(idx); ph = 2*pi*cumsum(fb)/fs;
    x    = exp(1j*(2*pi*P.fc*t + ph));
end

function y = delaySig(x, d)
    d = round(d);
    if d >= 0, y = [zeros(d,1); x(1:end-d)];
    else,      y = [x(1-d:end); zeros(-d,1)]; end
end

function y = addAWGN(x, snrdB)
    p = mean(abs(x).^2); n = p/(10^(snrdB/10));
    y = x + sqrt(n/2)*(randn(size(x)) + 1j*randn(size(x)));
end

function f = ifreq(x, fs)
    ph = unwrap(angle(x(:)));
    f  = [ph(2)-ph(1); diff(ph)] * fs/(2*pi);
end

function y = mixdown(x, fs, fc)
    t = (0:numel(x)-1).'/fs; y = x .* exp(-1j*2*pi*fc*t);
end

function y = lpf(x, fs, fcut)
    fcut = min(fcut, 0.45*fs);
    nord = min(256, 2*floor((numel(x)-1)/6)); nord = max(nord, 8);
    b = fir1(nord, fcut/(fs/2));
    y = filtfilt(b,1,real(x)) + 1j*filtfilt(b,1,imag(x));
end

function dom = resolveDomain(mode, ref, fs, fc, Rb)
    if ~strcmp(mode,'auto'), dom = mode; return; end
    N = numel(ref); f = (0:N-1).'/N*fs; f(f>=fs/2) = f(f>=fs/2) - fs;
    X = abs(fft(ref)).^2; bnd = max(2*Rb, 0.2*fc);
    eDC = sum(X(abs(f) < bnd));
    eFC = sum(X(abs(abs(f)-fc) < bnd));
    if eFC > eDC, dom = 'passband'; else, dom = 'baseband'; end
end

function dev = estimateFSKdev(ref, fs, P)
    dom = resolveDomain(P.domain, ref, fs, P.fc, P.Rb);
    if strcmp(dom,'passband'), x = mixdown(ref, fs, P.fc); else, x = ref; end
    x = lpf(x, fs, min(0.45*fs, 0.8*P.fc));     % image-safe wide filter
    f = ifreq(x, fs); sps = fs/P.Rb; tau = bestTiming(f, sps);
    s = decideSamples(f, sps, tau*sps); med = median(s);
    hi = median(s(s>=med)); lo = median(s(s<med));
    dev = 0.5*(hi - lo);
end

function lag = parabolic(c, im, lags)
    if im > 1 && im < numel(c)
        y0=c(im-1); y1=c(im); y2=c(im+1);
        lag = lags(im) + 0.5*(y0-y2)/(y0 - 2*y1 + y2 + eps);
    else
        lag = lags(im);
    end
end

function tau = bestTiming(f, sps)
    taus = 0:0.02:0.98; q = zeros(size(taus));
    for i = 1:numel(taus)
        s = decideSamples(f, sps, taus(i)*sps);
        if isempty(s), q(i)=0; else, q(i) = mean(abs(s - mean(s))); end
    end
    [~,j] = max(q); tau = taus(j);
end

function s = decideSamples(f, sps, delaySamp)
    Nb = floor((numel(f)-delaySamp)/sps) - 1; Nb = max(Nb,0);
    idx = round(((0:Nb-1)+0.5)*sps + delaySamp) + 1;
    idx = idx(idx>=1 & idx<=numel(f));
    s = f(idx);
end

function bits = decideBits(f, sps, delaySamp)
    s = decideSamples(f, sps, delaySamp);
    if isempty(s), bits = false(0,1); return; end
    thr = 0.5*(prctile(s,10) + prctile(s,90));
    bits = s > thr;
end

function [a,b] = overlapShift(tx, rx, s)
    tx = tx(:); rx = rx(:);
    if s >= 0, a = tx(1+s:end); b = rx(1:numel(a));
    else,      b = rx(1-s:end); a = tx(1:numel(b)); end
    L = min(numel(a), numel(b)); a = a(1:L); b = b(1:L);
end

function [ber, ne, n, off, inv] = berPattern(rx, pat)
    rx = logical(rx(:)); pat = logical(pat(:)); Lp = numel(pat); Lr = numel(rx);
    if Lr < Lp, ber = NaN; ne = 0; n = 0; off = 0; inv = false; return; end
    best = inf; ne = 0; n = Lr; off = 0; inv = false;
    for p = 0:Lp-1
        b = circshift(pat, -p); ref = repmat(b, ceil(Lr/Lp)+1, 1); ref = ref(1:Lr);
        for fl = [false true]
            r = ref; if fl, r = ~r; end
            e = sum(r ~= rx);
            if e < best, best = e; off = p; inv = fl; end
        end
    end
    ne = best; ber = best/Lr;
end

function r = berPatFixed(rx, base)
    Lr = numel(rx); if Lr < 10, r = 0.5; return; end
    ref = repmat(base(:), ceil(Lr/numel(base))+1, 1); ref = ref(1:Lr);
    r = sum(logical(ref) ~= logical(rx(:)))/Lr;
end

function ber = berRef(tx, rx, S)
    best = inf;
    for fl = [false true]
        rxf = logical(rx); if fl, rxf = ~rxf; end
        for s = -S:S
            [a,b] = overlapShift(tx, rxf, s);
            if numel(a) < 10, continue; end
            r = sum(logical(a)~=logical(b))/numel(a);
            if r < best, best = r; end
        end
    end
    ber = best;
end

function w = marginWidth(off, ber, thr)
    [~,j] = min(ber); lo = j; hi = j;
    while lo > 1 && ber(lo-1) <= thr, lo = lo - 1; end
    while hi < numel(ber) && ber(hi+1) <= thr, hi = hi + 1; end
    w = off(hi) - off(lo);
end

function eyeDiagram(f, sps, nsym)
    spsI = max(2, round(sps)); seg = round(nsym*spsI);
    Nseg = floor(numel(f)/seg); if Nseg < 1, return; end
    M = reshape(f(1:Nseg*seg), seg, Nseg);
    plot((0:seg-1)/spsI, M, 'Color',[0 0.3 0.8], 'LineWidth',0.1);
end

function vline(x, c)
    yl = ylim; plot([x x], yl, c);
end

function printSummary(R, P, S, SELFTEST)
    fprintf('\n==================== RESULTS SUMMARY ====================\n');
    if SELFTEST
        fprintf('(self-test) injected delay = %d samp = %.1f us | true dev = %d Hz\n', ...
                 S.trueDelay, S.trueDelay/S.fs*1e6, S.devTest);
    end
    fprintf('FSK deviation (used) = %.0f Hz  (h = %.2f)\n', P.fdev, 2*P.fdev/P.Rb);
    fprintf('Pattern = %s  (%d bits = %.1f ms period at %d bit/s)\n', ...
             P.pattern, numel(P.patBits), numel(P.patBits)/P.Rb*1e3, P.Rb);
    if isfield(R,'sine')
        s = R.sine;
        fprintf('CW  : phase delay = %+.3f us  (ambiguous every %.1f us = 1/fc)\n', ...
                 s.phaseDelay_s*1e6, s.ambiguity_s*1e6);
    end
    for mm = {'fsk','msk'}
        m = mm{1};
        if isfield(R,m)
            d = R.(m);
            fprintf(['%-4s: dead time = %.1f us (IQ) / %.1f us (inst-freq) | ' ...
                     'BER = %.3g (%d/%d, pattern) | BERref = %.3g | margin = %.2f ms | dev = %.0f Hz\n'], ...
                     upper(m), d.deadtime_s*1e6, d.deadtime_if_s*1e6, ...
                     d.BER, d.bitErrors, d.bitsCompared, d.BER_ref, ...
                     d.timingMargin_s*1e3, d.dev_meas);
        end
    end
    fprintf('--------------------------------------------------------\n');
    fprintf('Dead time (IQ xcorr) is the latency number. CW gives a finer\n');
    fprintf('sub-cycle value, ambiguous every %.1f us; they should agree mod that.\n', 1/P.fc*1e6);
    fprintf('========================================================\n');
end
