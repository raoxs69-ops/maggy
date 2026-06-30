%% mag_coil_model.m
%  3-axis coil model for a Rb-87 / K magnetometer.
%     Field model:     B = A*I + b0
%     Larmor relation: f = 7 Hz/nT * |B|   <=>   |B| = f / 7
%  Works in BOTH directions and always reports magnitude AND frequency.
%  Needs only base MATLAB.
clear; clc; close all;

GAMMA = 7;                                  % Hz per nT (gyromagnetic ratio)
D  = readmatrix('mag_data.csv','NumHeaderLines',1);   % Ix Iy Iz Mag X Y Z
I  = D(:,1:3);  mag = D(:,4);  B = D(:,5:7);

%% ---- clean (your file has some shifted / typo rows) ----
good = abs(mag - vecnorm(B,2,2)) <= max(0.05*abs(mag), 800);
for it = 1:5
    coef = [I(good,:) ones(sum(good),1)] \ B(good,:);
    res  = vecnorm(B - [I ones(size(I,1),1)]*coef, 2, 2);
    md   = median(res(good)); thr = max(md + 3*1.4826*median(abs(res(good)-md)), 1500);
    ng   = good & (res < thr); if isequal(ng,good), break; end
    good = ng;
end
fprintf('Using %d of %d rows.\n', sum(good), numel(mag));

%% ---- fit  B = A*I + b0 ----
coef = [I(good,:) ones(sum(good),1)] \ B(good,:);
A  = coef(1:3,:).';   b0 = coef(4,:).';
fprintf('Baseline b0 = [%.0f %.0f %.0f] nT  |b0|=%.0f nT  (= %.0f Hz)\n', ...
        b0, norm(b0), GAMMA*norm(b0));
disp('Coupling A (nT per current unit):'); disp(A)

field_from_currents = @(Ivec) b0 + A*Ivec(:);          % currents -> field
currents_from_field = @(Bvec) A \ (Bvec(:) - b0);      % field    -> currents

%% =====================================================================
%  (1) FORWARD: enter CURRENTS  ->  field vector, magnitude, frequency
Ix = 72;   Iy = -32;   Iz = -10;        % <-- your currents
%% =====================================================================
Bv   = field_from_currents([Ix;Iy;Iz]);
Bmag = norm(Bv);
fprintf('\n--- From currents [%.1f %.1f %.1f] ---\n', Ix,Iy,Iz);
fprintf('  field      Bx = %.0f   By = %.0f   Bz = %.0f  nT\n', Bv);
fprintf('  MAGNITUDE  |B| = %.0f nT\n', Bmag);
fprintf('  frequency  f = 7*|B| = %.0f Hz\n', GAMMA*Bmag);

%% =====================================================================
%  (2) TUNE: enter desired FREQUENCY -> required magnitude -> currents
%      Orientation held at  By ~ 0,  Bz slightly below zero.
freq      = 25000;     % desired operating frequency [Hz]
Bz_target = -100;      % Bz slightly below zero [nT]   <-- set your value
BxSign    = -1;        % your data shows Bx is negative; use +1 to flip
%% =====================================================================
reqMag    = freq / GAMMA;                                 % required |B|
Bx_target = BxSign*sqrt(max(reqMag^2 - Bz_target^2, 0));
Btune     = [Bx_target; 0; Bz_target];
Itune     = currents_from_field(Btune);
fprintf('\n--- Tune to %.0f Hz ---\n', freq);
fprintf('  REQUIRED MAGNITUDE  |B| = %.0f nT\n', reqMag);
fprintf('  target field   Bx = %.0f   By = 0   Bz = %.0f  nT\n', Bx_target, Bz_target);
fprintf('  set currents   Ix = %.1f   Iy = %.1f   Iz = %.1f\n', Itune);
if any(abs(Itune) > 100)
    warning('A current exceeds +/-100 (beyond calibrated range).');
end

%% ---- relationship plot: field vs each current axis (others = 0) ----
figure('Color','w','Position',[80 80 1320 430]); names = {'I_x','I_y','I_z'};
for ax = 1:3
    other = setdiff(1:3, ax); sel = good & all(I(:,other)==0, 2);
    x = I(sel,ax); Bs = B(sel,:); [x,si] = sort(x); Bs = Bs(si,:);
    xr = linspace(min(x),max(x),60).'; Im = zeros(60,3); Im(:,ax) = xr;
    Bm = [Im ones(60,1)]*coef;
    subplot(1,3,ax);
    plot(x,Bs(:,1),'o', x,Bs(:,2),'s', x,Bs(:,3),'^', x,vecnorm(Bs,2,2),'d');
    hold on; set(gca,'ColorOrderIndex',1);
    plot(xr,Bm(:,1),'-', xr,Bm(:,2),'-', xr,Bm(:,3),'-', xr,vecnorm(Bm,2,2),'-','LineWidth',1.3);
    grid on; xlabel([names{ax} ' (current units)']); ylabel('field (nT)');
    legend({'B_x','B_y','B_z','|B|'},'Location','best'); title(['Field vs ' names{ax}]);
end
sgtitle('Magnetic field vs coil current');
