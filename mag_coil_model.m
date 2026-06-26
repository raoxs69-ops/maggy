%% mag_coil_model.m
%  3-axis coil model.  Relates coil currents (Ix,Iy,Iz) to the magnetic
%  field (Bx,By,Bz) via a linear model:   B = A*I + b0
%     b0 = the field at zero current (your ambient / Earth field)
%     A  = 3x3 coupling matrix (how each current moves each field axis)
%  Gives both directions and the relationship plot. Needs only base MATLAB.
clear; clc; close all;

D  = readmatrix('mag_data.csv','NumHeaderLines',1);   % Ix Iy Iz Mag X Y Z
I  = D(:,1:3);  mag = D(:,4);  B = D(:,5:7);

%% ---- clean the data (your file has some shifted / typo rows) ----
good = abs(mag - vecnorm(B,2,2)) <= max(0.05*abs(mag), 800);   % |B| must equal sqrt(X^2+Y^2+Z^2)
for it = 1:5                                                    % then drop large-residual rows
    coef = [I(good,:) ones(sum(good),1)] \ B(good,:);
    res  = vecnorm(B - [I ones(size(I,1),1)]*coef, 2, 2);
    md   = median(res(good));
    thr  = max(md + 3*1.4826*median(abs(res(good)-md)), 1500);
    ng   = good & (res < thr);
    if isequal(ng,good), break; end
    good = ng;
end
fprintf('Using %d of %d rows after cleaning.\n', sum(good), numel(mag));

%% ---- fit  B = A*I + b0 ----
coef = [I(good,:) ones(sum(good),1)] \ B(good,:);
A  = coef(1:3,:).';     % rows = [Bx;By;Bz],  cols = [Ix Iy Iz]
b0 = coef(4,:).';       % field at I = 0

R2 = 1 - sum(sum((B(good,:)-[I(good,:) ones(sum(good),1)]*coef).^2)) ...
        / sum(sum((B(good,:)-mean(B(good,:))).^2));
fprintf('\nBaseline field at I=0:  b0 = [%.0f  %.0f  %.0f] nT   |b0| = %.0f nT\n', b0, norm(b0));
disp('Coupling matrix A (nT per current unit):'); disp(A)
fprintf('Fit R^2 = %.4f   cond(A) = %.2f\n', R2, cond(A));

%% ---- the two directions ----
field_from_currents = @(Ivec) b0 + A*Ivec(:);          % currents -> field
currents_from_field = @(Bvec) A \ (Bvec(:) - b0);      % field    -> currents

%% =====================================================================
%  ENTER YOUR TARGET FIELD HERE (nT).  Example below nulls the field.
Bx = 0;   By = 0;   Bz = 0;
%% =====================================================================
Ireq   = currents_from_field([Bx;By;Bz]);
Bcheck = field_from_currents(Ireq);
fprintf('\nTarget field  [%.0f %.0f %.0f] nT   (|B| = %.0f nT)\n', Bx,By,Bz, norm([Bx;By;Bz]));
fprintf('  -> set currents   Ix = %.1f    Iy = %.1f    Iz = %.1f\n', Ireq);
fprintf('  (model check: those currents give [%.0f %.0f %.0f] nT)\n', Bcheck);
if any(abs(Ireq) > 100)
    warning('A current exceeds +/-100 -- beyond your calibrated range (extrapolation).');
end

%% ---- relationship plot: field vs each current axis (others = 0) ----
figure('Color','w','Position',[80 80 1320 430]);
names = {'I_x','I_y','I_z'};
for ax = 1:3
    other = setdiff(1:3, ax);
    sel = good & all(I(:,other)==0, 2);
    x = I(sel,ax); Bs = B(sel,:); [x,si] = sort(x); Bs = Bs(si,:);
    xr = linspace(min(x),max(x),60).';  Im = zeros(60,3); Im(:,ax) = xr;
    Bm = [Im ones(60,1)]*coef;
    subplot(1,3,ax);
    plot(x,Bs(:,1),'o', x,Bs(:,2),'s', x,Bs(:,3),'^', x,vecnorm(Bs,2,2),'d');
    hold on; set(gca,'ColorOrderIndex',1);
    plot(xr,Bm(:,1),'-', xr,Bm(:,2),'-', xr,Bm(:,3),'-', xr,vecnorm(Bm,2,2),'-','LineWidth',1.3);
    grid on; xlabel([names{ax} ' (current units)']); ylabel('field (nT)');
    legend({'B_x','B_y','B_z','|B|'},'Location','best');
    title(['Field vs ' names{ax}]);
end
sgtitle('Magnetic field vs coil current');
