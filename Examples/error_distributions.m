function result = error_distributions(outputDir)
%ERROR_DISTRIBUTIONS Direct-product and right-invariant pose uncertainty.
% Left : SO(3) x R^6   (direct product)  -> ellipsoids (no banana)
% Right: SE_2(3)       (right-invariant) -> bananas that grow with ||p_hat||

if nargin < 1, outputDir = []; end
[outputDir, cleanup] = experiment_setup(outputDir); %#ok<ASGLU>

% --------- GLOBAL FIGURE FONT SETTINGS (LaTeX-friendly) ---------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FONT_AXES   = 20;
FONT_LABEL  = 20;
LINE_WIDTH  = 2;
bg_color    = [1.0, 0.9, 0.9];   % red!10 style
COLOR_ELLIPSOID = [0.0 0.6 0.0];   % green
COLOR_BANANA    = [0.0 0.2 0.8];   % blue

%
% Parameters
%
p_list   = [-30 -20 -10 0 10 20 30];      % distances along x
N        = 20000;             % MC samples
Ns_total = 1600;              % total ticks per panel (distributed across p_list)
tick_len = 1.0;               % tick length (m)

%
% Monte Carlo (same Xi reused for all p_hat)
%
rng(7,'twister');

Rhat = eye(3);
vhat = [0;0;0];

deg = pi/180;

sig_yaw = 10*deg;     % yaw
sig_dv  = 0.3;        % velocity (not visualized)
sig_px  = 0.8;        % x position
sig_py  = 1.6;        % y position
sig_pz  = 1e-5;       % suppress z

Sigma = diag([ ...
    (1e-6)^2, (1e-6)^2, sig_yaw^2, ...   % roll, pitch ~ 0, yaw = 10°
    sig_dv^2 *ones(1,3), ...
    sig_px^2, sig_py^2, sig_pz^2 ...
]);

L  = chol(Sigma,'lower');
Xi = L * randn(9,N);

% Precompute yaw for the "direct product" attitude part (depends only on phi)
yaw_base = zeros(1,N);
for i = 1:N
    phi = Xi(1:3,i);
    Rerr = so3_exp(phi)';  % exp(-phi^) = exp(phi^)'
    yaw_base(i) = atan2(Rerr(2,1), Rerr(1,1));
end

% Subsample indices (same across all p_hat for fairness)
Ns_each = max(50, floor(Ns_total/numel(p_list)));
idx_all = round(linspace(1,N,Ns_each));

% Unit directions for ticks (based on yaw)
u_base = cos(yaw_base(idx_all));
v_base = sin(yaw_base(idx_all));

%
% Plot (LaTeX-ready)
%
figure('Color', bg_color,'Units','pixels','Position',[100 100 1200 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ===== Left: direct product =====
nexttile; hold on; grid on; axis equal;
ax = gca; ax.FontSize = FONT_AXES; ax.LineWidth = 1.2;
xlabel('$p_x$ (m)','FontSize',FONT_LABEL);
ylabel('$p_y$ (m)','FontSize',FONT_LABEL);
title({'Distribution on $SO(3)\times R^6$'}, ...
      'FontSize', FONT_LABEL);

% ===== Right: SE_2(3) =====
nexttile; hold on; grid on; axis equal;
ax = gca; ax.FontSize = FONT_AXES; ax.LineWidth = 1.2;
xlabel('$p_x$ (m)','FontSize',FONT_LABEL);
ylabel('$p_y$ (m)','FontSize',FONT_LABEL);
title({'Distribution on $SE_2(3)$'}, ...
      'FontSize', FONT_LABEL);

%
% Draw multiple clouds (ticks) for each p_hat
%
for j = 1:numel(p_list)
    phat = [p_list(j); 0; 0];
    Xhat = Gamma(Rhat, vhat, phat);

    % --- Left panel: p = phat - dp, yaw = yaw_base ---
    p_prod = zeros(3, Ns_each);
    for k = 1:Ns_each
        i  = idx_all(k);
        dp = Xi(7:9,i);
        p_prod(:,k) = phat - dp;
    end

    % --- Right panel: SE_2(3) mapping ---
    p_se23 = zeros(3, Ns_each);
    yaw_se = zeros(1, Ns_each);

    for k = 1:Ns_each
        i  = idx_all(k);
        xi = Xi(:,i);

        phi = xi(1:3);
        dv  = xi(4:6);
        dp  = xi(7:9);

        E = Gamma(so3_exp(phi)', ...
                  -so3_left_jacobian(-phi)*dv, ...
                  -so3_left_jacobian(-phi)*dp);

        X = E * Xhat;

        p_se23(:,k) = X(1:3,5);
        Rse = X(1:3,1:3);
        yaw_se(k) = atan2(Rse(2,1), Rse(1,1));
    end

    u_se = cos(yaw_se);
    v_se = sin(yaw_se);

    % --- Plot on left tile ---
    nexttile(1);
    q = quiver(p_prod(1,:), p_prod(2,:), ...
           tick_len*u_base, tick_len*v_base, 0, ...
           'LineWidth', 1.0, ...
           'Color', COLOR_ELLIPSOID);
    q.MaxHeadSize = 0;

    plot(phat(1), phat(2), 'kx', 'MarkerSize', 10, 'LineWidth', LINE_WIDTH);

    % --- Plot on right tile ---
    nexttile(2);
    q = quiver(p_se23(1,:), p_se23(2,:), ...
        tick_len*u_se, tick_len*v_se, 0, ...
        'LineWidth', 1.0, ...
        'Color', COLOR_BANANA);

    q.MaxHeadSize = 0;

    plot(phat(1), phat(2), 'kx', 'MarkerSize', 10, 'LineWidth', LINE_WIDTH);
end

%
% Link axes + padding
%
axs = findall(gcf,'Type','axes');
linkaxes(axs,'xy');

pad_frac = 0.08;
xl = axs(1).XLim;
dx = xl(2) - xl(1);
set(axs,'XLim', xl + [-pad_frac, pad_frac]*dx, ...
        'YLim', [-20 20]);

% Export
exportgraphics(gcf,fullfile(outputDir,'distribution_growth.pdf'), 'ContentType','vector','Padding',12);

result = struct('means',p_list,'samples',Xi,'display_indices',idx_all);
save(fullfile(outputDir,'error_distributions.mat'),'result');

end
