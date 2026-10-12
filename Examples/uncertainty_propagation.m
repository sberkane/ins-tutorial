function result = uncertainty_propagation(outputDir)
%UNCERTAINTY_PROPAGATION Monte Carlo propagation and covariance contours.
%   X0 = I, Sigma0 = 0
%   omega_nom = 0, specific force = [a0;0;0] - g
%   gyro noise only on z: delta_omega_z ~ N(0, sig_w^2)
%   The specified accelerometer noise covariance is zero.
%
% Compare:
%   - Monte Carlo trajectories + endpoint cloud in (x,y)
%   - 3σ "principal great circles" of covariance ellipsoid in algebra,
%     mapped to the group and projected to (x,y)

if nargin < 1, outputDir = []; end
[outputDir, cleanup] = experiment_setup(outputDir); %#ok<ASGLU>

% --------- LaTeX-friendly settings ---------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FONT_AXES   = 25;
FONT_LABEL  = 25;
FONT_LEGEND = 21;
LINE_WIDTH  = 2;
bg_color    = [1.0, 0.9, 0.9];

%
% Parameters
%
dt     = 0.05;           % time step
T      = 15.0;            % total time
K      = round(T/dt);    % number of steps
Nmc    = 500;           % Monte Carlo trials
Ngc    = 400;            % points per great circle
sigmaN = 3;              % 3σ boundary for great circles

g  = [0;0;9.81];         % gravity (world)
a0 = 1.0;                % constant specific force along body x
a  = [a0;0;0]-g;

deg   = pi/180;
sig_w = 20*deg;           % std of gyro noise on z (rad/s)

% Lateral accelerometer noise standard deviation
sig_ay = 0;           % m/s^2 

rng(7,'twister');

%
% Initial mean state: X0 = I
%
R0 = eye(3);
v0 = zeros(3,1);
p0 = zeros(3,1);
Xbar = Gamma(R0,v0,p0);

% Initial covariance: Sigma0 = 0 (in algebra coords [phi; dv; dp])
Sigma = zeros(9);

% Process noise covariance (6x6) for [gyro; accel]
% Gyro noise acts on z; the specified accelerometer noise is zero.
Qk = diag([0, 0, sig_w^2,  0, sig_ay^2, 0]);
% Small diagonal regularization permits Cholesky sampling of singular Qk.
Lq = chol(Qk + 1e-18*eye(6), 'lower');
%
% 1) Propagate mean + covariance (nominal omega = 0)
%
omega_nom = [0;0;0];

for k = 1:K
    % Preserve the random-stream offset used for the tutorial samples.
    randn(6,1);

    % A_k and L_k
    Ak = Ak_RIEKF(g, dt);
    Lk = Lk_RIEKF(Xbar, g, dt);
    % mean propagation
    Xbar = se23_propagate(Xbar, omega_nom, a, g, dt);
    

    % covariance propagation
    Sigma = Ak*Sigma*Ak' + Lk*Qk*Lk';
end

% final mean position
[~, ~, pf] = se23_extract(Xbar);

%
% 2) Monte Carlo propagation (store trajectories + endpoints)
%
Pmc_traj = zeros(2, K+1, Nmc);
Pmc_end  = zeros(2, Nmc);

for i = 1:Nmc
    X = Gamma(R0,v0,p0);

    % initial position
    p = X(1:3,5);
    Pmc_traj(:,1,i) = p(1:2);

    for k = 1:K
        w = Lq * randn(6,1);

        domega = w(1:3);          % yaw noise plus numerical regularization
        da     = w(4:6);          % numerical regularization only

        omega = domega;           % nominal omega is zero
        a_i   = a + da;

        X = se23_propagate(X, omega, a_i, g, dt);

        p = X(1:3,5);
        Pmc_traj(:,k+1,i) = p(1:2);
    end

    Pmc_end(:,i) = Pmc_traj(:,end,i);
end
%
% 3) "Principal great circles" of 3σ covariance ellipsoid, mapped to xy
%
Sigma = (Sigma + Sigma')/2;
[U, D] = eig(Sigma);
lam = diag(D);
[lam, idx] = sort(lam, 'descend');
U = U(:,idx);

theta = linspace(0,2*pi,Ngc);
pairs = [1 2; 1 3; 2 3];
curves = cell(3,1);

for c = 1:3
    i = pairs(c,1); j = pairs(c,2);
    XY = zeros(2,Ngc);

    for kk = 1:Ngc
        xi = sigmaN*( sqrt(max(lam(i),0))*cos(theta(kk))*U(:,i) + ...
                      sqrt(max(lam(j),0))*sin(theta(kk))*U(:,j) );

        phi = xi(1:3);
        dv  = xi(4:6);
        dp  = xi(7:9);

        E = Gamma(so3_exp(phi)', ...
                  -so3_left_jacobian(-phi)*dv, ...
                  -so3_left_jacobian(-phi)*dp);

        Xg = E * Xbar;
        p  = Xg(1:3,5);
        XY(:,kk) = p(1:2);
    end

    curves{c} = XY;
end

%
% Plot (trajectories + endpoint cloud + covariance curves)
%
figure('Color', bg_color,'Units','pixels','Position',[100 100 850 550]); hold on; grid on;
ax = gca; ax.FontSize = FONT_AXES; ax.LineWidth = 1.2;
% Light gray trajectories (background)
traj_color = [0.65 0.65 0.65];
for i = 1:Nmc
    plot(Pmc_traj(1,:,i), Pmc_traj(2,:,i), ...
         'Color', traj_color, ...
         'LineWidth', 0.5,'HandleVisibility','off');
end
% Endpoint cloud
plot(Pmc_end(1,:), Pmc_end(2,:), '.', 'MarkerSize', 15, 'Color', 0.5*traj_color,'HandleVisibility','off');

% Empirical Cartesian covariance, centered at the nominal endpoint
plot_cov_ellipse(pf(1:2), cov(Pmc_end.'), 3, [1.0 0.5 0.0], 2.5);   % orange
% Covariance mapped great circles
plot(curves{1}(1,:), curves{1}(2,:), '-', 'LineWidth', LINE_WIDTH, 'Color', [0 0.7 1]);
plot(curves{2}(1,:), curves{2}(2,:), '-', 'LineWidth', LINE_WIDTH, 'Color', [0 0.7 1]);
plot(curves{3}(1,:), curves{3}(2,:), '-', 'LineWidth', LINE_WIDTH, 'Color', [0 0.7 1]);
% Nominal endpoint
plot(pf(1), pf(2), 'o', ...
    'Color',[0 0.7 0], ...
    'MarkerFaceColor',[0 0.7 0], ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'HandleVisibility','off');
xlabel('$p_x$ (m)','FontSize',FONT_LABEL);
ylabel('$p_y$ (m)','FontSize',FONT_LABEL);
xlim([-10 130])
ylim([-55 55])
pbaspect([1.7 1 1])
legend({'Euclidean covariance','Lie-group covariance'},'Location','northwest','FontSize',FONT_LEGEND);

exportgraphics(gcf,fullfile(outputDir,'uncertainty_propagation.pdf'),'ContentType','vector','Padding',12);
result = struct('mean_state',Xbar,'covariance',Sigma, ...
    'trajectories',Pmc_traj,'endpoints',Pmc_end,'dt',dt,'T',T);
save(fullfile(outputDir,'uncertainty_propagation.mat'),'result');

end

function plot_cov_ellipse(mu, C, k, colorRGB, lw)
% Plot k-sigma covariance ellipse for 2D Gaussian N(mu, C)
% mu: 2x1, C: 2x2 symmetric PSD

C = (C + C')/2;             % force symmetry
[V,D] = eig(C);

% sort eigenvalues descending
[d, idx] = sort(diag(D), 'descend');
V = V(:,idx);

% protect against tiny negatives from numerics
d = max(d, 0);

t = linspace(0, 2*pi, 300);
circle = [cos(t); sin(t)];

% ellipse transform: mu + k * V*sqrt(D)*circle
A = V * diag(sqrt(d));
E = mu + k * (A * circle);

plot(E(1,:), E(2,:), '--', 'Color', colorRGB, 'LineWidth', lw);
end
