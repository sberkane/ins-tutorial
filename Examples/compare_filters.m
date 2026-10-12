function results = compare_filters(outputDir)
%COMPARE_FILTERS Reproduce the GPS-aided inertial navigation comparison.
%   RESULTS = COMPARE_FILTERS() writes figures and numerical results to the
%   repository results directory. COMPARE_FILTERS(OUTPUTDIR) uses OUTPUTDIR.
%
%   The initial attitude error is a 30-degree pitch rotation. Each sensing
%   configuration uses 50 paired noise realizations, with seeds 101:150.
%   This scenario was selected in exploratory runs using seeds 1:5.
%   Initial physical covariances are matched locally to first order. The
%   MEKF uses classical mean propagation; the InvEKF uses exact propagation.
%   Neither filter applies covariance resets.
%
%   The returned structure contains parameters, ground truth, per-run
%   metrics, error norms, seeds, and the first realization in each scenario.

thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(rootDir,'results');
end
outdir = outputDir;
if ~exist(outdir,'dir'), mkdir(outdir); end

session.rng = rng;
session.path = path;
session.textInterpreter = get(groot,'defaultTextInterpreter');
session.tickInterpreter = get(groot,'defaultAxesTickLabelInterpreter');
session.legendInterpreter = get(groot,'defaultLegendInterpreter');
cleanup = onCleanup(@() restore_session(session));
addpath(fullfile(rootDir,'Tools'));
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% ------------------------ Shared settings ------------------------
params.T   = 60;
params.dt  = 0.01;
params.t   = (0:params.dt:params.T)';
params.N   = numel(params.t);

params.gps_dt   = 0.1;
params.gps_step = max(1, round(params.gps_dt/params.dt));
params.gps_idx  = 1:params.gps_step:params.N;

params.g   = [0;0;9.81];

params.r0    = 20;
params.speed = 3;
params.w0    = params.speed/params.r0;

% IMU PSD convention (continuous-time intensities)
params.sig_w = deg2rad(0.2);          % rad/s/sqrt(Hz)
params.sig_a = 0.05;                  % m/s^2/sqrt(Hz)
params.Qc    = diag([params.sig_w^2*ones(3,1); params.sig_a^2*ones(3,1)]);  % 6x6

% GPS
params.sig_gps = 0.5;
params.Rgps    = (params.sig_gps^2)*eye(3);

% MAG
params.mag_dt   = 0.05;
params.mag_step = max(1, round(params.mag_dt/params.dt));
params.mag_idx  = 1:params.mag_step:params.N;

params.mI = [0.4; 0.0; 0.9]; params.mI = params.mI/norm(params.mI); % inertial mag direction
params.sig_mag = 0.05;
params.Rmag    = (params.sig_mag^2)*eye(3);

% Initial errors
params.dtheta0 = deg2rad([0;30;0]);
params.dv0     = zeros(3,1);
params.dp0     = zeros(3,1);

% Cov init
params.P0 = eye(9);

% Figure style
style.FONT_AXES   = 22;
style.FONT_LABEL  = 25;
style.FONT_LEGEND = 22;
style.LINE_WIDTH  = 2;

%% ------------------------ Precompute GT (scalar helper) ------------------------
gt.R = zeros(3,3,params.N);
gt.v = zeros(3,params.N);
gt.p = zeros(3,params.N);
for k = 1:params.N
    [Rk, vk, pk]   = circular_ground_truth(params.t(k), params.speed, params.r0);
    gt.R(:,:,k) = Rk;  gt.v(:,k) = vk;  gt.p(:,k) = pk;
end


%% Repeat each sensing configuration with paired noise realizations.
seed_list = 101:150;
all_norm = zeros(2,3,params.N,numel(seed_list),2);
results = zeros(numel(seed_list)*6,12);
row = 0;
for j = 1:numel(seed_list)
    seed = seed_list(j);
    rng(seed,'twister'); state = rng;
    res_no = run_case(false,params,gt);
    rng(state);
    res_mag = run_case(true,params,gt);
    if j == 1
        example_no = res_no; example_mag = res_mag;
    end
    cases = {res_no,res_mag};
    for scenario = 1:2
        z = cases{scenario};
        assert(all(isfinite(z.err_m(:))) && all(isfinite(z.err_i(:))));
        for q = 1:3
            idx = (q-1)*3+(1:3);
            m = vecnorm(z.err_m(idx,:),2,1);
            i = vecnorm(z.err_i(idx,:),2,1);
            all_norm(1,q,:,j,scenario) = reshape(m,1,1,[]);
            all_norm(2,q,:,j,scenario) = reshape(i,1,1,[]);
            early = params.t<=10; late = params.t>=50;
            row = row+1;
            results(row,:) = [seed scenario-1 q ...
                sqrt(mean(m(early).^2)) sqrt(mean(i(early).^2)) ...
                sqrt(mean(m.^2)) sqrt(mean(i.^2)) ...
                sqrt(mean(m(late).^2)) sqrt(mean(i(late).^2)) ...
                max(m) max(i) 1];
        end
    end
    if mod(j,10)==0, fprintf('Completed %d/%d paired seeds\n',j,numel(seed_list)); end
end
headers = {'seed','magnetometer','quantity_1attitude_2velocity_3position', ...
    'early_RMSE_MEKF','early_RMSE_Inv','full_RMSE_MEKF','full_RMSE_Inv', ...
    'late_RMSE_MEKF','late_RMSE_Inv','peak_MEKF','peak_Inv','finite'};
writetable(array2table(results,'VariableNames',headers),fullfile(outdir,'paired_run_metrics.csv'));
save(fullfile(outdir,'validated_results.mat'),'params','gt','seed_list','all_norm', ...
    'results','headers','example_no','example_mag','-v7.3');
export_comparison(all_norm,example_no,example_mag,gt,params,style,outdir);
for mag=0:1
    for q=1:3
        r=results(results(:,2)==mag & results(:,3)==q,:);
        scale=1; if q==1, scale=180/pi; end
        fprintf('MAG=%d quantity=%d early median RMSE [MEKF Inv]=[%.5f %.5f]; Inv wins %d/%d\n', ...
            mag,q,scale*median(r(:,4:5)),sum(r(:,5)<r(:,4)),numel(seed_list));
    end
end
fprintf('Figures and reproducible results saved to %s\n',outdir);

results = struct('results',results,'headers',{headers},'all_norm',all_norm, ...
    'seed_list',seed_list,'example_no',example_no,'example_mag',example_mag, ...
    'params',params,'gt',gt,'output_dir',outdir);
end

function restore_session(session)
    rng(session.rng);
    path(session.path);
    set(groot,'defaultTextInterpreter',session.textInterpreter);
    set(groot,'defaultAxesTickLabelInterpreter',session.tickInterpreter);
    set(groot,'defaultLegendInterpreter',session.legendInterpreter);
end

function res = run_case(useMag, params, gt)
    I9 = eye(9);

    % Indices
    gps_idx = params.gps_idx;
    Ngps = numel(gps_idx);

    mag_idx = params.mag_idx;
    Nmag = numel(mag_idx);

    % Measurements
    y_gps = gt.p(:,gps_idx) + params.sig_gps*randn(3,Ngps);

    % Generate MAG noise in both cases so the IMU/GPS noise stays paired.
    y_mag = zeros(3,Nmag);
    for j = 1:Nmag
        k = mag_idx(j);
        y_mag(:,j) = gt.R(:,:,k)'*params.mI + params.sig_mag*randn(3,1);
    end

    % Initial condition
    Rhat0 = so3_exp(params.dtheta0) * gt.R(:,:,1);
    vhat0 = gt.v(:,1) + params.dv0;
    phat0 = gt.p(:,1) + params.dp0;

    % MEKF state (components)
    Rhat_m = Rhat0; vhat_m = vhat0; phat_m = phat0;

    % InvEKF state (group)
    Xhat_i = Gamma(Rhat0, vhat0, phat0);

    % Covariances: params.P0 specifies the common local physical prior
    % in MEKF coordinates [theta; v-vhat; p-phat], Rhat*R' = Exp(theta).
    % Map to xi = Log(Xhat*X^{-1})^vee at the common initial estimate.
    P_m = params.P0;
    T0 = [eye(3),         zeros(3), zeros(3);
          so3_hat(vhat0), -eye(3),  zeros(3);
          so3_hat(phat0), zeros(3), -eye(3)];
    P_i = T0 * P_m * T0';   % Matching physical covariance to first order

    % Storage
    N = params.N;
    err_m = zeros(9,N);
    err_i = zeros(9,N);
    p_m_hist = zeros(3,N); v_m_hist = zeros(3,N);
    p_i_hist = zeros(3,N); v_i_hist = zeros(3,N);

    % Counters
    gps_counter = 1;
    mag_counter = 1;

    for k = 1:N

        if k > 1
            % Propagate from t_{k-1} to t_k
            [omega_true, acc_true] = imu_circular_motion(params.t(k-1), params.speed, params.r0, params.g);
            omega_m = omega_true + (params.sig_w/sqrt(params.dt))*randn(3,1);
            acc_m   = acc_true   + (params.sig_a/sqrt(params.dt))*randn(3,1);

            % Covariance propagation
            % MEKF
            A_m = [ zeros(3),               zeros(3), zeros(3);
                    so3_hat(Rhat_m*acc_m),  zeros(3), zeros(3);
                    zeros(3),               eye(3),   zeros(3) ];
            Phi_m = I9 + A_m*params.dt;

            G_m  = [ Rhat_m,    zeros(3);
                     zeros(3), -Rhat_m;
                     zeros(3),  zeros(3) ];
            P_m = Phi_m*P_m*Phi_m' + (G_m*params.Qc*G_m')*params.dt;

            % InvEKF
            Phi_i = Ak_RIEKF(params.g, params.dt);
            G_i   = Lk_RIEKF(Xhat_i, params.g, params.dt)/params.dt;
            P_i = Phi_i*P_i*Phi_i' + (G_i*params.Qc*G_i')*params.dt;

            % Both filters receive the same IMU sample.
            Xhat_i = se23_propagate(Xhat_i, omega_m, acc_m, params.g, params.dt);

            Xhat_m = se23_propagate_classical(Gamma(Rhat_m, vhat_m, phat_m), omega_m, acc_m, params.g, params.dt);
            [Rhat_m, vhat_m, phat_m] = se23_extract(Xhat_m);

        end

        [~, ~, pi] = se23_extract(Xhat_i);

        % GPS update
        if (gps_counter <= Ngps) && (k == gps_idx(gps_counter))
            yk = y_gps(:,gps_counter);
            gps_counter = gps_counter + 1;

            % MEKF: z = y - p
            z_m = yk - phat_m;
            H_m = [zeros(3), zeros(3), eye(3)];
            S_m  = H_m*P_m*H_m' + params.Rgps;
            K_m  = (P_m*H_m')/S_m;
            dx_m = K_m*z_m;

            dth = dx_m(1:3);  dv = dx_m(4:6);  dp = dx_m(7:9);
            Rhat_m = so3_exp(-dth) * Rhat_m;
            vhat_m = vhat_m + dv;
            phat_m = phat_m + dp;

            P_m = (I9 - K_m*H_m)*P_m;
            P_m = 0.5*(P_m + P_m');

            % InvEKF right-invariant residual
            z_i = pi - yk;
            H_i = [-so3_hat(yk) zeros(3) eye(3)];
            S_i = H_i*P_i*H_i' + params.Rgps;
            K_i = (P_i*H_i')/S_i;
            Xhat_i = se23_exp(-K_i*z_i) * Xhat_i;

            P_i = (I9 - K_i*H_i)*P_i;
            P_i = 0.5*(P_i + P_i');
        end

        % MAG update (optional)
        if useMag && (mag_counter <= Nmag) && (k == mag_idx(mag_counter))
            yk = y_mag(:,mag_counter);
            mag_counter = mag_counter + 1;

            % Magnetometer residual: z = Rhat*y - mI
            z_m = Rhat_m*yk - params.mI;
            Hm  = [-so3_hat(params.mI), zeros(3), zeros(3)];
            Sm  = Hm*P_m*Hm' + params.Rmag;
            Km  = (P_m*Hm')/Sm;
            dxm = Km*z_m;

            dth = dxm(1:3);
            Rhat_m = so3_exp(-dth) * Rhat_m;
            vhat_m = vhat_m + dxm(4:6);
            phat_m = phat_m + dxm(7:9);

            P_m = (I9 - Km*Hm)*P_m;
            P_m = 0.5*(P_m + P_m');

            % InvEKF
            [Ri,~,~] = se23_extract(Xhat_i);
            z_i = Ri*yk - params.mI;
            Hi  = [-so3_hat(params.mI), zeros(3), zeros(3)];
            Si  = Hi*P_i*Hi' + params.Rmag;
            Ki  = (P_i*Hi')/Si;
            Xhat_i = se23_exp(-Ki*z_i) * Xhat_i;

            P_i = (I9 - Ki*Hi)*P_i;
            P_i = 0.5*(P_i + P_i');
        end

        % True state
        Rtrue = gt.R(:,:,k);
        vtrue = gt.v(:,k);
        ptrue = gt.p(:,k);

        % Store estimates
        p_m_hist(:,k) = phat_m;  v_m_hist(:,k) = vhat_m;

        [Ri, vi, pi]  = se23_extract(Xhat_i);
        p_i_hist(:,k) = pi;      v_i_hist(:,k) = vi;

        % Exact rotation-vector errors (valid also for large angles).
        th_m = rotation_error_vector(Rhat_m * Rtrue');
        th_i = rotation_error_vector(Ri * Rtrue');
        err_m(:,k) = [th_m; (vhat_m - vtrue); (phat_m - ptrue)];
        err_i(:,k) = [th_i; (vi     - vtrue); (pi     - ptrue)];

    end

    % Pack output
    res.p_m_hist = p_m_hist; res.v_m_hist = v_m_hist;
    res.p_i_hist = p_i_hist; res.v_i_hist = v_i_hist;
    res.err_m = err_m; res.err_i = err_i;
end

%% ======================================================================
function phi = rotation_error_vector(R)
    u = 0.5*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    sn = norm(u);
    cs = max(-1, min(1, (trace(R)-1)/2));
    angle = atan2(sn, cs);
    if angle < 1e-8
        phi = u;
    elseif pi-angle < 1e-6
        [V,D] = eig(0.5*(R+R'));
        [~,j] = max(diag(D));
        axis = V(:,j);
        if axis'*u < 0, axis = -axis; end
        phi = angle*axis;
    else
        phi = (angle/sn)*u;
    end
end

%% ======================================================================

function export_comparison(all_norm,example_no,example_mag,gt,params,style,outdir)
% Error panels summarize all paired seeds. Trajectories use the first seed.
% Full-time view uses a logarithmic y-axis; startup view uses linear y-axes.
% Percentile bands and every displayed trajectory point set the axis limits.
    colors=[0.00 0.4470 0.7410; 0.8500 0.3250 0.0980];
    titles={'$\|\tilde{\theta}\|$ [rad]','$\|\tilde v\|$ [m/s]','$\|\tilde p\|$ [m]','3D Trajectory'};
    names={'No Compass','With Compass'};
    t=params.t;
    med=zeros(2,3,numel(t),2); lo=med; hi=med;
    for r=1:2
        for q=1:3
            for f=1:2
                y=squeeze(all_norm(f,q,:,:,r));
                med(f,q,:,r)=median(y,2);
                lo(f,q,:,r)=point_percentile(y,0.10);
                hi(f,q,:,r)=point_percentile(y,0.90);
            end
        end
    end
    cases={example_no,example_mag};
    for view_id=1:2
        if view_id==1
            window=[0 params.T]; tag='Example_circular_log'; scale='log';

        else
            window=[0 10]; tag='Example_circular'; scale='linear';

        end
        visible=t>=window(1) & t<=window(2);
        fig=figure('Color',[1 .9 .9],'Units','pixels','Position',[50 50 1800 800]);
        tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','loose');

        upper=zeros(1,3); lower=zeros(1,3);
        for q=1:3
            yy=reshape(hi(:,q,visible,:),[],1);
            upper(q)=max(yy)*1.10;
            if strcmp(scale,'log')
                yy=reshape(lo(:,q,visible,:),[],1); yy=yy(yy>1e-12);
                lower(q)=min(yy)/1.3;
                upper(q)=upper(q)*1.2;
            end
        end
        for r=1:2
            for q=1:3
                ax=nexttile(tl,(r-1)*4+q); hold(ax,'on');
                set(ax,'YScale',scale,'Color','w','FontSize',style.FONT_AXES,'LineWidth',1);
                handles=gobjects(1,2);
                for f=1:2
                    a=squeeze(lo(f,q,:,r)); b=squeeze(hi(f,q,:,r)); m=squeeze(med(f,q,:,r));
                    if strcmp(scale,'log')
                        % Omit zero and floating-point roundoff at exact initialization.
                        a(a<=1e-12)=NaN; b(b<=1e-12)=NaN; m(m<=1e-12)=NaN;
                    end
                    fill(ax,[t;flipud(t)],[a;flipud(b)],colors(f,:), ...
                        'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
                    handles(f)=plot(ax,t,m,'Color',colors(f,:),'LineWidth',style.LINE_WIDTH);
                end
                xlim(ax,window); ylim(ax,[lower(q),upper(q)]);
                if strcmp(scale,'log')
                    ax.YTick=10.^(ceil(log10(lower(q))):floor(log10(upper(q))));
                end
                grid(ax,'on'); box(ax,'on'); ax.YMinorGrid='off';
                if r==1, title(ax,titles{q},'FontSize',style.FONT_LABEL,'Interpreter','latex'); end
                xlabel(ax,'$t$ [s]','FontSize',style.FONT_LABEL,'Interpreter','latex');
                if r==1, ax.XTick=[]; ax.XTickLabel={}; xlabel(ax,''); end
                if q==1, ylabel(ax,names{r},'FontSize',style.FONT_LABEL,'Interpreter','latex'); end
                if r==1 && q==1
                    lg=legend(ax,handles,{'MEKF','InvEKF'},'Orientation','horizontal','FontSize',style.FONT_LEGEND);
                    lg.Layout.Tile='north';
                end
            end
            ax=nexttile(tl,r*4); hold(ax,'on'); z=cases{r};
            plot3(ax,z.p_m_hist(1,:),z.p_m_hist(2,:),z.p_m_hist(3,:),'Color',colors(1,:),'LineWidth',style.LINE_WIDTH);
            plot3(ax,z.p_i_hist(1,:),z.p_i_hist(2,:),z.p_i_hist(3,:),'Color',colors(2,:),'LineWidth',style.LINE_WIDTH);
            plot3(ax,gt.p(1,:),gt.p(2,:),gt.p(3,:),'k--','LineWidth',2.5);
            set(ax,'FontSize',style.FONT_AXES,'Color','w'); grid(ax,'on'); box(ax,'on'); view(ax,45,25);
            xlabel(ax,''); ylabel(ax,''); zlabel(ax,'');
            ax.YTick=[]; ax.YTickLabel={};
            if r==1, ax.XTick=[]; ax.XTickLabel={}; end
            zextent=max(2,ceil(max(abs([z.p_m_hist(3,:),z.p_i_hist(3,:),gt.p(3,:)]))));
            zlim(ax,[-zextent zextent]);
            if r==1, title(ax,titles{4},'FontSize',style.FONT_LABEL,'Interpreter','latex'); end
            % Auto limits include every estimated and true trajectory point.

        end
        exportgraphics(fig,fullfile(outdir,[tag '.pdf']),'ContentType','vector','Padding',12);
        exportgraphics(fig,fullfile(outdir,[tag '.png']),'Resolution',170,'Padding',12);
    end
end

function p=point_percentile(y,fraction)
% Linear interpolation of order statistics; no toolbox dependency.
    ys=sort(y,2); n=size(y,2); u=1+fraction*(n-1);
    a=floor(u); b=ceil(u);
    p=ys(:,a)*(b-u)+ys(:,b)*(u-a);
    if a==b, p=ys(:,a); end
end

