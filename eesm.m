% =========================================================================
%  EESM_fit_script_tensile_v1.m
%  Fits EESM parameters to experimental TENSILE test data (dogbone/probe).
%  PLAIN SCRIPT — all variables visible in Workspace after run.
%
%  INPUT MATRIX FORMAT (bgir):
%    column 1 -> anything (not used)
%    column 2 -> displacement/extension, mm  (POSITIVE = stretching)
%    column 3 -> force, kN  (POSITIVE = tensile pulling force)
%
%  Set L0 (mm), A (m^2), and E (Pa, independently measured linear
%  modulus) in Section 2.
%
%  =====================  WHY THIS VERSION (tensile)  =====================
%  Every previous EESM_fit_script (v3, v5) assumed COMPRESSION:
%       lambda = 1 - eps        (stretch DECREASES as displacement grows)
%  For a genuine TENSILE test, stretch INCREASES instead:
%       lambda = 1 + eps
%  This isn't just a sign flip — the LOCKING behavior also moves to the
%  opposite side of lambda=1:
%    - In compression, I1(lambda)=lambda^2+2/lambda blows up as
%      lambda -> 0 (the SMALL-lambda root of I1=target is the locking
%      point), which is what eps_max_safe_of_N solved for previously.
%    - In tension, I1(lambda) blows up as lambda -> infinity instead (the
%      LARGE-lambda root is the locking point). eps_max_safe_of_N has
%      been RE-DERIVED below (Newton-Raphson started from lambda=2 instead
%      of lambda=0.3) to find this large-lambda root correctly.
%
%  Both roots are real and distinct: I1(lambda) has a single minimum at
%  lambda=1 (I1=3) and increases monotonically as lambda moves away from
%  1 in EITHER direction — so I1(lambda)=target always has exactly one
%  root below 1 (compression locking) and one root above 1 (tension
%  locking). Using the wrong root silently gives a physically meaningless
%  locking-strain bound, which is why this needed a real re-derivation,
%  not just changing a sign.
%
%  STABLE DESIGN RETAINED FROM v5 (not v4's unstable derived-A1 approach):
%  A1 is a BOUNDED free parameter (no division-based amplification risk).
%  E-consistency enforced via a strong soft penalty (w_E=2000), not a
%  hard constraint.
%
%  ALL OTHER FIXES RETAINED:
%  FIX-F: mu anchored near E/3.
%  FIX-I: A2 magnitude penalized relative to mu*N.
%  FIX-J: Cancellation check swept across the entire strain range.
% =========================================================================

fprintf('=============================================================\n');
fprintf('   EESM PARAMETER FITTING — TENSILE TEST  (lambda=1+eps)      \n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 1 — DATA
%  Replace with your tensile experimental matrix [idx, disp_mm, force_kN]
% =========================================================================

fprintf('--- TENSILE TEST DATA ---\n');
fprintf('  Expected format: [idx, disp_mm, force_kN]  (both POSITIVE for tension)\n');
fprintf('  You can either:\n');
fprintf('    (a) type an existing workspace variable name\n');
fprintf('    (b) type the matrix directly, e.g.  [0 0 0; 1 0.1 0.05; 2 0.2 0.11]\n');
bgir = input('  bgir = ');
while ~isnumeric(bgir) || size(bgir,2)<3
    fprintf('  *** bgir must be a numeric matrix with at least 3 columns.\n');
    bgir = input('  Re-enter bgir = ');
end

% =========================================================================
%  SECTION 2 — SPECIMEN GEOMETRY & LINEAR MODULUS  (edit/enter these three)
% =========================================================================

L0 = input('  L0 (mm)   initial gauge length of the tensile probe: ');
A  = input('  A (m^2)   cross-sectional area of the gauge section: ');
E  = input('  E (Pa)    independently measured linear (small-strain) modulus: ');

fprintf('  L0 = %.4f mm\n', L0);
fprintf('  A  = %.4e m^2\n', A);
fprintf('  E  = %.4e Pa   (anchor for mu/f/A1 consistency, strong penalty)\n\n', E);

% =========================================================================
%  SECTION 3 — CONVERT TO STRAIN AND FORCE  (TENSION: both positive)
% =========================================================================

disp_mm  = bgir(:,2);
force_kN = bgir(:,3);

eps_data = disp_mm / L0;       % nominal strain, POSITIVE = tension/stretching
F_data   = force_kN * 1000;    % kN -> N  (POSITIVE = tensile force)

abs_F_data = abs(F_data);

fprintf('  Data points   : %d\n', numel(eps_data));
fprintf('  Strain range  : %.4f to %.4f (%.1f%% to %.1f%%)\n', ...
        min(eps_data), max(eps_data), min(eps_data)*100, max(eps_data)*100);
fprintf('  |Force| range : %.3f to %.3f N\n', min(abs_F_data), max(abs_F_data));
fprintf('  lambda range  : %.4f up to %.4f   (TENSION: lambda=1+eps, increasing)\n\n', ...
        1+min(eps_data), 1+max(eps_data));

if min(eps_data) > 1e-6
    fprintf('  Adding anchor (eps=0, F=0) — data does not start at zero.\n\n');
    eps_data    = [0;    eps_data];
    abs_F_data  = [0;    abs_F_data];
else
    fprintf('  Data starts at eps=0. Good.\n\n');
end

eps_max_data = max(eps_data);

% =========================================================================
%  SECTION 3.5 — DATA-DEPENDENT LOWER BOUND ON N   (TENSION locking branch)
%
%  eps_max_safe(N) = TENSILE strain at which the chain network reaches
%  95% of its locking stretch. Unlike compression, this is found from the
%  LARGE-lambda root of I1(lambda)=target (see header note above) — the
%  Newton-Raphson search below starts from lambda=2 (tension side),
%  not lambda=0.3 (compression side) as in the original scripts.
% =========================================================================

margin = 1.3;   % locking strain must exceed max tested strain by 30%
target_eps_max_safe = margin * eps_max_data;

fprintf('  Solving for minimum N such that eps_max_safe >= %.1fx eps_max_data ...\n', margin);
N_lo = 0.3; N_hi = 200;   % bisection bracket
if eps_max_safe_of_N_tension(N_lo) >= target_eps_max_safe
    N_min_physical = N_lo;
else
    for bi = 1:60
        N_mid = 0.5*(N_lo+N_hi);
        if eps_max_safe_of_N_tension(N_mid) < target_eps_max_safe
            N_lo = N_mid;
        else
            N_hi = N_mid;
        end
    end
    N_min_physical = N_hi;
end
fprintf('  -> N_min_physical = %.4f  (gives eps_max_safe = %.4f, data max = %.4f)\n\n', ...
        N_min_physical, eps_max_safe_of_N_tension(N_min_physical), eps_max_data);

% =========================================================================
%  SECTION 4 — BOUNDS   (A1 bounded/stable, v5 design)
%
%       N                    f       A1          A2          mu
% =========================================================================
mu_center = E/3;   % Arruda-Boyce small-strain limit for incompressible solid

lb = [  10           0.05     -5e5            -5e5        0.1*mu_center ];
ub = [  40          0.99      5e12         5e12         2.0*mu_center ];

fprintf('  Bounds:\n');
fprintf('    N  in [%.3f, %.3f]   (physical lower bound from data, TENSION branch)\n', lb(1), ub(1));
fprintf('    f  in [%.3f, %.3f]\n', lb(2), ub(2));
fprintf('    A1 in [%.3e, %.3e] Pa   (bounded, stable — no division instability)\n', lb(3), ub(3));
fprintf('    A2 in [%.3e, %.3e] Pa   (also soft-penalized vs mu*N)\n', lb(4), ub(4));
fprintf('    mu in [%.3e, %.3e] Pa   (anchored near E/3 = %.3e)\n\n', lb(5), ub(5), mu_center);

eps_max_check = eps_max_data;

% =========================================================================
%  SECTION 5 — GLOBAL SEARCH  (Genetic Algorithm or multistart fallback)
% =========================================================================

fprintf('Stage 1: Global search ...\n');
objfun = @(p) eesm_objective_tension(p, eps_data, abs_F_data, A, eps_max_check, E);

if exist('ga','file')==2 || exist('ga','file')==6
    ga_opts = optimoptions('ga', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 400, ...
        'FunctionTolerance', 1e-10, ...
        'Display', 'iter', ...
        'UseParallel', false);
    [p_ga, cost_ga] = ga(objfun, 5, [], [], [], [], lb, ub, [], ga_opts);
else
    fprintf('  ga() not found. Using multistart fmincon fallback.\n\n');
    [p_ga, cost_ga] = multistart_fmincon(objfun, lb, ub, 60);
end

fprintf('\n  Global search cost = %.6e\n', cost_ga);
c_ga = compute_c(p_ga(1), p_ga(5));
fprintf('  N=%.4f  f=%.4f  A1=%.4e  A2=%.4e  mu=%.4e  c=%.4e\n\n', ...
        p_ga(1),p_ga(2),p_ga(3),p_ga(4),p_ga(5), c_ga);

% =========================================================================
%  SECTION 6 — LOCAL REFINEMENT  (bounded fmincon)
% =========================================================================

fprintf('Stage 2: Local refinement (bounded fmincon) ...\n');
fmc_opts = optimoptions('fmincon', ...
    'Display','off', ...
    'MaxIterations', 3000, ...
    'MaxFunctionEvaluations', 20000, ...
    'OptimalityTolerance', 1e-12, ...
    'StepTolerance', 1e-12);
[p_final, cost_final] = fmincon(objfun, p_ga, [], [], [], [], lb, ub, [], fmc_opts);
fprintf('  fmincon cost = %.6e\n\n', cost_final);

% Unpack
N_fit  = p_final(1);
f_fit  = p_final(2);
A1_fit = p_final(3);
A2_fit = p_final(4);
mu_fit = p_final(5);
c_fit  = compute_c(N_fit, mu_fit);

fprintf('=============================================================\n');
fprintf('  FINAL FITTED EESM PARAMETERS  (from TENSILE data)\n');
fprintf('=============================================================\n');
fprintf('  mu = %.6e  Pa\n', mu_fit);
fprintf('  N  = %.6f\n',     N_fit);
fprintf('  A1 = %.6e  Pa\n', A1_fit);
fprintf('  A2 = %.6e  Pa\n', A2_fit);
fprintf('  f  = %.6f\n',     f_fit);
fprintf('  c  = %.6e  Pa   [computed: W_iso(lam=1)=0]\n', c_fit);
fprintf('\n  These describe the MATERIAL — usable for ANY geometry made of\n');
fprintf('  the same resin, unlike a fit against a specific specimen''s\n');
fprintf('  structural (geometry-diluted) response.\n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 7 — VALIDATION
% =========================================================================

F_pred_data = F_model_abs_tension(eps_data, A, N_fit, f_fit, A1_fit, A2_fit, mu_fit);
ss_res = sum((abs_F_data - F_pred_data).^2);
ss_tot = sum((abs_F_data - mean(abs_F_data)).^2);
r2 = 1 - ss_res / max(ss_tot, 1e-12);

eps_fine = linspace(0.001, eps_max_data*1.05, 400)';
F_fine   = F_model_abs_tension(eps_fine, A, N_fit, f_fit, A1_fit, A2_fit, mu_fit);
dF       = diff(F_fine);
is_mono  = all(dF >= -1e-3);
min_F    = min(F_fine);

% ---- sigma sign check at maximum strain (TENSION: lambda=1+eps) --------
eps_check_sig = eps_max_data;
lam_check     = 1 + eps_check_sig;   % TENSION: + not -
sig_check     = eesm_sigma_check(lam_check, mu_fit, N_fit, ...
                                  A1_fit, A2_fit, c_fit, f_fit);
sigma_ok      = isfinite(sig_check) && sig_check > 0;

% ---- F_eq28 check at maximum strain --------------------------------------
WT_check  = W_T_scalar(lam_check, N_fit, f_fit, A1_fit, A2_fit, mu_fit);
F_eq28    = 2 * A * WT_check / max(eps_check_sig, 1e-9);
F_target  = max(abs_F_data);
force_ok  = isfinite(F_eq28) && abs(F_eq28/max(F_target,1e-9)-1) < 0.30;

% ---- Cancellation ratio SWEPT ACROSS FULL RANGE (TENSION) ---------------
eps_sweep = linspace(max(eps_max_data*0.05,1e-3), eps_max_data, 40)';
cancel_sweep = zeros(size(eps_sweep));
for kk = 1:numel(eps_sweep)
    lam_k = 1 + eps_sweep(kk);   % TENSION: + not -
    [Wiso_k, Waniso_k] = W_components(lam_k, mu_fit, N_fit, A1_fit, A2_fit, c_fit, f_fit);
    WT_k = (1-f_fit)*Wiso_k + f_fit*Waniso_k;
    cancel_sweep(kk) = abs(WT_k) / max(abs(Wiso_k), 1e-12);
end
cancel_ratio_min = min(cancel_sweep);
[cancel_worst_val, cancel_worst_idx] = min(cancel_sweep);
cancel_worst_eps = eps_sweep(cancel_worst_idx);

% ---- E-consistency check -------------------------------------------------
E_eff       = 3*(1-f_fit)*mu_fit + 2*f_fit*A1_fit;
E_mismatch  = abs(E_eff - E)/E;
E_ok        = E_mismatch < 0.15;

% ---- A2 relative-magnitude check ----------------------------------------
A2_scale    = 2*mu_fit*N_fit;
A2_ratio    = abs(A2_fit)/max(A2_scale,1e-9);
A2_ok       = A2_ratio <= 1.0;

fprintf('=============================================================\n');
fprintf('  FIT VALIDATION\n');
fprintf('=============================================================\n');
fprintf('  R^2                     = %.5f\n', r2);
fprintf('  Monotonic (up to eps_max) = %d\n', is_mono);
fprintf('  Min |F| over range      = %.3f N\n', min_F);

fprintf('\n  --- Physical checks at eps = %.4f (max strain, TENSION) ---\n', eps_check_sig);
fprintf('  lam (=1+eps)            = %.6f\n', lam_check);
fprintf('  sigma (lam*dW/dlam)     = %.4f Pa\n', sig_check);
fprintf('  F via Eq.28             = %.4f N\n', F_eq28);
fprintf('  F from data (max)       = %.4f N\n', F_target);

fprintf('\n  --- Cancellation, swept across full strain range ---\n');
fprintf('  Worst |WT|/|Wiso|        = %.2f%%   at eps = %.4f\n', ...
        cancel_worst_val*100, cancel_worst_eps);

fprintf('\n  --- E-consistency check (w_E=2000) ---\n');
fprintf('  3*(1-f)*mu + 2*f*A1     = %.4e Pa\n', E_eff);
fprintf('  Measured E              = %.4e Pa\n', E);
fprintf('  Mismatch                = %.2f%%\n', E_mismatch*100);

fprintf('\n  --- A2 relative magnitude check ---\n');
fprintf('  |A2| / (2*mu*N)         = %.3f  (want <= 1.0)\n', A2_ratio);

fprintf('\n  SIGMA CHECK  : %s\n', tf_str(sigma_ok, 'PASS — sigma > 0 (EESM valid in tension)', ...
    'FAIL — sigma <= 0. RE-FIT required (try lower f upper bound).'));
fprintf('  FORCE CHECK  : %s\n', tf_str(force_ok, 'PASS — F_eq28 within 30%% of data', ...
    'FAIL — F_eq28 far from data. Parameters may not match material.'));
fprintf('  CANCEL CHECK : %s\n', tf_str(cancel_ratio_min>0.4, ...
    'PASS — cancellation stays above 40%% everywhere in range', ...
    'WARN — cancellation drops below 40%% somewhere. Et unreliable there.'));
fprintf('  E CHECK      : %s\n', tf_str(E_ok, ...
    'PASS — model small-strain modulus matches measured E within 15%%', ...
    'WARN — model disagrees with measured E by >15%%. Try raising w_E further.'));
fprintf('  A2 CHECK     : %s\n', tf_str(A2_ok, ...
    'PASS — A2 within its physical scale (<= 2*mu*N)', ...
    'WARN — A2 dominates the network term.'));

if ~sigma_ok
    fprintf('\n  *** ACTION: re-run fitting with lower f upper bound.\n');
    fprintf('  *** Current f = %.4f (upper bound was %.2f)\n', f_fit, ub(2));
    fprintf('  *** Try setting ub(2) = %.2f and re-running.\n', min(f_fit-0.1, 0.70));
end
if ~E_ok
    fprintf('\n  *** ACTION: fit disagrees with measured E by %.1f%% even with\n', E_mismatch*100);
    fprintf('  *** w_E=2000. Try raising w_E further (5000-10000) in\n');
    fprintf('  *** eesm_objective_tension, or double-check E was measured/entered\n');
    fprintf('  *** correctly (this should be the TRUE material E, e.g. from a\n');
    fprintf('  *** solid coupon''s linear region — never a lattice''s apparent E).\n');
end
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 7.5 — SAVE FITTED PARAMETERS
% =========================================================================

cancel_ok = cancel_ratio_min > 0.4;

eesm_mu = mu_fit;
eesm_N  = N_fit;
eesm_A1 = A1_fit;
eesm_A2 = A2_fit;
eesm_f  = f_fit;
eesm_c  = c_fit;

eesm_sigma_ok   = sigma_ok;
eesm_E_ok       = E_ok;
eesm_cancel_ok  = cancel_ok;
eesm_A2_ok      = A2_ok;
eesm_fit_all_ok = sigma_ok && E_ok && cancel_ok && A2_ok;

eesm_fit_eps_max   = eps_max_data;
eesm_fit_E_used    = E;
eesm_fit_A_used    = A;
eesm_fit_L0_used   = L0;
eesm_fit_r2        = r2;
eesm_fit_timestamp = datestr(now);
eesm_fit_source    = 'tensile';   % NEW: records that this came from a tensile test

out_file = 'eesm_fitted_params.mat';
save(out_file, 'eesm_mu','eesm_N','eesm_A1','eesm_A2','eesm_f','eesm_c', ...
               'eesm_sigma_ok','eesm_E_ok','eesm_cancel_ok','eesm_A2_ok','eesm_fit_all_ok', ...
               'eesm_fit_eps_max','eesm_fit_E_used','eesm_fit_A_used','eesm_fit_L0_used', ...
               'eesm_fit_r2','eesm_fit_timestamp','eesm_fit_source');

fprintf('=============================================================\n');
fprintf('  SAVED: %s\n', out_file);
fprintf('=============================================================\n');
fprintf('  The FEM solvers will load this file automatically when you\n');
fprintf('  select the EESM model — no manual parameter entry needed.\n');
if eesm_fit_all_ok
    fprintf('  All validation checks PASSED. Safe to use directly.\n');
else
    fprintf('  *** NOT ALL CHECKS PASSED — the solver will warn/block on load. ***\n');
end
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 8 — PLOTS
% =========================================================================

figure('Name','EESM Fit (Tensile)','NumberTitle','off');
plot(eps_data, abs_F_data, 'ko', 'MarkerFaceColor','k', ...
     'MarkerSize',6, 'DisplayName','Experimental |F| (tensile)'); hold on;
plot(eps_fine, F_fine, 'r-', 'LineWidth',2, 'DisplayName','EESM fit');
xline(eps_max_data,'--','Color',[.5 .5 .5],'Label','data range end');
xlabel('Strain \epsilon'); ylabel('|Force| (N)');
title(sprintf('EESM tensile fit  R^2=%.4f  mono=%d  sigma=%s  E-match=%s', ...
    r2, is_mono, tf_str(sigma_ok,'YES','NO'), tf_str(E_ok,'YES','NO')));
legend('Location','northwest'); grid on;

lam_range = linspace(0.999, 1+eps_max_data*1.1, 200);   % TENSION: >1, not <1
sig_range  = zeros(size(lam_range));
for kk=1:numel(lam_range)
    sig_range(kk) = eesm_sigma_check(lam_range(kk), mu_fit, N_fit, ...
                                      A1_fit, A2_fit, c_fit, f_fit);
end
figure('Name','sigma vs lambda (Tensile)','NumberTitle','off');
plot(lam_range, sig_range/1e3, 'b-','LineWidth',2); hold on;
yline(0,'--r','LineWidth',1.5,'Label','sigma=0');
xline(lam_check,'--k','Label',sprintf('lam at max eps=%.3f',lam_check));
xlabel('\lambda'); ylabel('\sigma (kPa)');
title('Cauchy stress vs stretch (must be positive for tension \lambda>1)');
grid on;

figure('Name','Cancellation ratio vs strain (Tensile)','NumberTitle','off');
plot(eps_sweep, cancel_sweep*100, 'm-','LineWidth',2); hold on;
yline(40,'--r','LineWidth',1.2,'Label','40% floor');
xlabel('Strain \epsilon'); ylabel('|WT|/|Wiso|  (%)');
title('Cancellation ratio across full strain range');
grid on;

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

% ------------------------------------------------------------------
%  eps_max_safe_of_N_tension — TENSION locking-strain estimate.
%  Finds the LARGE-lambda root of I1(lambda)=target (lambda>1), unlike
%  the compression version which finds the SMALL-lambda root (lambda<1).
%  Newton-Raphson started from lambda=2 (tension side).
% ------------------------------------------------------------------
function eps_max_safe = eps_max_safe_of_N_tension(N)
    tI1 = 3*N*0.95^2; ls = 2.0;
    for kk=1:50
        fl  = ls^2 + 2/ls - tI1;
        dfl = 2*ls - 2/ls^2;
        ls  = max(1.001, ls - fl/dfl);
        if abs(fl)<1e-12; break; end
    end
    eps_max_safe = ls - 1;   % TENSION: lambda = 1+eps, so eps = lambda-1
end

% ------------------------------------------------------------------
%  compute_c — c from W_iso(lam=1)=0 condition (unchanged — this
%  condition is symmetric about lambda=1, so it's the same for tension
%  and compression).
% ------------------------------------------------------------------
function c = compute_c(N, mu)
    lam_r0 = 1/sqrt(N);
    beta0  = 3*lam_r0 / (1 - lam_r0^3);
    c = -mu*(N*(beta0*lam_r0 + log(beta0/sinh(beta0))) - log(beta0/lam_r0));
end

% ------------------------------------------------------------------
%  W_T_scalar — WT at a single lam value (unchanged — takes lam
%  directly, caller supplies lam=1+eps for tension)
% ------------------------------------------------------------------
function WT = W_T_scalar(lam, N, f, A1, A2, mu)
    c    = compute_c(N, mu);
    lam  = min(max(lam, 1e-4), 200);   % widened upper clamp for tension
    lam2 = 1/sqrt(lam);
    I1   = lam^2 + 2*lam2^2;
    I3   = (lam*lam2*lam2)^2;
    lam_r = min(sqrt(max(I1/(3*N),1e-9)), 0.999998);
    beta  = min(max(3*lam_r/(1-lam_r^3), 1e-6), 600);
    W_iso   = mu*(N*(beta*lam_r + log(beta/sinh(beta))) - log(beta/lam_r)) + c;
    W_aniso = (A1/3)*(I1-3) + (A2/9)*(I1-3)^2 - (A1/3)*log(max(I3,1e-12));
    WT = (1-f)*W_iso + f*W_aniso;
end

% ------------------------------------------------------------------
%  W_components — return Wiso and Waniso separately (for diagnostics)
% ------------------------------------------------------------------
function [Wiso, Waniso] = W_components(lam, mu, N, A1, A2, c, f) %#ok<INUSL>
    lam  = min(max(lam, 1e-4), 200);
    lam2 = 1/sqrt(lam);
    I1   = lam^2 + 2*lam2^2;
    I3   = (lam*lam2*lam2)^2;
    lam_r = min(sqrt(max(I1/(3*N),1e-9)), 0.999998);
    beta  = min(max(3*lam_r/(1-lam_r^3),1e-6), 600);
    Wiso   = mu*(N*(beta*lam_r + log(beta/sinh(beta))) - log(beta/lam_r)) + c;
    Waniso = (A1/3)*(I1-3) + (A2/9)*(I1-3)^2 - (A1/3)*log(max(I3,1e-12));
end

% ------------------------------------------------------------------
%  W_T_fun — vectorised WT over lam vector (used in F_model_abs_tension)
% ------------------------------------------------------------------
function WT = W_T_fun(lam, N, f, A1, A2, mu)
    c    = compute_c(N, mu);
    lam  = min(max(lam(:), 1e-4), 200);
    lam2 = 1./sqrt(lam);
    I1   = lam.^2 + 2*lam2.^2;
    I3   = (lam.*lam2.*lam2).^2;
    lam_r = min(sqrt(max(I1/(3*N), 1e-9)), 0.999998);
    beta  = min(max(3*lam_r./(1-lam_r.^3), 1e-6), 600);
    W_iso   = mu*(N*(beta.*lam_r + log(beta./sinh(beta))) - log(beta./lam_r)) + c;
    W_aniso = (A1/3)*(I1-3) + (A2/9)*(I1-3).^2 - (A1/3)*log(max(I3,1e-12));
    WT = (1-f)*W_iso + f*W_aniso;
end

% ------------------------------------------------------------------
%  F_model_abs_tension — |F| prediction (vectorised over eps), TENSION:
%  lam = 1+eps, not 1-eps.
% ------------------------------------------------------------------
function F = F_model_abs_tension(eps, A, N, f, A1, A2, mu)
    eps      = eps(:);
    lam      = 1 + eps;   % TENSION: + not -
    WT       = W_T_fun(lam, N, f, A1, A2, mu);
    eps_safe = max(abs(eps), 1e-6);
    F        = abs(2*A*WT ./ eps_safe);
end

% ------------------------------------------------------------------
%  eesm_sigma_check — Cauchy stress lam*dW/dlam (scalar, for validation).
%  Takes lam directly — caller supplies lam=1+eps for tension.
% ------------------------------------------------------------------
function s = eesm_sigma_check(lam, mu, N, A1, A2, c, f)
    dh = 1e-7;
    wp = W_T_scalar(lam+dh, N, f, A1, A2, mu);
    wm = W_T_scalar(lam-dh, N, f, A1, A2, mu);
    if ~isfinite(wp)||~isfinite(wm); s=NaN; return; end
    s = lam*(wp-wm)/(2*dh);
end

% ------------------------------------------------------------------
%  eesm_objective_tension — SSE + monotonicity/positivity + physical
%  penalties. p = [N, f, A1, A2, mu]. TENSION: lam=1+eps throughout.
% ------------------------------------------------------------------
function cost = eesm_objective_tension(p, eps_data, abs_F_data, A, eps_max_check, E)
    N=p(1); f=p(2); A1=p(3); A2=p(4); mu=p(5);

    F_pred = F_model_abs_tension(eps_data, A, N, f, A1, A2, mu);
    if any(~isfinite(F_pred)); cost=1e12; return; end
    resid = (F_pred - abs_F_data(:)) / (max(abs_F_data)+1e-9);
    sse   = sum(resid.^2);

    % Monotonicity and positivity over data range
    eps_check = linspace(0.005, eps_max_check, 80)';
    F_check   = F_model_abs_tension(eps_check, A, N, f, A1, A2, mu);
    if any(~isfinite(F_check)); cost=1e12; return; end
    dF  = diff(F_check);
    pen = 100*(sum(max(-dF,0).^2) + sum(max(-F_check,0).^2));

    % Sigma sign + cancellation penalty, TENSION: lam=1+eps
    c_now      = compute_c(N, mu);
    eps_sweep  = linspace(max(eps_max_check*0.05,1e-3), eps_max_check, 15)';
    sig_pen    = 0;
    cancel_pen = 0;
    for kk = 1:numel(eps_sweep)
        lam_k = 1 + eps_sweep(kk);   % TENSION: + not -
        s_k   = eesm_sigma_check(lam_k, mu, N, A1, A2, c_now, f);
        if ~isfinite(s_k) || s_k<=0
            sig_pen = sig_pen + 1e4;
        end
        [Wiso_k, Waniso_k] = W_components(lam_k, mu, N, A1, A2, c_now, f);
        WT_k    = (1-f)*Wiso_k + f*Waniso_k;
        cancel_k = abs(WT_k)/max(abs(Wiso_k),1e-12);
        if cancel_k < 0.4
            cancel_pen = cancel_pen + 50*(0.4-cancel_k)^2;
        end
    end

    % E-consistency penalty (STRONG weight, v5 design)
    E_eff    = 3*(1-f)*mu + 2*f*A1;
    E_resid  = (E_eff - E)/max(E,1e-9);
    w_E      = 2000;
    E_pen    = w_E*E_resid^2;

    % A2 magnitude penalty relative to mu*N
    A2_scale = 2*mu*N;
    A2_ratio = abs(A2)/max(A2_scale,1e-9);
    A2_pen   = 20*max(A2_ratio-1, 0)^2;

    cost = sse + pen + sig_pen + cancel_pen + E_pen + A2_pen;
end

% ------------------------------------------------------------------
%  multistart_fmincon — fallback when ga is unavailable
% ------------------------------------------------------------------
function [p_best, cost_best] = multistart_fmincon(objfun, lb, ub, n)
    opts = optimoptions('fmincon','Display','off', ...
                        'MaxIterations',1000,'MaxFunctionEvaluations',5000);
    p_best=[]; cost_best=Inf; rng(42);
    for k=1:n
        p0 = lb + (ub-lb).*rand(1,5);
        try
            [pk,ck] = fmincon(objfun,p0,[],[],[],[],lb,ub,[],opts);
        catch; continue; end
        if ck<cost_best; cost_best=ck; p_best=pk; end
        if mod(k,10)==0
            fprintf('    multistart %d/%d  best=%.4e\n',k,n,cost_best);
        end
    end
end

% ------------------------------------------------------------------
%  tf_str — ternary string helper
% ------------------------------------------------------------------
function s = tf_str(cond, s_true, s_false)
    if cond; s=s_true; else; s=s_false; end
end