% symbolicDerivation.m
%
% 1) Symbolic setup: define k1, k2, ksw, D0, t and assume positivity.
%    Define initial conditions. 
%    Assemble the 4×4 rate matrix A under pseudo–first‑order conditions.
% 2) Compute c4(t) = W(t) via matrix exponential, c4_expm = (expm(A*t)*c0).
% 3) Compute c4(t) via Laplace -> partial fractions -> inverse Laplace -> c4_lap.
% 4) Numerically compare expm vs inverse‑Laplace for k1=k2=1, ksw=1.001, D0=100.
% 5) Compute degenerate limit k2->k_sw symbolically and overlay on same plot.
%
% Note: expm(X) uses MATLAB’s scaling‑and‑squaring algorithm.  
%       If X=V*D/V is diagonalizable, expm(X) = V*diag(exp(diag(D)))/V.
%
% Oliver Steinbock (FSU, 2025)

%% 1) Symbolic setup (P0 = 1)
syms k1 k2 ksw D0 t real
assumeAlso([k1 k2 ksw D0], 'positive')

% Build the 4×4 rate-matrix A under pseudo–first-order ([D]=D0=const)
A = [ -3*k1*D0,      0,        0,     0;
         k1*D0, -k2*D0,        0,     0;
          0,     k2*D0,  -ksw*D0,     0;
          0,         0,   ksw*D0,     0 ];

% Initial condition at t=0: P0 = 1, others = 0
c0 = [1; 0; 0; 0];

%% 2) Compute W(t) via matrix exponential
c_expm = expm(A*t) * c0;
c4_expm = simplify(c_expm(4));
fprintf('\n-----------\n')
disp('1) expm-derived W(t) =')
disp(c4_expm)

%% 3) Compute f(t) via Laplace -> partial fraction -> inverse Laplace
syms s
F    = laplace(c4_expm, t, s);
Fpf  = partfrac(F, s);
c4_lap = simplify(ilaplace(Fpf, s, t));
disp('2) inverse-Laplace (partial-fraction) W(t) =')
disp(c4_lap)

% Output the inverse-Laplace form as LaTeX
latex_str = latex(c4_lap);
disp('3) LaTeX of inverse-Laplace form:')
disp(latex_str)
fprintf('\n')

%% 4) Numeric plot comparison for example values
k1_val  = 1.0;
k2_val  = 1.0;
ksw_val = 1.001;
D0_val  = 100;

t_vals = linspace(0,0.15,300);
c4_expm_vals = double(subs(c4_expm, {t,k1,k2,ksw,D0}, {t_vals,k1_val,k2_val,ksw_val,D0_val}));
c4_lap_vals  = double(subs(c4_lap,  {t,k1,k2,ksw,D0}, {t_vals,k1_val,k2_val,ksw_val,D0_val}));

figure;
plot(t_vals, c4_expm_vals, 'k-', 'LineWidth',2.5);
hold on;
plot(t_vals, c4_lap_vals,  'r--', 'LineWidth',1.5);
xlabel('t');
ylabel('c_4(t)');
title(sprintf('c_4(t) for P_0=1, k1=%.3f, k2=%.3f, k_{sw}=%.3f, D_0=%g', ...
    k1_val,k2_val,ksw_val,D0_val));
% legend('expm(A t)','inverse-Laplace','Location','best');
grid on;

%% 5) l’Hôpital limit for k2 = ksw != 3 k1
temp = limit(c4_expm, k2, ksw);
c4_lim = simplify(temp, 'Steps', 100);
disp('4) Degenerate limit f(t) for k2=ksw~=3k1')
disp(c4_lim)
latex_lim = latex(c4_lim);
disp('5) LaTeX for the degenerate limit:')
disp(latex_lim)
fprintf('\n')

% Numeric evaluation of the degenerate limit
c4_lim_vals = double(subs(c4_lim, {t,k1,ksw,D0}, {t_vals,1,1,D0_val}));

% Overlay on the same plot
plot(t_vals, c4_lim_vals, 'g--', 'LineWidth',1.0);
legend('expm(A t)','inverse-Laplace','degenerate limit (k_2=k_{sw}=1\neq3k_1)','Location','southeast');
hold off;
