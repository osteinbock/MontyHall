% montyHall_sweepKst.m
% Sweep kst from 1e-3 to 1e3 for ksw=1, D0=100, and plot final win fraction.

% Parameters
k1    = 1;
k2    = 1;
ksw   = 1;            % switch‐branch rate constant
D0    = 100;          % initial door concentration
tspan = [0 10];       % integration interval (long enough to reach steady state)

% Sweep kst over log-scale
kst_vals = logspace(-3,3,50);
wins     = zeros(size(kst_vals));

% Loop over kst
for ii = 1:numel(kst_vals)
    kst = kst_vals(ii);

    % Initial conditions
    P0 = 1;
    y0 = [ P0;         % P
           D0; D0; D0; % D1,D2,D3
           zeros(3,1); % I1,I2,I3
           zeros(4,1); % I12,I13,I23,I32
           zeros(4,1); % I231,I321,I123,I132
           zeros(4,1)];% S121,S131,S232,S323

    % Integrate full network until tspan
    [~, Y] = ode45(@(t,y) odefun_full(t,y,k1,k2,ksw,kst), tspan, y0);

    yEnd = Y(end,:);
    % total wins = switch wins (I231+I321) + stay wins (I121+I131)
    wins(ii) = yEnd(12) + yEnd(13) + yEnd(16) + yEnd(17);
end

figure(1); clf; set(gcf,'Color','w');

% Steady‐state analytic overlay
W_inf = (2 + kst_vals) ./ (3*(1 + kst_vals));
semilogx(kst_vals, W_inf, '-k', 'LineWidth',1.5);
hold on

% Plot numerical results
semilogx(kst_vals, wins, 'o', 'Color',[0 0.7 0], 'LineWidth',1.5,'MarkerSize',5);

% Reference lines at 1/3 and 2/3
yline(1/3,'--b','1/3','LabelHorizontalAlignment','right','LineWidth',1.5);
yline(2/3,'--b','2/3','LabelHorizontalAlignment','right','LineWidth',1.5);

% sub-panel label (a or b; adjust manually)
text(-0.15, 1.05, '\bf a', ...
 'Units','normalized', 'FontSize', 22, ...
 'FontWeight','bold', 'VerticalAlignment','top');
    
hold off

% Formatting
legend('Analytical result','Numerical result','FontSize',10','Location','southwest');
ax = gca; ax.FontSize=14; ax.LineWidth=1.5; ax.Box='on';
xlabel('Stay-branch rate constant k_{st}','FontSize',14);
ylabel('Final total win concentration','FontSize',14);
xlim([min(kst_vals), max(kst_vals)]);
ylim([0.3, 0.7]);


% ------------------------------------------------------------------------------
function dydt = odefun_full(~, y, k1, k2, ksw, kst)
    % unpack
    P    = y(1);
    D1= y(2); D2= y(3); D3= y(4);
    I1= y(5); I2= y(6); I3= y(7);
    I12= y(8); I13= y(9); I23= y(10); I32= y(11);
    I231= y(12); I321= y(13); I123= y(14); I132= y(15);
    S121= y(16); S131= y(17); S232= y(18); S323= y(19);

    % dP
    dP = -k1*P*(D1 + D2 + D3);

    % dDj
    dD1 = -k1*P*D1 -               - ksw*(I23+I32)*D1 - kst*(I12+I13)*D1;
    dD2 = -k1*P*D2 - k2*(I1+I3)*D2 - ksw*(I13)*D2     - kst*(I23)*D2;
    dD3 = -k1*P*D3 - k2*(I1+I2)*D3 - ksw*(I12)*D3     - kst*(I32)*D3;

    % dI_j
    dI1 = k1*P*D1 - k2*I1*(D2+D3);
    dI2 = k1*P*D2 - k2*I2*D3;
    dI3 = k1*P*D3 - k2*I3*D2;

    % dI_{j,k}
    dI12 = k2*I1*D2 - ksw*I12*D3 - kst*I12*D1;
    dI13 = k2*I1*D3 - ksw*I13*D2 - kst*I13*D1;
    dI23 = k2*I2*D3 - ksw*I23*D1 - kst*I23*D2;
    dI32 = k2*I3*D2 - ksw*I32*D1 - kst*I32*D3;

    % dI_{j,k,ℓ} (switch products)
    dI231 = ksw*I23*D1;
    dI321 = ksw*I32*D1;
    dI123 = ksw*I12*D3;
    dI132 = ksw*I13*D2;

    % dS_{j,k,j} (stay products)
    dS121 = kst*I12*D1;
    dS131 = kst*I13*D1;
    dS232 = kst*I23*D2;
    dS323 = kst*I32*D3;

    % pack
    dydt = [dP;
            dD1; dD2; dD3;
            dI1; dI2; dI3;
            dI12; dI13; dI23; dI32;
            dI231; dI321; dI123; dI132;
            dS121; dS131; dS232; dS323];
end
