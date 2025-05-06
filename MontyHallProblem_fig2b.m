% montyHall_alwaysStay.m
% Integrate the Monty Hall network with always-stay (kst=1, ksw=0), k1=k2=1,
% [P]_0=1, [D_j]_0=100

% Rate constants
k1  = 1;
k2  = 1;
ksw = 0;   % no switching
kst = 1;   % always stay

% Initial conditions
P0   = 1;
D0   = 100;                  % same for D1,D2,D3
y0   = [ P0; repmat(D0,3,1);    ...  % P, D1–D3
         zeros(3,1);             ...  % I1–I3
         zeros(4,1);             ...  % I12,I13,I23,I32
         zeros(4,1);             ...  % I231,I321,I123,I132
         zeros(4,1) ];              % S121,S131,S232,S323

tspan = [0 0.1];
[t,y]  = ode45(@(t,y) odefun_full(t,y,k1,k2,ksw,kst), tspan, y0);

% Plot (pass parameters into the plotting function)
plotResults(t, y, k1, k2, D0, kst);


function dydt = odefun_full(~, y, k1, k2, ksw, kst)
    % unpack
    P    = y(1);
    D1=y(2); D2=y(3); D3=y(4);
    I1=y(5); I2=y(6); I3=y(7);
    I12=y(8); I13=y(9); I23=y(10); I32=y(11);
    I231=y(12); I321=y(13); I123=y(14); I132=y(15);
    S121=y(16); S131=y(17); S232=y(18); S323=y(19);

    % dP
    dP = -k1*P*(D1 + D2 + D3);

    % dDj (no door-1 reveal in stage-2)
    dD1 = -k1*P*D1                   - ksw*(I23 + I32)*D1 - kst*(I12 + I13)*D1;
    dD2 = -k1*P*D2 - k2*(I1 + I3)*D2 - ksw*I13*D2         - kst*I23*D2;
    dD3 = -k1*P*D3 - k2*(I1 + I2)*D3 - ksw*I12*D3         - kst*I32*D3;

    % dI_j
    dI1 = k1*P*D1 - k2*I1*(D2 + D3);
    dI2 = k1*P*D2 - k2*I2*D3;
    dI3 = k1*P*D3 - k2*I3*D2;

    % dI_{j,k}
    dI12 = k2*I1*D2 - ksw*I12*D3 - kst*I12*D1;
    dI13 = k2*I1*D3 - ksw*I13*D2 - kst*I13*D1;
    dI23 = k2*I2*D3 - ksw*I23*D1 - kst*I23*D2;
    dI32 = k2*I3*D2 - ksw*I32*D1 - kst*I32*D3;

    % dI_{j,k,ℓ} (switch products) – zero because ksw=0
    dI231 = ksw*I23*D1;
    dI321 = ksw*I32*D1;
    dI123 = ksw*I12*D3;
    dI132 = ksw*I13*D2;

    % dS_{j,k,j} (stay products)
    dS121 = kst*I12*D1;
    dS131 = kst*I13*D1;
    dS232 = kst*I23*D2;
    dS323 = kst*I32*D3;

    % pack all 19 derivatives
    dydt = [ dP;
             dD1; dD2; dD3;
             dI1; dI2; dI3;
             dI12; dI13; dI23; dI32;
             dI231; dI321; dI123; dI132;
             dS121; dS131; dS232; dS323 ];
end


function plotResults(t, y, k1, k2, D0, kst)
    % aggregated signals
    stayWins       = y(:,16) + y(:,17);       % I121 + I131
    intermediates1 = sum(y(:,5:7),2);        % I1 + I2 + I3
    intermediates2 = sum(y(:,8:11),2);       % I12 + I13 + I23 + I32
    totalLosses    = y(:,14) + y(:,15) + y(:,18) + y(:,19);

    figure(1); clf; set(gcf,'Color','w'); hold on;
    plot(t, y(:,1),          '-k', 'LineWidth',1.5); % P
    plot(t, stayWins,        '-',  'Color',[0 .7 0], 'LineWidth',2.5);
    plot(t, intermediates1,  ':b', 'LineWidth',1.5);
    plot(t, intermediates2,  ':m', 'LineWidth',1.5);
    plot(t, totalLosses,     '--r','LineWidth',1.5);

    % Analytical overlay for always-stay
    tvec = t;
    I121 =  1/6 ...
        - (k1*k2*kst) .* exp(-3*D0*k1 .* tvec) ...
          ./ ((3*k1 - 2*k2) .* (3*k1 - kst) * 3*k1) ...
        + (k1*k2*kst) .* exp(-2*D0*k2 .* tvec) ...
          ./ ((3*k1 - 2*k2) .* (2*k2 - kst) * 2*k2) ...
        - (k1*k2*kst) .* exp(-D0*kst .* tvec) ...
          ./ ((3*k1 - kst) .* (2*k2 - kst) * kst);
    W_exact = 2 * I121;
    plot(tvec, W_exact, '-.k', 'LineWidth',1.0);

    % reference lines at 1/3 and 2/3
    yline(1/3, '--b', '1/3', ...
          'LabelHorizontalAlignment','right', ...
          'LabelVerticalAlignment','bottom', 'LineWidth',1);
    yline(2/3, '--b', '2/3', ...
          'LabelHorizontalAlignment','right', ...
          'LabelVerticalAlignment','bottom', 'LineWidth',1);
    ylim([-0.02,1.02]);

    text(-0.15,1.05,'\bf b','Units','normalized', ...
         'FontSize',22,'FontWeight','bold','VerticalAlignment','top');

    hold off;
    ax = gca; ax.FontSize=14; ax.LineWidth=1.5; ax.Box='on';
    xlabel('Time','FontSize',14);
    ylabel('Concentration','FontSize',14);
    legend('Player species P','Stay wins','Stage‑1 Intermediates','Stage‑2 Intermediates','Stay losses', ...
           'Analytical result','Location','northeast','FontSize',10,'NumColumns',2);
end
