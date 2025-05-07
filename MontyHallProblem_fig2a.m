% MontyHallProblem042525_fig2a.m
% Integrate the full Monty Hall network (always–switch) 
% with k1=k2=ksw=1, [P]=1, [D_j]=100, and overlay the symmetric‐cascade analytical curve.
%
% Oliver Steinbock (FSU, 2025)

%% 0) Parameters
k1  = 1;
k2  = 1;
ksw = 1;
kst = 0;

P0  = 1;
D0  = 100;          % initial [D1]=[D2]=[D3]
tspan = [0 0.1];

%% 1) Integrate full network
y0 = [ P0; D0; D0; D0; ...    % P, D1–D3
       zeros(3,1); ...        % I1–I3
       zeros(4,1); ...        % I12, I13, I23, I32
       zeros(4,1); ...        % I231, I321, I123, I132
       zeros(4,1) ];          % S121, S131, S232, S323

[t, y] = ode45(@(t,y) odefun_full(t,y,k1,k2,ksw,kst), tspan, y0);

%% 2) Plot results with analytical overlay
plotResults(t, y, k1, k2, ksw, kst, D0);

%% ODE function
function dydt = odefun_full(~, y, k1, k2, ksw, kst)
    % Unpack
    P    = y(1);
    D1   = y(2); D2 = y(3); D3 = y(4);
    I1   = y(5); I2 = y(6); I3 = y(7);
    I12  = y(8); I13 = y(9); I23 = y(10); I32 = y(11);
    I231 = y(12); I321 = y(13); I123 = y(14); I132 = y(15);
    S121 = y(16); S131 = y(17); S232 = y(18); S323 = y(19);

    % dP
    dP = -k1 * P * (D1 + D2 + D3);

    % dD 
    dD1 = -k1*P*D1                 - ksw*(I23+I32)*D1; % door-1 is not revealed
    dD2 = -k1*P*D2 - k2*(I1+I3)*D2 - ksw*I13*D2;
    dD3 = -k1*P*D3 - k2*(I1+I2)*D3 - ksw*I12*D3;

    % dI
    dI1 =  k1*P*D1 - k2*I1*(D2+D3);
    dI2 =  k1*P*D2 - k2*I2*D3;
    dI3 =  k1*P*D3 - k2*I3*D2;

    % dI_{j,k}
    dI12 = k2*I1*D2 - ksw*I12*D3;
    dI13 = k2*I1*D3 - ksw*I13*D2;
    dI23 = k2*I2*D3 - ksw*I23*D1;
    dI32 = k2*I3*D2 - ksw*I32*D1;

    % dI_{j,k,ℓ} (switch products)
    dI231 = ksw*I23*D1;
    dI321 = ksw*I32*D1;
    dI123 = ksw*I12*D3;
    dI132 = ksw*I13*D2;

    % Pack derivatives
    dydt = [ dP;
             dD1; dD2; dD3;
             dI1; dI2; dI3;
             dI12; dI13; dI23; dI32;
             dI231; dI321; dI123; dI132;
             zeros(4,1) ];  % no stay branch when kst=0
end

%% Plot function
function plotResults(t, y, k1, k2, ksw, kst, D0)
    % Aggregated traces
    swWins        = y(:,12) + y(:,13);       % I231 + I321
    intermediates1 = sum(y(:,5:7),2);        % I1 + I2 + I3
    intermediates2 = sum(y(:,8:11),2);       % I12 + I13 + I23 + I32
    losses        = y(:,14) + y(:,15) + y(:,18) + y(:,19);

    figure(1); clf; set(gcf,'Color','w'); hold on;

    % Simulation curves
    plot(t, y(:,1),         '-k', 'LineWidth',1.5);             % P
    plot(t, swWins,         'Color',[0 0.9 0], 'LineWidth',2.5);% switch wins
    plot(t, intermediates1, ':b', 'LineWidth',1.5);            % stage-1 intermediates
    plot(t, intermediates2, ':m', 'LineWidth',1.5);            % stage-2 intermediates
    plot(t, losses,         '--r','LineWidth',1.5);             % losses

    % Analytical overlay (degenerate‐limit expression)
    tvec = t;
    c4_sym = ( ksw.^2 .* exp(-3*D0*k1.*tvec) .* exp(-D0*ksw.*tvec) ...
              .* ( exp(3*D0*k1.*tvec) - exp(D0*ksw.*tvec) ) ) ...
            ./ ( 3*(3*k1 - ksw).^2 ) ...
          - ( D0*ksw.^2 .* tvec .* exp(-D0*ksw.*tvec) ) ...
            ./ ( 3*(3*k1 - ksw) ) ...
          - ( exp(-D0*ksw.*tvec) ...
              .* ( D0*ksw.*tvec - exp(D0*ksw.*tvec) + 1 ) ) / 3;
    
    W_exact = 2 * c4_sym;
    plot(tvec, W_exact, '-.k', 'LineWidth', 1.0);


    % Reference lines
    yline(1/3, '--b', '1/3', ...
        'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','LineWidth',1);
    yline(2/3, '--b', '2/3', ...
        'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','LineWidth',1);
    ylim([-0.02,1.02]);

    % Sub-panel label
    text(-0.15, 1.05, '\bf a', ...
        'Units','normalized','FontSize',22,'FontWeight','bold','VerticalAlignment','top');

    hold off;
    ax = gca; ax.FontSize = 14; ax.LineWidth = 1.5; ax.Box = 'on';
    xlabel('Time','FontSize',14);
    ylabel('Concentration','FontSize',14);
    legend('Player species P','Switch wins','Stage-1 intermediates ','Stage-2 intermediates', ...
        'Switch losses','Analytical result','Location','northeast', ...
        'FontSize',10,'NumColumns',2 );
end
