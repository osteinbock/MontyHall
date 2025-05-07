% montyHall_sweepD_compare.m
% Sweep initial [D_j] from 0.1 to 5 and compare final win concentrations
% for Always‐Stay vs Always‐Switch strategies.

% time span for integration
tspan = [0 2];

% range of initial door concentrations
Dj_vals = linspace(0.0,6,30);

% preallocate
wins_stay   = zeros(size(Dj_vals));
wins_switch = zeros(size(Dj_vals));

% common parameters
k1 = 1;
k2 = 1;

for ii = 1:numel(Dj_vals)
    Dj = Dj_vals(ii);
    % initial conditions
    P0   = 1;
    D0   = Dj*[1;1;1];
    I0   = zeros(3,1);
    I2_0 = zeros(4,1);
    I3_0 = zeros(4,1);
    S0   = zeros(4,1);
    y0   = [P0; D0; I0; I2_0; I3_0; S0];

    % Always‐Stay: ksw=0, kst=1
    [~,Y] = ode45(@(t,y) odefun_full(t,y,k1,k2,0,1), tspan, y0);
    yEnd = Y(end,:);
    % total wins (ℓ=1 products)
    wins_stay(ii) = yEnd(12) + yEnd(13) + yEnd(16) + yEnd(17);

    % Always‐Switch: ksw=1, kst=0
    [~,Y] = ode45(@(t,y) odefun_full(t,y,k1,k2,1,0), tspan, y0);
    yEnd = Y(end,:);
    % total wins (ℓ=1 products)
    wins_switch(ii) = yEnd(12) + yEnd(13);
end

% Plot comparison
figure(1);
set(gcf,'Color','w');  % white background
plot(Dj_vals, wins_stay,   '-ob','LineWidth',1.5,'MarkerSize',6); hold on
plot(Dj_vals, wins_switch, '-sk','LineWidth',1.5,'MarkerSize',6); hold off

% horizontal reference lines, labels on right
yline(1/3, '--b', '1/3',  ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom', ...
    'LineWidth',1.0);
yline(2/3, '--b', '2/3',  ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom', ...
    'LineWidth',1.0);

xlabel('Initial door concentration [D_j]_0','FontSize',14);
ylabel('Final total win concentration','FontSize',14);
set(gca,'FontSize',12,'Box','on');
xlim([min(Dj_vals) max(Dj_vals)]);
ylim([-0.02 2/3+0.02]);

legend('Always–stay','Always–switch','Location','SouthEast','FontSize',10);

% axes formatting
ax = gca;
ax.FontSize   = 14;     % increase font size
ax.LineWidth  = 1.5;    % thicker frame
ax.Box        = 'on';   % keep box
    
% sub-panel label (a or b; adjust manually)
text(-0.15, 1.05, '\bf b', ...
 'Units','normalized', 'FontSize', 22, ...
 'FontWeight','bold', 'VerticalAlignment','top');
    


function dydt = odefun_full(~, y, k1, k2, ksw, kst)
    % unpack
    P    = y(1);
    D1= y(2); D2= y(3); D3= y(4);
    I1= y(5); I2= y(6); I3= y(7);
    I12= y(8); I13= y(9); I23= y(10); I32= y(11);
    I231= y(12); I321= y(13); I123= y(14); I132= y(15);
    S121= y(16); S131= y(17); S232= y(18); S323= y(19);

    % dP
    dP = -k1*(P*D1 + P*D2 + P*D3);

    % dD1
    dD1 = -k1*P*D1 ...
          - ksw*(I23*D1 + I32*D1) ...
          - kst*(I12*D1 + I13*D1);

    % dD2
    dD2 = -k1*P*D2 ...
          - k2*(I1*D2 + I3*D2) ...
          - ksw*(I13*D2) ...
          - kst*(I23*D2);

    % dD3
    dD3 = -k1*P*D3 ...
          - k2*(I1*D3 + I2*D3) ...
          - ksw*(I12*D3) ...
          - kst*(I32*D3);

    % dI1, dI2, dI3
    dI1 = k1*P*D1 - k2*(I1*D2 + I1*D3);
    dI2 = k1*P*D2 - k2*(I2*D3);
    dI3 = k1*P*D3 - k2*(I3*D2);

    % dI_{j,k}
    dI12 = k2*I1*D2 - ksw*(I12*D3) - kst*(I12*D1);
    dI13 = k2*I1*D3 - ksw*(I13*D2) - kst*(I13*D1);
    dI23 = k2*I2*D3 - ksw*(I23*D1) - kst*(I23*D2);
    dI32 = k2*I3*D2 - ksw*(I32*D1) - kst*(I32*D3);

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