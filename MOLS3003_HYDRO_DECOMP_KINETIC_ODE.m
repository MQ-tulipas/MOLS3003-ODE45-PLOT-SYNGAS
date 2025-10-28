%% Chemical Kinetics of Syngas Production Model
% PET Hydrolysis + Acid Decomposition

%IM -> Intermediate

function dxdt = pet_syngas(t, x, k1, k2, n1, n2)
    % x = [PET, IM, H2, CO, CH4, CO2]
    CP = x(1); CI = x(2); 

    % Reaction 1: PET hydrolysis - base catalysed
    R1 = k1 * CP^n1;

    % Reaction 2: Acid decomposition - syngas formation
    R2 = k2 * CI^n2;

    % Stoichiometric coefficients - normalised to total syngas yield
    total = 582 + 276 + 4 + 7;
    a = 582 / total;  % H2
    b = 276 / total;  % CO
    c = 4 / total;    % CH4
    d = 7 / total;    % CO2

    % Differential equations
    dRPdt   = -R1;
    dRIdt   =  R1 - R2;
    dH2dt   =  a * R2;
    dCOdt   =  b * R2;
    dCH4dt  =  c * R2;
    dCO2dt  =  d * R2;

    dxdt = [dRPdt; dRIdt; dH2dt; dCOdt; dCH4dt; dCO2dt];
end

%% ODE Solution

k1 = 1e-3;   % PET hydrolysis rate constant
k2 = 5e-4;   % Syngas formation rate constant
n1 = 1;      % First-order hydrolysis
n2 = 0;      % Zero-order decomposition

x0 = [1, 0, 0, 0, 0, 0];   % Initial PET = 1 mol/g, others = 0
tspan = [0 1200];           % seconds

[t, x] = ode45(@(t,x) pet_syngas(t,x,k1,k2,n1,n2), tspan, x0);

%% Plot of Model

figure;
plot(t, x(:,3), 'r', t, x(:,4), 'b', t, x(:,5), 'g', t, x(:,6), 'k', 'LineWidth', 1.5);
ylim([-0.05 0.45]); % Sets y-limit to -0.05 to 0.15
ax = gca; % Sets magnitude of decimal points to 0.XX
ax.YAxis.TickLabelFormat = '%.2f';
xlabel('Time (s)');
ylabel('Concentration (mol/g)');
legend('$H_{2}$','$CO$','$CH_{4}$','$CO_{2}$', 'interpreter', 'latex');
title('Kinetic Model of PET → Na_2TPA + EG → Syngas');
grid on;
