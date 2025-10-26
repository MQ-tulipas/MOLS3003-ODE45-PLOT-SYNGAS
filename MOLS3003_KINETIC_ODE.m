%% Chemical Kinetics of Syngas

%% Kinetic ODE system for acid-catalyzed syngas formation
% Units: mol per g (of catalyst) per s
% k is rate constant, defined by reaction order; in this case, we have a
% Zero- / First- reaction order (mol/g/s)

function dxdt = syngas_kinetics(t, x, k, n)
    % x = [R, H2, CO, CH4, CO2]
    CA  = x(1);
    dRdt = -k * CA^n;

    % Stoichiometric coefficients; normalised yields
    a = 582/869;  % H2
    b = 276/869;  % CO
    c = 4/869; % CH4
    d = 7/869; % CO2

    dH2dt  =  a * k * CA^n;
    dCOdt  =  b * k * CA^n;
    dCH4dt =  c * k * CA^n;
    dCO2dt =  d * k * CA^n;

    dxdt = [dRdt; dH2dt; dCOdt; dCH4dt; dCO2dt];
end


%% ODE Plot (Zero-order 1st)

k = 1e-4;    % trial rate constant [s^-1]
n = 0; % Zero-order

x0 = [1, 0, 0, 0, 0];   % initial concentrations
tspan = [0 1200];      % seconds; *tspan in ~20 mins

[t, x] = ode45(@(t,x) syngas_kinetics(t,x,k,n), tspan, x0);

figure;
plot(t, x(:,2), 'r', t, x(:,3), 'b', t, x(:,4), 'g', t, x(:,5), 'k');
ylim([-0.02 0.15]); % Sets y-limit to -0.02 to 0.15
ax = gca; % Sets magnitude of decimal points to 0.XX
ax.YAxis.TickLabelFormat = '%.2f';
xlabel('Time (s)');
ylabel('Concentration (mol/g)');
legend('$H_{2}$','$CO$','$CH_{4}$','$CO_{2}$', 'interpreter', 'latex');
title('Zero-Order Reaction Kinetics of Syngas');

%% First-Order

hold on;

clear n;

n = 1; % First-order

[t, x] = ode45(@(t,x) syngas_kinetics(t,x,k,n), tspan, x0);

figure;
plot(t, x(:,2), 'r', t, x(:,3), 'b', t, x(:,4), 'g', t, x(:,5), 'k');
ylim([-0.02 0.15]); % Sets y-limit to -0.02 to 0.15
ax = gca; % Sets magnitude of decimal points to 0.XX
ax.YAxis.TickLabelFormat = '%.2f';
xlabel('Time (s)');
ylabel('Concentration (mol/g)');
legend('$H_{2}$','$CO$','$CH_{4}$','$CO_{2}$', 'interpreter', 'latex');
title('First-Order Reaction Kinetics of Syngas');