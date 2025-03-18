%% Zach Freyman - AEE342 Project 2
% Source Panel Method for NACA 0012 Airfoil

%% Initialization of key symbols and variables
close all;
clear


%% Generating Airfoil Plot
% Equations sourced from https://turbmodels.larc.nasa.gov/naca0012_val.html

p = 100; % Arbitrary number of panels

x = 0:(2/p):1; % This array will be used again later for panel location hence the setup using panel number
y_upper = 0.594689181.*(0.298222773 .* sqrt(x) - 0.127125232.* x - 0.357907906 .* x.^2 + 0.291984971 .* x.^3 - 0.105174606 .* x.^4);
y_lower = -0.594689181.*(0.298222773 .* sqrt(x) - 0.127125232.* x - 0.357907906 .* x.^2 + 0.291984971 .* x.^3 - 0.105174606 .* x.^4);

figure(1)
hold on
plot(x, y_upper)
plot(x, y_upper, 'bs', 'MarkerSize', 4)
plot(x, y_lower)
plot(x, y_lower, 'rs', 'MarkerSize', 4)
grid on
axis equal
hold off
title('NACA 0012 Airfoil')
legend('Upper Surface','Upper Panel Locations', 'Lower Surface', 'Lower Panel Locations')

%% Creating Influence Coefficient Matrices (I and J)

