%% Zach Freyman - AEE342 Project 2
% Source Panel Method for NACA 0012 Airfoil

%% Initialization of key symbols and variables
close all;
clear


%% Generating Airfoil Plot
% Equations sourced from https://turbmodels.larc.nasa.gov/naca0012_val.html

% Analytical Solution of Airfoil Shape (Independent of panel count)
n = 250;
x1 = linspace(0, 1, n/2);
y_upper = 0.594689181.*(0.298222773 .* sqrt(x1) - 0.127125232.* x1 - 0.357907906 .* x1.^2 + 0.291984971 .* x1.^3 - 0.105174606 .* x1.^4);
y_lower = -0.594689181.*(0.298222773 .* sqrt(x1) - 0.127125232.* x1 - 0.357907906 .* x1.^2 + 0.291984971 .* x1.^3 - 0.105174606 .* x1.^4);

% Setting up panels
p = 10; % Arbitrary number of panels
alpha = 0; %Adding the option for an angle of attack

x = 0:2/p:1; % This array will be used again later for panel location hence the setup using panel number
f = 0.594689181.*(0.298222773 .* sqrt(x) - 0.127125232.* x - 0.357907906 .* x.^2 + 0.291984971 .* x.^3 - 0.105174606 .* x.^4);
Y = zeros(p, 1);
X = linspace(0, 1, p/2);

for i = 1:p + 1
    if i<= 0.5 * p + 1
        Y(i) = f(i);
        X(i) = x(i);
    else
        Y(i) = -1 .* f(p + 2 - i);
        X(i) = x(p + 2 - i);
    end
end

figure(1)
hold on
plot(x1, y_upper)
plot(x1, y_lower)
plot(X, Y)
plot(X, Y, 'bs', 'MarkerSize', 4)
grid on
axis equal
hold off
title('NACA 0012 Airfoil')
legend('Airfoil Surface','Panel Endpoints')

%% Setting up control points and establishing panel lengths
dx = zeros(p, 1);
dy = zeros(p, 1);
Xc = zeros(p, 1);
Yc = zeros(p, 1);
Length = zeros(p, 1);
phi = zeros(p, 1);
beta = zeros(p, 1);

for i = 1:p
    dx(i) = X(i+1) - X(i);
    dy(i) = Y(i+1) - Y(i);
    Xc(i) = X(i) + 0.5 * dx(i);
    Yc(i) = Y(i) + 0.5 * dy(i);
    Length(i) = sqrt(dx(i).^2 + dy(i).^2);
    phi(i) = atan2(dy(i), dx(i));
    beta(i) = phi(i) + pi/2 - alpha;
end

figure(2)
hold on
plot(x1, y_upper)
plot(x1, y_lower)
plot(X, Y)
plot(Xc, Yc, 'bs', 'MarkerSize', 4)
grid on
axis equal
hold off
title('NACA 0012 Airfoil')
legend('Airfoil Surface','Panel Control Points')

%% Creating Influence Coefficient Matrices (I and J)
% Setting up coefficents and Matrices
I = zeros(p);
J = zeros(p);
A = zeros(p);
B = zeros(p);
C = zeros(p);
D = zeros(p);
E = zeros(p);

for i = 1:p
    for j = 1:p
        A(i,j) = -1 * (Xc(i) - X(j)) * cos(phi(j)) - (Yc(i) - Y(j)) * sin(phi(j));
        B(i,j) = (Xc(i) - X(j))^2 + (Yc(i) - Y(j))^2;
        C(i,j) = sin(phi(i) - phi(j));
        D(i,j) = (Yc(i) - Y(j)) * cos(phi(i)) - (Xc(i) - X(j)) * sin(phi(i));
        E(i,j) = (Xc(i) - X(j)) * sin(phi(j)) - (Yc(i) - Y(j)) * cos(phi(j));
        %LengthS(i,j) = ;
        if i == j
            I(i,j) = 0.5;
        else
            I(i,j) = %Something
        end
    end
end



%imagesc(I)