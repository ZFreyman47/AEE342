%% Zach Freyman - AEE342 Project 2
% Source Panel Method for NACA 4 Digit Airfoils

%% Initialization and Refresh
close all;
clear


%% Generating Airfoil Plot

% Primary Simulation Settings
default = {'100', '0', '1', '12'}; % 
userinput = inputdlg({'Panel Count', 'Angle of Attack', 'Freestream Velocity', 'NACA Thickness (Last two digits)'}, 'Program Inputs', 1, default);
p = str2double(userinput{1});
alpha = str2double(userinput{2})*(2*pi)/360;
V_inf = str2double(userinput{3});
thickness = (str2double(userinput{4}))/100;

x = 0:2/p:1; % This array will be used again later for panel location hence the setup using panel number
f = 5 * thickness * (0.2969 .* sqrt(x) - 0.126 .* x - 0.3516 .* x.^2 + 0.2843 .* x.^3 - 0.1015 .* x.^4);
Y = zeros(p, 1);
X = linspace(0, 1, p/2);

% Analytical Solution of Airfoil Shape (Independent of panel count)
n = 500; % Number used as sufficent to generate accurate and visually appealing airfoil shape
x1 = linspace(0, 1, n/2);
y_upper = 5 * thickness * (0.2969 .* sqrt(x1) - 0.126 .* x1 - 0.3516 .* x1.^2 + 0.2843 .* x1.^3 - 0.1015 .* x1.^4);
y_lower = -5 * thickness * (0.2969 .* sqrt(x1) - 0.126 .* x1 - 0.3516 .* x1.^2 + 0.2843 .* x1.^3 - 0.1015 .* x1.^4);

for i = 1:p + 1
    if i<= 0.5 * p + 1
        Y(i) = f(i);
        X(i) = x(i);
    else
        Y(i) = -1 .* f(p + 2 - i);
        X(i) = x(p + 2 - i);
    end
end

% Plotting Airfoil with Panel Endpoints
figure(1)
hold on
plot(x1, y_upper)
plot(x1, y_lower)
plot(X, Y)
plot(X, Y, 'bs', 'MarkerSize', 4)
grid on
axis equal
hold off
title('NACA 0012 Airfoil and Panel Endpoints')
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

% Plotting Airfoil with Panel Control Points
figure(2)
hold on
plot(x1, y_upper)
plot(x1, y_lower)
plot(X, Y)
plot(Xc, Yc, 'bs', 'MarkerSize', 4)
grid on
axis equal
hold off
title('NACA 0012 Airfoil and Panel Control Points')
legend('Airfoil Surface','Panel Control Points')

%% Creating Influence Coefficient Matrix (I)
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
        if i == j
            I(i,j) = ((1)/(2*pi))*pi;
        else
            I(i,j) = ((1)/(2*pi))*((C(i,j)/2)*log((Length(j)^2 + 2 * A(i,j) * Length(j) + B(i,j))/(B(i,j))) + ((D(i,j) - A(i,j) * C(i,j))/(E(i,j))) * (atan((Length(j) + A(i,j))/(E(i,j))) - atan(A(i,j)/E(i,j))));
        end
    end
end

%% Solving for Source Strengths
n_velocity = -1 * V_inf * cos(beta);
t_velocity = V_inf * sin(beta);

lambda = mldivide(I, n_velocity);

%% Creating Influence Coefficient Matrix (J)

for i = 1:p
    for j = 1:p
        if i == j
            J(i,j) = 0;
        else
            J(i,j) = (D(i,j) - A(i,j) * C(i,j))/(2 * E(i,j)) * log((Length(j)^2 + 2 * A(i,j) * Length(j) + B(i,j)) / (B(i,j))) - C(i,j) * (atan((Length(j) + A(i,j)) / (E(i,j))) - atan((A(i,j))/(E(i,j))));
        end

    end
end

V_tan = J * lambda(:) * (1/(2*pi)) + t_velocity;

Cp = zeros(1,p);
for i = 1:p
    Cp(i) = 1 - ((V_tan(i) / V_inf).^2);
end

mass_flux = 0;
for i = 1:p
    mass_flux = mass_flux + lambda(i) .* Length(i); 
end

disp(mass_flux)

%% General Plotting

% Sanity Checks for Matrices
figure(3)
imagesc(I)
title('Sanity Check for I Matrix')
figure(4)
imagesc(J)
title('Sanity Check for J Matrix')

% Chord vs. Source Strength
figure(5)
subplot(2,1,1)
plot(Xc(1:p/2), lambda(1:p/2))
title('Chord vs. Source Strength (Upper Surface)')
xlabel('X/c')
ylabel('Lambda')
subplot(2,1,2)
plot(Xc(p/2 +1:p), lambda(p/2 +1:p))
title('Chord vs. Source Strength (Lower Surface)')
xlabel('X/c')
ylabel('Lambda')


% Chord vs. Coefficient of Pressure
Cpdata_x_0 = importdata('0_x.txt');
Cpdata_c_0 = importdata('0_c.txt');
Cpdata_x_10 = importdata('10_x.txt');
Cpdata_c_10 = importdata('10_c.txt');


figure(6)
hold on
plot(Xc, Cp)
if alpha == 0
    plot(Cpdata_x_0, Cpdata_c_0, 'o')
elseif alpha == 10*pi/180
    plot(Cpdata_x_10, Cpdata_c_10, 'o')
else
end

set(gca, "YDir", "reverse")
title("Chord vs. Coefficient of Pressure at angle of attack of " + alpha*(360/(2*pi)) + " degrees ")
xlabel('X/c')
ylabel('Coefficient of Pressure')
hold off

%% Streamline Plot

[x_grid, y_grid] = meshgrid(-0.1:0.01:1.1, -0.4:0.001:0.4);


u = V_inf * cos(alpha);
v = V_inf * sin(alpha);

for i = 1:p
    u = u + ((lambda(i))./(2 * pi)) .* ((x_grid - Xc(i))./((x_grid - Xc(i)).^2 + (y_grid - Yc(i)).^2)) .* Length(i);
    v = v + ((lambda(i))./(2 * pi)) .* ((y_grid - Yc(i))./((x_grid - Xc(i)).^2 + (y_grid - Yc(i)).^2)) .* Length(i);
end

figure(7)
hold on
streamslice(x_grid, y_grid, u, v, 7)
plot(x1, y_upper, 'k', 'LineWidth', 2)
plot(x1, y_lower, 'k', 'LineWidth', 2)
title("Flow Streamlines over NACA 00" + thickness*100 + " at angle of attack of " + alpha*(360/(2*pi)) + " degrees ")
xlim([-0.1, 1.1])
ylim([-0.4, 0.4])
hold off


massimbalance_plot = importdata('mass_imbalance_20.50.100.200.txt');
panel_count = [20, 50, 100, 200];
figure(8)
plot(panel_count, massimbalance_plot)
title('Mass Imbalance Versus Panel Count')
xlabel('Panel Count')
ylabel('Simulation Mass Imbalance')