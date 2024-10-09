% This Code goes over the exact harmonic wave in 2D of the Haiti 2010
% Earthquake, where the graph of the wave at a given time is plotted and
% using a 2nd order central finite difference scheme, an approximate
% plotting of the exact harmonic wave is obtained, and the error between the
% exact and approximate solution is plotted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-Dimension general Wave Equation: 
%               u_{tt} = c^{2} * (u_{xx} + u_{yy}) 

% Initial conditions:  
%               u(x,y,0) = sin(c1 * pi * x) * sin(c2 * pi * y), 0<x<1, 0<y<1

% Boundary conditions: 
%                   u(0,y,t) = 0 , t>0
%                   u(1,y,t) = 0 , t>0
%                   u(x,0,t) = 0 , t>0
%                   u(x,1,t) = 0 , t>0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear previous variables, graphs, figures etc.
clc; clear all; close all;

% *** Make sure xlsx doc is in the same folder ***
% Table read in
T = readtable('haiti_ground_station_data.xlsx');
% time was from 8 - 18 seconds, now 1-10 seconds
time_data = xlsread('haiti_ground_station_data.xlsx', 'G:G');
% x_data = xlsread('haiti_ground_station_data.xlsx', 'K:K');
% lcolumn = xlsread('haiti_ground_station_data.xlsx', 'L:L');
amplitude_data = xlsread('haiti_ground_station_data.xlsx', 'K:K');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Area of impact of Haiti Earthquake - 120 km^2
% Radius = 6,180.387232371 m 
% 35-60 seconds ~ 40 seconds
% Estimated Velocity = 154.509680809 m/s
% Harmonic Wave Equation = a*cos(omega*(t-x/c) - phi)
% u_exact = -v^2*amplitude_data.*(sin(omega*(time_data - x_data./v) - phi) + cos(omega*(time_data - x_data./v) - phi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for u_exact
f = 1; %Hz
v = 154.509680809; % m/s
omega = 2*pi*f;
phi = 0;

% Parameters for u_approximate
c = 3;  
C_cfl = 1/sqrt(2); % Courant-Friedrich Levy Condition Value
dx = 0.01;
dy = dx;
h = 1/sqrt(2); 
dt = h*(dx / c);

% Parameters for mesh
t = 0:dt:0.9; 
x = 0:dx:1; 
y = 0:dy:1; 

u_app = zeros(length(x),length(y),length(t));

% Extracting exact amplitude for plotting of the exact wave
amp = zeros(length(t),1);

% Based on our changed dt, grabbing our data points to the correct size
% matching length(t)
for f=1:length(t)
    amp(f) = amplitude_data(floor(57601/length(t))*f);
end

% u_exact from the Haiti 2010 Data sheet
u_ex = zeros(length(x),length(y),length(t));

% Setting Initial Conditions
for i = 2 : length(x) - 1 
    for j = 2 : length(y) - 1
        u_app(:, :, 1) = -v^2*amp(1).*(sin(omega*(1 - j./v) - phi) + cos(omega*(1 - i./v) - phi));
        u_ex(:, :, 1) = -v^2*amp(1).*(sin(omega*(1 - j./v) - phi) + cos(omega*(1 - i./v) - phi));
    end
end

% Setting Boundary Conditions
for i = 2 : length(x) - 1 
    for j = 2 : length(y) - 1
        u_ex(i, j, 2) = -v^2*amp(2).*(sin(omega*(2 - j./v) - phi) + cos(omega*(2 - i./v) - phi));
        u_app(i, j, 2) = -v^2*amp(2).*(sin(omega*(2 - j./v) - phi) + cos(omega*(2 - i./v) - phi));
    end
end

% Approximating using 2nd order central finite difference scheme
for k = 2 : length(t) - 1
    for i = 2 : length(x) - 1
        for j = 2 : length(y) - 1
            u_ex(i, j, k+1) = -v^2*amp(k+1).*(sin(omega*(k+1 - j./v) - phi) + cos(omega*(k+1 - i./v) - phi));
            u_app(i, j, k+1)= (h^2)*(u_app(i+1,j,k)-2*u_app(i,j,k)+u_app(i-1,j,k))...
                +(C_cfl^2)*(u_app(i,j+1,k)-2*u_app(i,j,k)+u_app(i,j-1,k)) + 2*u_app(i,j,k) - u_app(i,j,k-1);
        end
    end
end

% Obtaining the Error
e_sol = abs((u_app/60) - u_ex);
max_e_sol = max(max(max(e_sol)));

% Creating The mesh-grid for the plot
[x_m,y_m] = meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the wave plots
for j = 1 : length(t)
       % Exact wave plot
       subplot(3, 1, 1)
       view(45,45)
       colormap(cool);
       hold on
       function1_ex = surf(x_m, y_m, u_ex(:, :, j));
       hold off
       title(sprintf('Haiti Earthquake bouy data at t = %1.2f',t(j)));
       axis ([0 0.4 0 1 -200 250]);
       xlabel('x'); 
       ylabel('y');
       zlabel(sprintf('u_exact(x,y,t = %1.2f)',t(j)));
       
       % Numerical method approximated wave plot
       subplot(3, 1, 2)
       view(45,45)
       colormap(cool);
       hold on
       function1_app = surf(x_m, y_m, u_app(:, :, j)./60);
       hold off
       title(sprintf('Finite Elements Approximation at t = %1.2f',t(j)));
       axis ([0 0.4 0 1 -200 250]);
       xlabel('x'); 
       ylabel('y');
       zlabel(sprintf('u(x,y,t = %1.2f)',t(j)));
       
       % Error between numerical solution and exact solution wave plot
       subplot(3, 1, 3)
       function_error = surf(x_m, y_m, e_sol(:, :, j)); 
       title(sprintf('Error -> |Approx - Exact| at t = %1.2f',t(j)));
       axis ([0 1 0 1 0 max_e_sol]);
       xlabel('x'); 
       ylabel('y');
       zlabel(sprintf('e_sol(x,y,t = %1.2f)',t(j)));
       pause(0.06);
       
       if(j ~= length(t))
           delete(function1_app);
           delete(function1_ex);
           delete(function_error);
       end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Screen captures of the wave plots

t_set = 382*[0.08, 0.13, 0.86, 0.88];
for n = 1 : 4
    figure(n)
    subplot(2, 1, 1)
    view(45, 45)
    colormap(cool);
    hold on
    function1_ex = surf(x_m, y_m, u_ex(:, :, floor(t_set(n))));
    hold off
    title(sprintf('Haiti Earthquake bouy data at t = %1.2f',t(ceil(t_set(n)))));
    axis ([0 0.4 0 1 -200 250]);
    xlabel('x'); 
    ylabel('y');
    zlabel(sprintf('u_exact(x,y,t = %1.2f)',t(floor(t_set(n)))));

    subplot(2, 1, 2)
    view(45, 45)
    colormap(cool);
    hold on
    function1_app = surf(x_m, y_m, u_app(:, :, floor(t_set(n)))/60);
    hold off
    title(sprintf('Finite Elements Approximation at t = %1.2f',t(ceil(t_set(n)))));
    axis ([0 0.4 0 1 -200 250]);
    xlabel('x'); 
    ylabel('y');
    zlabel(sprintf('u(x,y,t = %1.2f)',t(floor(t_set(n)))));
end

% Screen captures of the error plots
t_set2 = [0.08, 0.13, 0.86, 0.88];
for n = 1 : 4
    figure(4 + n)
    function_error = surf(x_m, y_m, e_sol(:, :, floor(382*t_set2(n)))); 
    title(sprintf('Error -> |Approx - Exact| at t = %1.2f',t(floor(382*t_set2(n)))));
    axis ([0 0.4 0 1 0 max_e_sol]);
    xlabel('x'); 
    ylabel('y');
    zlabel(sprintf('e_sol(x,y,t = %1.2f)', t(floor(382*t_set2(n)))));
    pause(0.06);
    
    if(j ~= length(t))
       delete(function_error);
    end
end

y_set2 = [0.33, 0.33, 0.33, 0.39];
for n = 1 : 4
    figure(8 + n)
    error_array = e_sol(:, floor(101*y_set2(n)), floor(382*t_set2(n)));
    x = linspace(0, 1, 101);
    plot(x, error_array);
    title(sprintf('Error in y compared to x t = %1.2f', t(floor(382*t_set2(n)))));
    xlabel('x'); 
    ylabel('y');

end