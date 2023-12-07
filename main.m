%% Main - attempt 2
clear
clc

%% Define initial conditions and necessary variables
% Initial conditions
I = diag([100 100 50]);
eul = [0 pi/2 0]';
W0 = [0.2 0.2 1]';

% Place initial conditions in a single matrix
X0 = [W0; eul];

% Settling time and damping ratio desired for the controller
ts_desired = 60;
zeta = 0.75;

% Define timespan for integration - change multiplier as desired
T = ts_desired;
tspan = linspace(0,T,1000);

% Integration options
tol = 1e-13;
options = odeset('RelTol', tol, 'AbsTol', tol, 'Events', @eventsFunction);

%% Controller design
% Calculate the K values for the controllers
[K1, K2] = CalculateGain(I, ts_desired, zeta);

%% Drive axial velocity to zero
% Integrations for velocity controllers
[txy,xy] = ode45(@(t, xy) diffEq(t, xy, K1, I), tspan, X0, options);

% Find values at end of integration
final_value_xy = xy(end,:);
final_t_xy = txy(end);

% Define time span for axial velocity integration
tz_span = linspace(final_t_xy, 2.5*T, 1000);

%% Drive transverse velocity to zero
% Integration for axial velocity controller
[tz, z] = ode45(@(t, z) diffEq_z(t, z, K2, I), tz_span, final_value_xy, options);

% Find values at end of integration
final_value_z = z(end,:);
final_t_z = tz(end);

% Define time span for phi integration
tphi_span = linspace(final_t_z, 5*T, 1000);

%% Drive phi to zero
% Integration for phi controller
[tphi,phi] = ode45(@(t,phi) diffEq_phi(t, phi, K2, I), tphi_span, final_value_z, options);

%% Display results
% Combine data from both integrations
t_combined = [txy; tz; tphi];
xy_combined = [xy; z; phi];

% Display gains
disp(K1)
disp(K2)

% Plot results
plot2d(t_combined,xy_combined)


%% Functions Below
%
%
%
%

function [K1,K2] = CalculateGain(I, ts_desired, zeta)
% Calculates the gains for each controller
%   

Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);

Axy = [0 ((-Ix - Iz) / Iy); ((-Ix - Iz)/Iy) 0];
Bxy = [(1 / Ix) 0; 0 (1 / Iy)];
Bz = 1 / Iz;

w_n = 4 / (zeta * ts_desired);
tau = ts_desired / 4;

real_pole = - zeta * w_n;
imag_pole = w_n * sqrt(1 - zeta^2);

poles = [real_pole+1j*imag_pole; real_pole-1j*imag_pole];

K1 = place(Axy, Bxy, poles);
K2 = 1 / (Bz * tau);

end

function dxdt = diffEq(t, X0, K, I)
% Integrates transverse velocity proportional gain controller



I_x = I(1,1);
I_y = I(2,2);
I_z = I(3,3);


        x = X0(1:2);
        u = -K * x;
        
        A = [0   ((-I_x-I_z)/I_y); ((-I_x-I_z)/I_y) 0];
        B = [1/I_x 0; 0 1/I_y];  
        
        phi_dot = (cos(X0(6)) * x(2) + sin(X0(6)) * x(1))/(sin(X0(5)));
        theta_dot = (cos(X0(6)) * x(1)) - (sin(X0(6)) * x(2));
        psi_dot = X0(3) - (cos(X0(6)) * x(2) + sin(X0(6)) * x(1)) * cot(X0(5));
       

        dxdt(1,1:2) = A*x + B*u;
        dxdt(3) = 0;
        dxdt(4) = phi_dot;
        dxdt(5) = theta_dot;
        dxdt(6) = psi_dot;


        dxdt = dxdt';



function dxdt = diffEq_z(t, X0, K2, I)
% Integrates axial velocity proportional gain controller

x = X0(3);

r = 0;
u = K2 * (r - x);

I_z = I(3,3);
  
    B = 1/I_z;

    phi_dot = (cos(X0(6)) * X0(2) + sin(X0(6)) * X0(1))/(sin(X0(5)));
    theta_dot = (cos(X0(6)) * X0(1)) - (sin(X0(6)) * X0(2));
    psi_dot = X0(3) - (cos(X0(6)) * X0(2) + sin(X0(6)) * X0(1)) * cot(X0(5));

    dxdt(1,1:2) = 0;
    dxdt(3) = B * u;
    dxdt(4) = phi_dot;
    dxdt(5) = theta_dot;
    dxdt(6) = psi_dot;
     
    dxdt = dxdt';

end




function dxdt = diffEq_phi(t, X0, K2, I)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = X0(4);

r = 0;
u = K2 * (r - x);

I_z = I(3,3);
  
    B = 1/I_z;

    dxdt(1,1:2) = 0;
    dxdt(3) = B*u;
    dxdt(4) = X0(3);
    dxdt(5) = 0;
    dxdt(6) = 0;

    dxdt = dxdt';

end




function [value,isterminal,direction] = eventsFunction(t,x)
    % Example condition: stop if the first two elements of x (xy velocity) are within ±0.02 of zero
    value = all(abs(x(1:2)) <= 0.001) - 1; % Subtracting 1 to change the 'zero crossing' event
    isterminal = 1; % Stop the integration
    direction = 0; % The zero can be approached from any direction
end





function converted_euler = eulerConvert(angles)
% Takes total euler rotation and wraps them in radians to [-pi, pi]
%   Output of integration is total development of the 313 euler
%   rotations, but controller only needs to orient the spacecraft
%   based on the offset between [-pi, pi] rads

converted_euler = mod(angles + pi, 2 * pi) - pi;

end





function plot2d(t, xy)

% Define angular velocity matrices for Wx, Wy, and Wz
Wx = xy(:,1);
Wy = xy(:,2);
Wz = xy(:,3);

%% Plot angular velocities over time
figure(1)
hold on
    
plot(t,Wx,'g');
plot(t,Wy,'b');
plot(t,Wz,'r');
title('Angular velocity components over time' );
    
xlabel('Time(s)');
ylabel('Angular velocity (degrees/sec)');
legend('Wx','Wy','Wz');
    
hold off

%% Plot 313 Euler rotation development (total change)
figure(2)
hold on

plot(t,xy(:,4),'g');
plot(t,xy(:,5),'b');
plot(t,xy(:,6),'r');
title('313 Euler Angle Change Over Time')

xlabel('Time(s)');
ylabel('Euler [rad]');
legend('phi','theta','psi');

hold off
end

%% Plot 313 Euler rotation development (pi -> -pi)

phi = xy(:,4);
theta = xy(:,5);
psi = xy(:,6);

phi_converted = eulerConvert(phi);
theta_converted = eulerConvert(theta);
psi_converted = eulerConvert(psi);

figure(3)
hold on

plot(t,phi_converted(:), 'g');
plot(t,theta_converted(:),'b');
plot(t,psi_converted(:), 'r');
title('converted')

hold off

end


















