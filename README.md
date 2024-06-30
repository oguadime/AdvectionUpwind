# AdvectionUpwind

A MATLAB implementation of the advection equation using the upwind scheme.

## Description

This repository contains a MATLAB code that demonstrates the solution of the advection equation using the upwind scheme. The Heaviside function is used to initialize the problem, and the numerical solution is compared to the exact solution at each time step.

## Code Overview

The main function `AdvEqnUpwindb` advances the solution of the advection equation using the upwind scheme. The code includes:
- Initialization of the spatial grid and initial condition.
- A time-stepping loop that advances the solution and plots the numerical and exact solutions at each step.
- Calculation of error norms (L2, L1, and Linf) between the numerical and exact solutions.

### Main Function

```matlab
function AdvEqnUpwindb = AdvEqnUpwindb(M, nu)

format long 

%% set up initial data
A = -0.5; 
B = 1;
dx = (B-A)/M;
x = (A:dx:B); 
uprev = uinit(x); 
unew = 0*uprev; 
t = 0; 
n = 0; 
a = 1;
dt = nu*dx/a;

%% Plot initial condition
plot(x, uprev); 
title(sprintf('Initial condition at t=%g', t)); 
pause; 

while t < 0.5
    %% advance to new time step
    t = t + dt;
    exvec = exfun(x, t);
    for j = 2:M+1
         unew(j) = uprev(j) - nu*(uprev(j) - uprev(j-1));
    end
    unew(1) = 1; %% Left B.C.
 
    %% plot, compare with true solution
    plot(x, unew, 'bo-', x, exvec, 'r'); 
    title(sprintf('Numerical solution at t=%g', t)); 
    axis([-1 3 -0.1 1.5]); 
    legend('numerical', 'exact');
    pause;

    %% calculate error 
    uprev = unew;
    gridL2 = sqrt(dx) * norm(exvec - unew, 2);
    gridL1 = dx * norm(exvec - unew, 1);
    gridLinf = norm(exvec - unew, inf);
end

fprintf('Grid Error L2 = %g\n', gridL2);
fprintf('Grid Error L1 = %g\n', gridL1);
fprintf('Grid Error Linf = %g\n', gridLinf);

end

```
%% Helper Functions
```
function H = SFHeaviside(x)
    H = 0*x; 
    H(x > 0) = 1;
    for j = 1:length(x)        
        if x(j) < 0
            H(j) = 0;
        else
            H(j) = 1;
        end
    end
end

```
%%Initial condition
```
function v = uinit(x)
    v = 1 - SFHeaviside(x);
end
```
%%Exac solution
```
function v = exfun(x, t)
    v = uinit(x - t);
end

```
%% Usage
To run the code, call the AdvEqnUpwindb function with the desired number of grid points M and the Courant number nu. For example:
```
M = 100; % number of grid points
nu = 0.5; % Courant number
AdvEqnUpwindb(M, nu);

```
%%License
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.

