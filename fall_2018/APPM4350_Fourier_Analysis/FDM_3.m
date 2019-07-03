% APPM 4350 - Fourier Analysis
% Finite Difference Method 2
% Approximate solution to a PDE of the form:
% d(v)/dt = d^2(v)/dx^2 - v + H(v-theta) + J
%
% Author: Matthew Ryan
% Created: Dec. 6, 2018

% Output:   u = solution approximation
% Input:    x = set of x points to evaluate at
%           t = set of t points to evaluate at
%           BC1 = first boundary condition function handle (f(0,t))
%           BC2 = second boundary condition function handle (f(end,t))
%           IC = Initial condition function handle (f(x,0))
%           inh = inhomogenous function handle (f(x,t))
function u = FDM_3(x, t, BC1, BC2, IC, inh, theta)
    nx = length(x); % number of x points
    nt = length(t); % number of t points
    h = (x(end) - x(1))/(nx-1); % x step
    q = (t(end) - t(1))/(nt-1); % t step
    
    u = zeros(nx,nt);
    u(:,1) = feval(IC, x); % initial conditions
    u(1,:) = feval(BC1, t); % first x boundary conditions
    u(end,:) = feval(BC2, t); % last x boundary conditions
    
    v = feval(inh, x, t); % forcing function
    
    h_sq = h*h;
    k = .5*q/h_sq;
    a_1 = 2*k + 1;
    a_2 = 1 - 2*k - q;
    

    A = zeros(nx-2,nx);
    for ii = 1:(nx-2)
        A(ii,ii) = -k;
        A(ii,ii+1) = a_1;
        A(ii,ii+2) = -k;
    end
    
    B = zeros(nx-2,nx);
    for ii = 1:(nx-2)
        B(ii,ii) = k;
        B(ii,ii+1) = a_2;
        B(ii,ii+2) = k;
    end
    
    for jj = 1:(nt-1) % time loop
        u(:,jj+1) = A\(B*u(:,jj) + q*heaviside(u(2:end-1,jj) - theta*(1 + 0.5*cos(x(2:end-1)'))) + q*v(2:end-1,jj));
    end
end
