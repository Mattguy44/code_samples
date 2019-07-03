% ASEN 2012
% Assignment 7
% Problem 3
%
% Author: Matthew Ryan
% Created: 12/3/2017

%% Housekeeping
close; clear; clc;

%% Logistics
filename = 'A7P3output.txt';
fid = fopen(filename, 'w');
fprintf(fid, [  'ASEN 2012\n'...
                'Assignment 7\n'...
                'Problem 3\n'...
                'Author: Matthew Ryan\n'...
                'Date: 12/4/2017\n\n'       ]);

%% Problem 3
% Verify the handwritten solution to the system of equations

% System of equations in the form Ax = b
A = [7 3 -17; -4 0 2; 4 3 -9];
b = [13; -2; -5];
Reps = 1000000;

% Gaussian Elimination
% Note: tic toc is inaccurate for computation times less than 0.1 seconds.
%       Therefore computations have to be repeated a large number of times
%       and averaged for an accurate timing profile.
t01 = tic;
% repeat the operation 100 times
for i = 1:100
    x0 = gaussHouse(A, b);
end
% divide by 100 for an average computation time
t02 = toc(t01)/100;

% Verification using \ operator
t11 = tic;
% repeat the operation Reps number of times
for i = 1:Reps
    x1 = A\b;
end
% divide by Reps for an average computation time
t12 = toc(t11)/Reps;

% Verification using ^-1 operator
t21 = tic;
% repeat the operation Reps number of times
for i = 1:Reps
    x2 = A^-1*b;
end
% divide by Reps for an average computation time
t22 = toc(t21)/Reps;

% Write results to output file
fprintf(fid, 'The Gaussian elimination method took %1.1e seconds and gave:\n\tx = %f\n\ty = %f\n\tz = %f\n\n', t02, x0);
fprintf(fid, 'The \\ operator verification took %1.1e seconds and gave:\n\tx = %f\n\ty = %f\n\tz = %f\n\n', t12, x1);
fprintf(fid, 'The ^-1 operator verification took %1.1e seconds and gave:\n\tx = %f\n\ty = %f\n\tz = %f\n\n', t22, x2);
fprintf(fid, 'The most computationally efficient system of equation solver is the \\ operator.\n');

%% Extra Functions
function x = gaussHouse(A, b)
    if (any(diag(A)) == 0)
        error('Dividing by zero. Rearrange A matrix.');
    end
    
    [row, col] = size(A);
    x = zeros(row,1);
    
    % Elimination
    for i = 1:row-1
        num = A(i,i);
        for j = i+1:row
            den = A(j,i);
            for k = 1:col
                A(j,k) = A(j,k)*num/den - A(i,k);
            end
            b(j) = b(j)*num/den - b(i);
        end
    end
    
    % Back-substitution
    for i = row:-1:1
        x(i) = b(i);
        for j = i+1:col
            x(i) = x(i) - A(i,j)*x(j);
        end
        x(i) = x(i)/A(i,i);
    end
end
            
            