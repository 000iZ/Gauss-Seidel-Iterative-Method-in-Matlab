%{ 
    MTE 204 Project 1
    Gauss-Seidel method for solving square, partially pivoted matrices
    2022/10/8
    By Isaac Zhang
    For Circuit #1
%}


close all;
clear all;
clc;

stopping_num_of_iters = 1000;  % to prevent endless looping
stopping_error = 1e-6;   % convergence criterion, symbol epsilon_s

relaxation_factor = 1;  % to adjust convergence rate, symbol omega
%{
    0.25 -> 124 iterations
    0.50 -> 55 iterations
    0.75 -> 30 iterations
    1.00 -> 318 iterations
    1.25 -> 854 iterations
    1.50 -> Did not converge
%}


% matrix representing unknown currents in the circuit
A = [35 0 -15 0 0 0 -5;
     0 -35 15 0 -7 0 0;
     1 -1 -1 0 0 0 0;
     0 0 0 1 0 -1 0;
     0 1 0 -1 -1 0 0;
     0 0 0 -3 7 -15 0;
     0 0 1 0 1 1 -1];

% solution vector corresponding to the matrix A
b = [180; 0; 0; 0; 0; 0; 0];

% variables names corresponding to the unknown currents
unknowns = ["i_12", "i_23", "i_25", "i_34", "i_35", "i_45", "i_56"];

% solving for the unknown currents using MATLAB's built-in solver
solved = A\b

% number of iterations thus far
num_of_iters = 0;

% the number of rows of the matrix A
num_of_rows = size(A, 1);

% a vector to store guesses, initialized as all zeros
guesses = zeros(1, num_of_rows);

%{
    a vector to store errors, which should all eventually reach below the
    convergence criterion; initialized as all ones
%}
errors = ones(1, num_of_rows);

%{
    convergence is not reached if there is at least one variable whose
    error is more than the convergence criterion; loop until
    convergence or the maximum number of iterations is reached
%}
while (any(errors > stopping_error) ...
            && (num_of_iters < stopping_num_of_iters))
    % store previous guesses
    prev_guess = guesses;

    for i = 1:num_of_rows
        %{
            a summation of elements in a row, multiplied by new or previous
            guesses as specified in the formula for a new guess
        %}
        tally = 0;

        % loop until just before the i-th element, the diagonal
        for j = 1:(i - 1)
            tally = tally + A(i,j) * guesses(j); 
        end

        % loop from just after the i-th element until the last
        for j = (i + 1):num_of_rows
            tally = tally + A(i,j) * prev_guess(j);
        end

        % the new guess, with relaxation
        guesses(i) = (relaxation_factor / A(i,i)) ...
                        * (b(i) - tally) ...
                        + (1 - relaxation_factor) ...
                        * prev_guess(i);

        % update the error for this current element
        errors(i) = abs((guesses(i) - prev_guess(i)) / guesses(i));
    end

    % update the number of iterations thus far
    num_of_iters = num_of_iters + 1;
end

% state if the algorithm converged
if (num_of_iters == stopping_num_of_iters)
    fprintf(['The solution did not converge after %d iterations, ' ...
        'using a relaxation of %.3f \n'], ...
        stopping_num_of_iters, relaxation_factor);
else
    % print to console the guesses for currents
    fprintf('The unknown currents are as follows... \n');
    for i = 1:num_of_rows
        fprintf('\t %s = %.3f \n', unknowns(i), guesses(i));
    end
    
    fprintf(['\nThe above are found after %d iterations of ' ...
        'Gauss-Seidel, using a relxation factor of %.3f \n'], ...
        num_of_iters, relaxation_factor);
end
