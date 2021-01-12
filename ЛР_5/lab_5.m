function [] = lab_5()

clc;
clear all;

% Finding equilibrium position
syms x y z;

s = solve(2.7*y + z == 0, -x + y^2 == 0, x + y == 0);

disp('Положения равновесия:');
disp([double(s.x), double(s.y),  double(s.z)]);

% Writing nonzero equilibrium position in list
nonZeroEquilibrium = [];
len = size(s.x);
for i = 1 : len(1)
    if (s.x(i) ~= '0') & (s.y(i) ~= '0') & (s.z(i) ~= '0')
        nonZeroEquilibrium(1) = double(s.x(i));
        nonZeroEquilibrium(2) = double(s.y(i));
        nonZeroEquilibrium(3) = double(s.z(i));      
    end
end

%----------------------- Check the existence of attractor -----------------

% small step
eps = 1e-3;

% conditionals for start point in integration 
% There is a point, which remain in infinite small distance from  equilibrium position
cond = [nonZeroEquilibrium(1) - eps, nonZeroEquilibrium(2) - eps, nonZeroEquilibrium(3) - eps];

[t, X] = ode45(@f, [0:0.01:1500], cond);

% attractor visualisation
figure;
plot3(X(: ,1),X(:, 2),X(:, 3));
rotate3d on;
grid on;

% Finding the Jacobian
F = [2.7*y + z; -x + y^2; x + y];

% getting symbolic Jacobian
J = jacobian(F, [x, y, z]);

% transfom symbolic Jacobian to matlab function
J_matFunc = matlabFunction(J, 'Vars', [x, y, z]); 
df = @(x) J_matFunc(x(1), x(2), x(3));

%----------------------- Curve and Lyapunov's values ----------------------

L = list_of_index(df, X, t);

figure;
hold on;
plot(t(2:end), L(1, :));
plot(t(2:end), L(2, :), 'color', 'r');
plot(t(2:end), L(3, :), 'color', 'black');
grid on;

valueL = get_Lyapunov_value(3, L);
disp('Показатели Ляпунова:')
disp(valueL);

end

function [dx] = f(t, x)
    % because MATLAB need a col vector
    dx = [0; 0; 0];
    
    % it's my system
    dx(1) = 2.7*x(2) + x(3);
    dx(2) = -x(1) + x(2)^2;
    dx(3) = x(1) + x(2);
end

function [L] = list_of_index(df, X, t)
    
    % enter all list for finding Lyapunov's index.
    A_i = [];
    eta = eye(3);
    gamma = X;
    eta1col = eta(:, 1);
    eta2col = eta(: ,2);
    eta3col = eta(:, 3);
    
    for i = 2 : size(t, 1)
        % computing matrix A in [t_i, t_{i+1}]
        A_i = (df(gamma(i, :)) + df(gamma(i - 1, :))) / 2;
        eta = expm(A_i * (t(i) - t(i - 1))) * eta;
    
        eta1col = [eta1col eta(:, 1)];
        eta2col = [eta2col eta(:, 2)];
        eta3col = [eta3col eta(:, 3)];
    end
    
    L1 = [];
    L2 = [];
    L3 = [];

    for i = 2 : size(t, 1)
        L1 = [L1 log(norm(eta1col(:, i))) / t(i)];
        L2 = [L2 log(norm(eta2col(:, i))) / t(i)];
        L3 = [L3 log(norm(eta3col(:, i))) / t(i)];
    end
    
    L = [L1; L2; L3];
end

function [low] = low_value(indexL)

    N = length(indexL);
    s = ceil(0.01 * N);
    high = [];

    for j = 1 : N - s
        high = [high max(indexL(j:N))];
    end
    low = min(high(1:N - s));
end

function [Lyapunov_value] = get_Lyapunov_value(shape, index)

    for i = 1:shape
        Lyapunov_value(i) = low_value(index(i, :));
    end
end


