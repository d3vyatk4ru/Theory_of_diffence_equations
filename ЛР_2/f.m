function [dx] = f(t, x)
    % because we MATLAB need a col vector
    dx = [0;0];
    % it's nonlinear system
    dx(1) = 3*x(1) + 3*x(2) + x(1)*x(2) +x(2)^2;
    dx(2) = x(1) + x(2)^2;
end