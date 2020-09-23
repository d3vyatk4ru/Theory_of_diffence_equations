
%  для селовых точек найти (численно) устойчивое и неустойчивое многообразие
function [] = paint_YCbI(saddle_vec, val_eig)
    % small step from vector
    eps = 0.01;

    % make normalize vector and next time must do small step in this
    % direction
    for i = 1:2
        saddle_vec_norm(:, i) = eps * saddle_vec(:, i) / norm(saddle_vec(:, i));
        if val_eig(i, i) > 0
            [~, X_new] = ode45(@f, [0 1.4], saddle_vec_norm(:, i));
            plot(X_new(:, 1), X_new(:, 2), 'r');
            [~, X_new] = ode45(@f, [0 1.4], -saddle_vec_norm(:, i));
            plot(X_new(:, 1), X_new(:, 2), 'r');
        elseif val_eig(i, i) < 0
            [~, X_new] = ode45(@f, [10 0], saddle_vec_norm(:, i));
            plot(X_new(:, 1), X_new(:, 2), 'r');
            [~, X_new] = ode45(@f, [7 0], -saddle_vec_norm(:, i));
            plot(X_new(:, 1), X_new(:, 2), 'r');
        end
        hold on     
    end 
end

% Фазовый портрет в окрестности положения равновесия
function [] = paintPhase(point)

        [X, Y] = meshgrid(point(1) - 1: 0.01 : point(1) + 1, ...
            point(2) - 1 : 0.01 : point(2) + 1);
        
        % my nonlinear system
        U = 3*X + 3*Y + X.*Y +Y.^2;
        V = X + Y.^2;
        
        % painting vector 
        streamslice(X, Y, U, V, 4);
end

function [dx] = f(t, x)
    % because MATLAB need a col vector
    dx = [0;0];
    % it's nonlinear system
    dx(1) = 3*x(1) + 3*x(2) + x(1)*x(2) +x(2)^2;
    dx(2) = x(1) + x(2)^2;
end

