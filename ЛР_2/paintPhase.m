function [] = paintPhase(point)

% Фазовый портрет в окрестности положения равновесия
        [X, Y] = meshgrid(point(1) - 1: 0.01 : point(1) + 1, ...
            point(2) - 1 : 0.01 : point(2) + 1);
        
        % my nonlinear system
        U = 3*X + 3*Y + X.*Y +Y.^2;
        V = X + Y.^2;
        
        % paining vector 
        streamslice(X, Y, U, V, 4);
end