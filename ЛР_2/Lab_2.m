function [] = Lab_2()
    % enter the variables
    syms x y

    s = solve(3*x + 3*y + x*y +y^2 == 0, x + y^2 == 0); 

    % for all stab point
    all = zeros(2, 3);

    disp('Equilibrium positions of the system: ')
    for i = 1:3
        fprintf('x = %d,  y = %d;\n', [s.x(i), s.y(i)]);
        all(1, i) = s.x(i);
        all(2, i) = s.y(i);
    end

    % Finding the Jacobian
    disp('Jacobian of systems: ')
    J = jacobian([3*x+3*y+x*y+y^2, x+y^2],[x, y]);

    for i = 1:3
   
        fprintf('For %d-�� equilibrium  positions:\n', i);
        disp('1) matrix of first approximation:');
    
        % calculating jacobian in points
        subs_J = subs(J, [x y] ,[s.x(i), s.y(i)]);
        disp(subs_J)
        disp('2) eigen values: ')
    
        % finding eigen values
        [v, r] = eig(subs_J);
    
        % display eigen values
        disp(diag(r));
        disp(diag(double(r)));
        
        disp('2) eigen vectors: ');
        disp(v);
    
    % ������� ������� � ����������� ��������� ����������
        figure(i)
        paintPhase(all(:, i))
        hold on
        
        % saving system for saddle point
        if r(1, 1) * r(2, 2) < 0
            paint_YCbI(double(v), double(r));
            xlim([-1, 1])
            ylim([- 1, 1])
        end
    end
end

function [] = paintPhase(point)

        [X, Y] = meshgrid(point(1) - 1: 0.01 : point(1) + 1, ...
            point(2) - 1 : 0.01 : point(2) + 1);
        
        % my nonlinear system
        U = 3*X + 3*Y + X.*Y +Y.^2;
        V = X + Y.^2;
        
        % painting vector 
        streamslice(X, Y, U, V, 4);
end

%  ��� ������� ����� ����� (��������) ���������� � ������������ ������������
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

function [dx] = f(t, x)
    % because MATLAB need a col vector
    dx = [0;0];
    % it's nonlinear system
    dx(1) = 3*x(1) + 3*x(2) + x(1)*x(2) +x(2)^2;
    dx(2) = x(1) + x(2)^2;
end