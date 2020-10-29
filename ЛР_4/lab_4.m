function [] = lab_4 ()
    
    % enter perid variable
    global T;
    % enter the variables
    syms x y
    
    % finding equilibrium positions
    s = solve(-x^2 + 17/9*x*y - 9/10*y^2 + 40/9*x - 19/3*y + 8/15 == 0,...
              -x^2 + 7/4*x*y - y^2 + 4*x - 9/2*y + 3/2 == 0); 

    % for all stab point
    all = zeros(2, 2);
    
    for i = 1:2
        fprintf('x = %d,  y = %d;\n', [s.x(i), s.y(i)]);
        all(1, i) = double(s.x(i));
        all(2, i) = double(s.y(i));
    end
    
    % ������ �������� ������ ���� ��������� ����������, ��� ���
    % ������ �������� ������ � ��� ������������� �� �����.
    % enter global variable because dont use paremetrs in fucntion
    global equilibriumY; 
    global equilibriumX;
    
    equilibriumX = all(1, 1);
    equilibriumY = all(2, 1);
    
    % Puankare();
    % painting Phase trajectory
    figure
    paintPhase([equilibriumX, equilibriumY]);
    hold on;
    
    % painting cicle and calculate period
    [X, Y] = build_curv();
    
    % �� ���������� ���������� ����������, ������ ���
    % ������� �� ��������� �� �������.
    % painting equilibrium position 
    plot(all(1, 1), all(2, 1) ,'o');

    multiplic_cicle(X, Y);
    
end

function [] = paintPhase(point)

        [X, Y] = meshgrid(point(1) -8: 0.05 : point(1) + 8  , ...
            point(2) - 8 : 0.05 : point(2) + 8);
        
        % my nonlinear system
        U = -X.^2 + (17/9)*X.*Y - (9/10)*Y.^2 + (40/9)*X - (19/3)*Y + 8/15;
        V = -X.^2 + (7/4)*X.*Y -Y.^2 + 4*X - (9/2)*Y + 3/2;
        
        % painting vector 
        streamslice(X, Y, U, V, 10);
end

function [loc_X, loc_Y] = build_curv()

    % feature of global variable in MATLAB
    global T;
    global equilibriumX;
    global equilibriumY
    
    step = 1e-2;
    eps = 1e-4;
    h = (0:step:5);
    % make cycle curve
    
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @stopInt);
    [~, X] = ode45(@f, h, [equilibriumX + 0.1; equilibriumY], opts);
    
    while abs(X(1, 1) - X(size(X, 1), 1)) >= eps
        
        % enter setting for curve
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', @stopInt);
        
        % enter start point for integrate in ode45
        cond = [X(size(X, 1), 1); equilibriumY];
        
        % integration
        [t, X] = ode45(@f, h, cond, opts);
       
    end
    
    % painting if accuracy lower the epsilon
    plot(X(:, 1) , X(:, 2), 'r', 'LineWidth', 1);
    disp('������ �����:');
    T = t;
    disp(t(size(t)));
    
    % enter limits
    xlim([equilibriumX - 4 equilibriumX + 6]);
    ylim([equilibriumY - 4 equilibriumY + 6]);
    
    loc_X = X(:, 1);
    loc_Y = X(:, 2);
end

function [value, isterminal, direction] = stopInt(t, X)
     % feature of global variable in MATLAB
    global equilibriumY;
    
    delta = 1;
    
    % + delta *(t < 1) --- ��� ���� ����� �� ���� ��������� ��������������
    % �� 1-2 ����, � ������ �������� � ��.
    value = X(2) - equilibriumY + delta *(t < 0.1);
    
    % ������������� ������ ��� �������, ���� ���������� ��� � ����� �����.
    % ��� -1 ��������� ������� ��� ����������� ������ ���� ��� �.
    isterminal = 1;    
    direction = 1;

end

function [] = multiplic_cicle(X, Y)
    
    global T
    
    syms x y
    J = jacobian([-x^2 + 17/9*x*y - 9/10*y^2 + 40/9*x - 19/3*y + 8/15 ...
              -x^2 + 7/4*x*y - y^2 + 4*x - 9/2*y + 3/2], [x y]);
    
    C = eye(2);
    for i = 1:size(T) - 1
        A_lastStep = double(subs(J, [x y] ,[X(i) Y(i)]));
        %
        A_nextStep = double(subs(J, [x y], [X(i + 1) Y(i + 1)]));
        
        A = (A_nextStep + A_lastStep) / 2;
        
        C = expm(A * (T(i + 1) - T(i))) * C;
    end
    
    [~, E] = eig(C);
    
    disp(E);

end

function [dx] = f(t, x)
    % because MATLAB need a col vector
    dx = [0;0];
    % it's nonlinear system
    dx(1) = -x(1)^2 + 17/9*x(1)*x(2) - 9/10*x(2)^2 + 40/9*x(1) - 19/3*x(2) + 8/15;
    dx(2) = -x(1)^2 + 7/4*x(1)*x(2) - x(2)^2 + 4*x(1) - 9/2*x(2) + 3/2;
    
end

function [] = Puankare()

    global equilibriumY;
    global equilibriumX;

% ������ ����������� ��������
step = 0.01; 
eps = 1e-4; 
k = 0; 
X_p = [];
Y_p = [];

for j = (equilibriumX) : step : (equilibriumX + 3)
    opts = odeset('RelTol', 1e-7,'AbsTol',1e-7,'Events', @stopInt);
    
    [~, X] = ode45(@f,(0:step:1),[j; equilibriumY], opts);
    
    k = k + 1;
    X_p(k) = j;
    Y_p(k) = X(size(X, 1), 1);
    if Y_p(k) - X_p(k) < 0.01
        temp = X_p(k);
    end
    % ������� �����, � ����������� ������� ��������� ����
end

plot(X_p - equilibriumX, Y_p - X_p);
hold on
% plot(temp, 0,'o','LineWidth', 1);
title('����������e ��������');
grid on

end