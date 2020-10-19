function [] = Lab_3 ()
    % equilibrium position (conditional task)
    eq_pos = [0 0];

    % enter symbols variable
    syms x y 
    % Заданная система
    system = [(1+x)*(1-y)-cos(x-y), exp(x-y)-cos(x+y)];
    
    % ----------------- FINDING LINEAR APPROXIMATE ------------------------
    
    % finding jacobian
    disp('Jacoby matrix:');
    J = jacobian(system, [x, y]);
    disp(J);
    
    % change variable into real number(on equilibrium position)
    disp('Linear approximate into [0, 0]:')
    J_subs = subs(J, [x, y], [0, 0]);
    disp(J_subs);
    A = J_subs;
    
    % eigen number
    [V, R] = eig(J_subs);
    disp('Eigen number into [0, 0]: ');
    disp(diag(R));
     
    % finging normalise form
   
    % Order в тейлоре
    
    % It's basis in our space
    basis = [[x^2; 0], [x*y; 0], [y^2; 0], [0; x^2], [0; x*y], [0; y^2]];
    % helpList = [x^2, x*y, y^2];
    tempList = [];
    
    A_x = A * [x; y];
    
    % finding commutator
    for i = 1:6 
       col = basis(:, i);
       
       % calculate coloumn of transform matix from commutator of coloumn
       % from basis (col) and A_x (matrix of linear approximate)
       colMatrix = expand(simplify(jacobian(col, [x y]) * A_x - ...
           jacobian(A_x)*col));
       
       % Расписываем получившийся коммутатор по базисным векторам basis.
       % Дальше записываем получишиеся координаты в матрицу, это и будет 1
       % столбец матрицы перехода. Одна строчка кода --- одна координата 
       % (один элемент столбца).
       
        trMatrix(1, i) = subs(colMatrix(1), [x*y y^2], eq_pos) / x^2;
        trMatrix(2, i) = subs(colMatrix(1), [x^2 y^2], eq_pos) / (x*y);
        trMatrix(3, i) = subs(colMatrix(1), [x^2 x*y], eq_pos) / y^2;
       
        trMatrix(4, i) = subs(colMatrix(2), [x*y y^2], eq_pos) / x^2;
        trMatrix(5, i) = subs(colMatrix(2), [x^2 y^2], eq_pos) / (x*y);
        trMatrix(6, i) = subs(colMatrix(2), [x^2 x*y], eq_pos) / y^2;
    end
    
    % display matrix
    disp('Матрица линейного оператора: ')
    disp(trMatrix);
    
    % finding ker
    
    E6 = eye(6);
    trMatrE = rref([trMatrix E6]);
    
    S = [trMatrix(:, 1:4) E6(:, 1:2)];
    
    disp('Матрица перехода: ')
    disp(S);
    
    % вектор
    f2 = simplify(taylor(system, [x y], 'order', 3) - ...
                            taylor(system, [x y], 'order', 2));
                        
    F2(1) = subs(f2(1), [x*y y^2], eq_pos) / x^2;
    F2(2) = subs(f2(1), [x^2 y^2], eq_pos) / (x*y);
    F2(3) = subs(f2(1), [x^2 x*y], eq_pos) / y^2;
    F2(4) = subs(f2(2), [x*y y^2], eq_pos) / x^2;
    F2(5) = subs(f2(2), [x^2 y^2], eq_pos) / (x*y);
    F2(6) = subs(f2(2), [x^2 x*y], eq_pos) / y^2;
    
    G2 = S\F2';
    G2(1:4) = zeros(4, 1);
    
    g2 = S*G2;
    
    paintPhase(eq_pos);
end

function [] = paintPhase(point)

        [X, Y] = meshgrid(point(1) - 1: 0.01 : point(1) + 1, ...
            point(2) - 1 : 0.01 : point(2) + 1);
        
        % my nonlinear system
        U = (1 + X).*(1 - Y) - cos(X - Y);
        V = exp(X - Y) - cos(X + Y);
        
        % painting vector 
        streamslice(X, Y, U, V, 4);
end