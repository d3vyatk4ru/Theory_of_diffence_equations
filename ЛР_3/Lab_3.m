function [] = Lab_3 ()
    % equilibrium position (conditional task)
    bal_pos = [0, 0];

    % enter symbols variable
    syms x y 
    
    % ----------------- FINDING LINEAR APPROXIMATE ------------------------
    
    % finding jacobian
    disp('Jacoby matrix:');
    J = jacobian([(1+x)*(1-y)-cos(x-y), exp(x-y)-cos(x+y)], [x, y]);
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
   
    % Order � �������
    
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
       
       % ����������� ������������ ���������� �� �������� �������� basis.
       % ������ ���������� ����������� ���������� � �������, ��� � ����� 1
       % ������� ������� ��������. ���� ������� ���� --- ���� ���������� 
       % (���� ������� �������).
       for j = 1:6
            tempList(j) = subs(colMatrix(1), [x*y y^2], [0 0]) / x^2;
            tempList(j) = subs(colMatrix(1), [x^2 y^2], [0 0]) / (x*y);
            tempList(j) = subs(colMatrix(1), [x^2 x*y], [0 0]) / y^2;
       
            tempList(j) = subs(colMatrix(2), [x*y y^2], [0 0]) / x^2;
            tempList(j) = subs(colMatrix(2), [x^2 y^2], [0 0]) / (x*y);
            tempList(j) = subs(colMatrix(2), [x^2 x*y], [0 0]) / y^2;
       end
    end
    
    trMatrix(:, i) = tempList;
    
    % display matrix
    disp('Transform matrix: ')
    disp(trMatrix);
    
    % finding ker
    
    E6 = eye(6);
    traMatrE = rref([trMatrix E6]);
    
    S = [trMatrix(:, 1:4) E6(:, 1:2)];
    
    disp(S);
end