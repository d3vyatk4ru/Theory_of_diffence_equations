function [] = Lab_3 ()
    % palance posotion (conditional task)
    bal_pos = [0, 0];

    % enter symbols variable
    syms x y 
    
    % ----------------- FINDING LINEAR APPROXIMATE ------------------------
    
    % finding jacobian
    disp('Jacoby matrix:');
    J = jacobian([(1+x)*(1-y)-cos(x-y), exp(x-y)-cos(x+y)], [x, y]);
    disp(J);
    
    % change variable into real number(on balance position)
    disp('Linear approximate into [0, 0]:')
    J_subs = subs(J, [x, y], [0, 0]);
    disp(J_subs);
    
    % eigen number
    [V, R] = eig(J_subs);
    disp('Eigen number into [0, 0]: ');
    disp(diag(R));
    
end