
function [] = lab_1(A)
    % system matrix's in numerical form
    A = [-7 9 -4;
         -4 6 -4;
        1 -1 -2;];
 
    % changing numerical form to symbols
    As = sym(A);

% ------------------------- finding solution ------------------------------

    % adding variables type syms  for dsolve(...)&solve function
    syms x(t) y(t) z(t) 

    % general solution
    vars = [x(t); y(t); z(t)];
    s = dsolve(diff(vars) == A * vars);

    disp('General solution: ')

    disp(s.x);
    disp(s.y);
    disp(s.z);

    % ------------------------- finding the singular point  -------------------
    
    % enter syms variable for solve system
    syms x y z

    disp('Точка покоя: ');

    s = solve(A *[x; y; z] == 0);

    disp(s.x);
    disp(s.y);
    disp(s.z);

    % ------------------------- finding the eigenvectors and eigenvalues ------
    [V, D] = eig(As);

    % display to ComWin
    disp('Еigen vectors: ');
    disp(V);

    % display to ComWin
    disp('Еigen values: ');
    disp(D);

    % creating transform matrix
    disp('Transform matrix: ');
    disp(V);

    % number of split point
    N = 100; 
    % height
    H = 100;  
    % radius
    R = 10;                  

    draw_pic(A, H, R, N, double(V));
    hold on
    
    % painting line
    fplot3(V(1, 3) *t, V(2, 3)*t, V(3, 3), [-60 60], 'r')
    
    hold on
    
    % painting plane
    
    % creating grid
    [u, v] = meshgrid(-40:40);
    
    X = u*V(1, 1) + v*V(1, 2);
    Y = u*V(2, 1) + v*V(2, 2);
    Z = u*V(3, 1) + v*V(3, 2);
    
    plot3(X, Y, Z, 'y');
    plot3(0, 0, 0, '-o', 'Color', 'black')
    
    grid on
end
