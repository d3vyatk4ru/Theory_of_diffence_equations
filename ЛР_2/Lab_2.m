
% enter the variables
syms x y

s = solve(3*x + 3*y + x*y +y^2 == 0, x + y^2 == 0); 

disp('��������� ���������� �������: ')
for i = 1:3
    fprintf('x = %d,  y = %d;\n', [s.x(i), s.y(i)]);
end

%Finding the Jacobian
disp('������� �������: ')
J = jacobian([3*x+3*y+x*y+y^2, x+y^2],[x, y]);

for i = 1:3
   
    fprintf('��� %d-�� ��������� ����������:\n', i);
    disp('1) ������� ������� ������� �����������:');
    
    % calculating jacobian in points
    subs_J = subs(J, [x y] ,[s.x(i) s.y(i)]);
    disp(subs_J)
    disp('2) c���������� ��������: ')
    
    % finding eigen values
    r = eig(subs_J);
    fprintf('%d %d\n\n', [r(1,1), r(2,1)]);
end


n = 100;
m = 100;
p = 100;

[X, Y, Z] = meshgrid(1:n, 1:m, 1:p);
streamslice(X, Y, Z, 3*X + 3*Y + X.*Y +Y.*Y);
