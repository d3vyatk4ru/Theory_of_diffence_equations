
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
    disp('1) ������� ������� ������� �����������:')
    subs_J = subs(J, [x y] ,[s.x(i) s.y(i)]);
    disp(subs_J)
    disp('2) c���������� ��������: ')
    r = eig(subs_J);
    fprintf('%d %d\n\n', [r(1,1), r(2,1)]);
end
