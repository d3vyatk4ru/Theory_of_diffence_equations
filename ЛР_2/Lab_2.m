function [] = Lab_2()
    % enter the variables
    syms x y

    s = solve(3*x + 3*y + x*y +y^2 == 0, x + y^2 == 0); 

    % for all stab point
    all = zeros(2, 3);

    disp('��������� ���������� �������: ')
    for i = 1:3
        fprintf('x = %d,  y = %d;\n', [s.x(i), s.y(i)]);
        all(1, i) = s.x(i);
        all(2, i) = s.y(i);
    end

    % Finding the Jacobian
    disp('������� �������: ')
    J = jacobian([3*x+3*y+x*y+y^2, x+y^2],[x, y]);

    for i = 1:3
   
        fprintf('��� %d-�� ��������� ����������:\n', i);
        disp('1) ������� ������� ������� �����������:');
    
        % calculating jacobian in points
        subs_J = subs(J, [x y] ,[double(s.x(i)) double(s.y(i))]);
        disp(subs_J)
        disp('2) c���������� ��������: ')
    
        % finding eigen values
        [v, r] = eig(double(subs_J));
    
        % display eigen values
        disp(diag(r));
    
    % ������� ������� � ����������� ��������� ����������
        figure(i)
        paintPhase(all(:, i))
        hold on
        
        % saving system for saddle point
        if r(1, 1) * r(2, 2) < 0
            paint_YCbI(v, all, r);
            xlim([-1, 1])
            ylim([- 1, 1])
        end
    end

    
end