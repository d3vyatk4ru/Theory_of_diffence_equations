function  [] = draw_pic(A, H, R, N, U)
    % A --- system's matrix
    % H --- height of circle
    % R --- radius of circle
    % N --- number of points in circle for painting
    % U --- transform matrix

    % splitting 
    phi = (1 : N)/ N * 2 * pi; 
    % function in right part
    f = @(t, X) A*X;
    % conditions for built trajectory
    cond = U*[R*sin(phi); R*cos(phi); repmat(H, 1, N)];

    for i = 1:N
    
        [~, oX] = ode45(f, [0.7 0], cond(:, i));

        plot3(oX(:, 1), oX(:, 2), oX(:, 3), 'b');
        hold on
    
    end

end