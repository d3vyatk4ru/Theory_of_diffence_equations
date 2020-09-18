function [] = paint_YCbI(saddle_vec, point, val_eig)
    % small step from vector
    eps = 0.01;
    %    для селовых точек найти (численно) устойчивое и неустойчивое многообразие
    
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