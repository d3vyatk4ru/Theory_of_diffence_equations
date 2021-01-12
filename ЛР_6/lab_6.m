
% enter lambda function
p = @(C) 1 - C^2;
q = @(C) (C^2 - 1)/C - 2.7*C - C^2;

beta = @(C) -q(C)^2/(4*p(C));

rho = @(y, z, C) beta(C) - (C^2 - 1)/C*y - C*z;


% enter grid for y, z and C
y = -2.1:0.05:2.5;
z = -2.1:0.05:2.5;
[Y, Z] = meshgrid(y, z);

C = (-1 : 0.1 : 1);

P1 = zeros(size(Y, 1), size(Z, 1));

for i = 1:size(Y, 1)
   for j = 1:size(Z, 1)
      max = -100000;
      for k = 1:length(C)
        if rho(Y(i), Z(j), C(k)) > max
            max = rho(Y(i, 1), Z(j, 1), C(k));
        end
      end
    P1(i, j) = max;
   end
end

C = [-20:0.1:-1 - 0.1, 1 + 0.1:0.1:20];
P2 = zeros(size(Y, 1), size(Z, 1));

for i = 1:size(Y, 1)
   for j = 1:size(Z, 1)
     min = 1000000;
      for k = 1:length(C)
        if rho(Y(i), Z(j), C(k)) < min
            min = rho(Y(i), Z(j), C(k));
        end
      end
      P2(i, j) = min;
   end
end


for i = 1:size(Y, 1)
   for j = 1:size(Z, 1)
        if P1(i, j) < P2(i, j)
            P1(i, j) = NaN;
            P2(i, j) = NaN;
        end
   end
end


% disp(matrix);
mesh(P1,  Y, Z);
hold on
mesh(P2,  Y, Z);
hold on

shading interp
alpha .6

% small step
eps = 1e-3;

% conditionals for start point in integration 
% There is a point, which remain in infinite small distance from  equilibrium position
cond = [0 - eps, 0 - eps, 0 - eps];

[t, X] = ode45(@f, [0:0.01:1500], cond);

% attractor visualisation
plot3(X(: ,1),X(:, 2),X(:, 3), 'black');
rotate3d on;
grid on;
hold on

plot3(0, 0, 0, 'or');
plot3(1, -1, 2.7, 'or');




function [dx] = f(t, x)
    % because MATLAB need a col vector
    dx = [0; 0; 0];
    
    % it's my system
    dx(1) = 2.7*x(2) + x(3);
    dx(2) = -x(1) + x(2)^2;
    dx(3) = x(1) + x(2);
end