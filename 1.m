k = 10;
A = [0 1 0 0 0 0;
    -2*k 0 k 0 0 0;
    0 0 0 1 0 0;
    k 0 -2*k 0 k 0;
    0 0 0 0 0 1;
    0 0 k 0 -2*k 0];

b = [
    0 0 0
    1 0 0
    0 0 0
    0 1 0
    0 0 0
    0 0 1];

H_little = [0 1; 0 -1; 1 0; -1 0];

H = [H_little, zeros(4, 2), zeros(4, 2);
    zeros(4, 2),H_little,zeros(4, 2);
    zeros(4, 2), zeros(4, 2), H_little];

g = ones(12,1).*0.1;

t_0 = 0;
% t_f = 10;

t_f = 4.5;

N = 100; 
h = (t_f - t_0)/N;

x = zeros(6, N);
x(:,1) = 1;
x_0 = x(:,1);

F = @(t)expm(A*t);
% 
% for i = 1:N
%     x(:, i+1) = F(h) * x(:,i);
% end

for_d = @(t)H*F(t_f-t)*b;
get_d_h = @(s)integral(for_d,s,s+h,'ArrayValued', true);
d_h_values = zeros(12, 3, N);

for i = 1:N
    d_h_values(:,:,i) = get_d_h(t_0 + i*h-h);
end

Ale = zeros(12, 6*N);
for i = 1:N
    Ale(:,6*i - 5 : 6*i) = [d_h_values(:,:,i),-d_h_values(:,:,i)];
end

g_wave = g - H * F(t_f - t_0)* x_0;

c = ones(1, 6*N);

ub = 1*ones(6*N, 1);

z_and_v = linprog(c,Ale,g_wave,[],[],zeros(6*N, 1),ub);

%u for object of control number one, two and three
u_one = zeros(1, N); 
u_two = zeros(1, N);
u_three = zeros(1, N);
for i = 1:N 
    u_one(i) = z_and_v(6*i - 5) - z_and_v(6*i - 2);
    u_two(i) = z_and_v(6*i - 4) - z_and_v(6*i - 1);
    u_three(i) = z_and_v(6*i - 3) - z_and_v(6*i);
end

u = zeros(3, N);
for i = 1:N 
    u(:, i) = [u_one(i);u_two(i); u_three(i)];
end

for i = 1:N
    curr = t_0 + (i-1)*h;
    next = t_0 + i*h;
    x(:, i+1) = F(h) * x(:,i) +  integral(@(t) F(next - t)*b,curr, next, 'ArrayValued', true)*u(:, i);
end

figure('Name','x1_x2','NumberTitle','off');
plot(x(1, :), x(2, :));
grid on;

figure('Name','u_one','NumberTitle','off');
stairs(u_one);
grid on;
