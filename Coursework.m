function forGlobalArea

n = 6; r = 3; m = 2*n;

t_f = 4.5;

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

F = @(t)expm(A*t);

H = [eye(n); (-1).*eye(n)];

g = ones(m, 1).*0.1;

N_initial = 100;

h = t_f/N_initial;

z = ones(n,1);

tau = 0;

[u_res, x_res] = P(tau, z);

%showTrajectories(x_res);

% showControls(u_res);

function [u_opt_pos, x_opt_pos] = P(tau, z)
    t_0 = tau;
    x_0 = z;

    N = round((t_f - t_0)/h);%new N for new  t_0

    for_d = @(t)H*F(t_f-t)*b;
    get_d_h = @(s)integral(for_d,s,s+h,'ArrayValued', true);
    d_h_values = zeros(m, r, N);

    for i = 1:N
        d_h_values(:,:,i) = get_d_h(t_0 + i*h-h);
    end

    Ale = [];
    for i = 1:N
        Ale = [Ale  d_h_values(:,:,i) -d_h_values(:,:,i)];
    end

    g_wave = g - H * F(t_f - t_0)* x_0;

    c = ones(1, n*N);

    ub = ones(n*N, 1);

    z_and_v = linprog(c,Ale,g_wave,[],[],zeros(n*N, 1),ub);
    
    zv = reshape(z_and_v, n, N);
    u  = zv(1:3,:) - zv(4:6,:);% 3 /times N
    
    x = trajectory(x_0, t_0, t_f, u);
    
    u_opt_pos = u;
    x_opt_pos = x;
end
function x = trajectory(x0, t_begin, t_end, u)
    N = (t_end - t_begin)/h;
    x = zeros(n, N);
    x(:,1) = x0;
    for j = 1:N
        curr = t_begin + (j-1)*h;
        next = t_begin + j*h;
        x(:, j+1) = F(h) * x(:,j) +  integral(@(t) F(next - t)*b,curr, next, 'ArrayValued', true)*u(:, j);
    end
end
function showTrajectories(tr)
    figure('Name','x','NumberTitle','off'); 
    for i = 1:r
        plot(tr(2*i - 1, :), tr(2*i, :)); hold on; 
    end
    grid on;
end
function showControls(ctrl)
    for i = 1:r
        figure('Name','u','NumberTitle','on');
        stairs(ctrl(i,:));
        grid on;
    end
end
end

