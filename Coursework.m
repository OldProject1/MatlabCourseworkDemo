function forGlobalArea

clc, close all

% data for OC problem

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

L = 1;

% discretization

N_initial = 100;

h = t_f/N_initial;

% initial position

x_0 = ones(n,1);

t_0 = 0;

Ale = Form_LP(N_initial);

[u_res, x_res] = P(t_0, x_0);

showTrajectories(x_res);

showControls(u_res);

%----------------------------------------------
% functions
%--------------------------------------------------
function Ale = Form_LP(N)
        
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
end

%--------------------------------------------------
function [u_opt_pos, x_opt_pos] = P(tau, z)
    
    N = round((t_f - tau)/h);%new N for new  tau

    g_wave = g - H * F(t_f - tau) * z;

    c = ones(1, 2*r*N);

    ub = L*ones(2*r*N, 1);

    z_and_v = linprog(c,Ale,g_wave,[],[],zeros(2*r*N, 1),ub);
    
    zv = reshape(z_and_v, 2*r, N);
    u  = zv(1:r,:) - zv(r+1:2*r,:); % r /times N
    
    x = trajectory(z, tau, t_f, u);
    
    u_opt_pos = u;
    x_opt_pos = x;
end
%--------------------------------------------------
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
%--------------------------------------------------
function showTrajectories(tr)
    figure('Name','x','NumberTitle','off'); 
    for i = 1:n/2
        subplot(1, 3, i);
        plot(tr(2*i - 1, :), tr(2*i, :), 'Linewidth', 1); %hold on; 
        grid on;
    end
end
%--------------------------------------------------
function showControls(ctrl)
    figure('Name','u','NumberTitle','off');
    for i = 1:r
        subplot(1, 3, i);
        stairs(ctrl(i,:), 'Linewidth', 1);
        ylim([-L*1.1, L*1.1])
        grid on;
    end
end
end

