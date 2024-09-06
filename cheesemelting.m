
m = 10;    
n = 150;    
L0 = 0.001;    

mass = 5.23*10^(-4); 
g = 9.81;  
dt = 0.0001; 
N = 10000;

k_max = 4000;
k_min = 2000; 
mu_max = 0.01;
mu_min = 0.001; 

beta = 0.0001;

[x, y] = meshgrid(0:n-1, 0:m-1);
x = L0 * x; 
y = L0 * y;
vx = zeros(m, n); 
vy = zeros(m, n);
T = zeros(m,n);
x_centre = L0 * (n + 1) / 2;  
x = x - x_centre;
T(:, 1) = 100;      
T(:, end) = 100;     
T(1, :) = 100;       
T(end, :) = 100; 

figure;
hold on;
axis equal;
xlim([-0.15, 0.15]);
ylim([0-0.01, (m-1)*L0+0.02]);
 
set(gca,'TickLabelInterpreter','latex')
ylabel('height, $h$ m','Interpreter', 'latex')
xlabel('radius, $r$ m','Interpreter', 'latex')
cb = colorbar;
cb.Label.String = 'temperature, $T^\circ \mathrm{C}$';
cb.Label.Interpreter = 'latex';
set(cb,'TickLabelInterpreter','latex')
cmap = [
     linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(0, 0, 256)';  ];
  colormap(cmap);
 caxis([40 80]);
 h=scatter(x(:),y(:),36, T(:), 'filled');


for step = 1:N
    T_new = T;
    Fx = zeros(m, n);
    Fy = zeros(m, n) - mass * g;
    for i = 1:m
        for j = 1:n
            if (2 <= i) && (i <= m-1)
                if (2<= j) && (j <= n-1)
                    dx = T(i+1,j)-T(i,j);
                    dy = T(i,j+1)-T(i,j);
                    T_new(i, j) = T(i, j) + beta * (dt/L0^2) * ((T(i+1, j) - 2*T(i, j) + T(i-1, j))+ (T(i, j+1) - 2*T(i, j) + T(i, j-1)) );
                    if T_new(i,j) > 100
                        T_new(i,j) = 100;
                        T_new(i,j) = 100;
                    end
                end
            end
            for a = -1:1               
                for b = -1:1                 
                    if i+a <= m && j+b <= n && i+a>0 && j+b >0
                        if a ~= 0 && b ~= 0
                            if abs(a) + abs(b) == 2
                                L = sqrt(2)*L0;
                            else
                                L = L0;
                            end
                            T_avg = (T(i, j) + T(i+a, j+b)) / 2;
                            k = k_max - (k_max - k_min) * (T_avg / 100);
                            mu = mu_max - (mu_max - mu_min) * (T_avg / 100);
                            dx = x(i+a,j+b) - x(i,j);
                            dy = y(i+a,j+b) - y(i,j);
                            dist = sqrt(dx^2 + dy^2);
                            springForce = k * (dist - L) / dist;
                     
                            Fx(i,j) = Fx(i,j) + springForce * dx ;
                            Fy(i,j) = Fy(i,j) + springForce * dy ;
                            
                            dvx = vx(i+a, j+b) - vx(i, j);
                            dvy = vy(i+a, j+b) - vy(i, j);
                            Fx(i, j) = Fx(i, j) + mu*dvx;
                            Fy(i, j) = Fy(i, j) + mu*dvy;
                        end
                    end
                end
            end
        end
    end
    T(:, 1) = 100;      
    T(:, end) = 100;     
    T(1, :) = 100;       
    T(end, :) = 100; 
    T = T_new;

    vx = vx + (Fx / mass) * dt;
    vy = vy + (Fy / mass) * dt;
    
    x = x + vx*dt;
    y = y + vy*dt;
    
    boundary = y < 0;
    y(boundary) = 0;
    vy(boundary) = 0; 

   set(h, 'XData', x(:), 'YData', y(:),'CData', T(:));
   drawnow;
end
