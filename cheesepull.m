
m = 50;    
n = 50;    
L0 = 0.01;    

mass = 4.32 * 10^(-5); 
g = 9.81;  
dt = 0.001; 
N = 10000;

k = 0.1;
mu = 0.01;
pullspeed=0.5;

[x, y] = meshgrid(0:n-1, 0:m-1);
x = L0 * x; 
y = L0 * y;
vx = zeros(m, n); 
vy = zeros(m, n);

x_center = L0 * (n + 1) / 2;  
x = x - x_center;

figure;
hold on;
axis equal;
xlim([-2, 2]);
ylim([-0.1, (n-1)*L0+5]); 
set(gca,'TickLabelInterpreter','latex')
ylabel('height, $h$ (m)','Interpreter', 'latex')
xlabel('radius, $r$ (m)','Interpreter', 'latex')
h=scatter(x(:),y(:),36,'filled');


for step = 1:N
    Fx = zeros(m, n);
    Fy = zeros(m, n) - mass * g;
    for i = 1:m
        for j = 1:n
            for a = -1:1               
                for b = -1:1                 
                    if i+a <= m && j+b <= n && i+a>0 && j+b >0
                        if a ~= 0 && b ~= 0
                            if abs(a) == abs(b)
                                L = sqrt(2)*L0;
                            else
                                L = L0;
                            end
        
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

    vx(2:end-1, :) = vx(2:end-1, :) + (Fx(2:end-1, :) / mass) * dt;
    vy(2:end-1, :) = vy(2:end-1, :) + (Fy(2:end-1, :) / mass) * dt;
    
    x(2:end-1, :) = x(2:end-1, :) + vx(2:end-1, :) * dt;
    y(2:end-1, :) = y(2:end-1, :) + vy(2:end-1, :) * dt;
    
    boundary = y < 0;
    y(boundary) = 0;
    vy(boundary) = 0; 
    
    y(end, :) = y(end, :) + pullspeed * dt;  
    vy(end, :) = pullspeed;
    
    y(1, :) = y(1, :) + 0 * dt;  
    vy(1, :) = 0;
    
    set(h, 'XData', x(:), 'YData', y(:));
    drawnow;

end
 