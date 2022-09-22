clc
clear
close all
Nx =100;      % Value of N
h = 1/Nx;   
PI = 3.141592654;
time_step = 0.1
method = 'BTCS'
total_iteration = 1000/(time_step);
time = 0.0;
c = 1.0;
b = 0.01;
mu = c/b;

%Declaration
u_exact = zeros(Nx+1,1);
u= zeros(Nx+1,1);
u_new = zeros(Nx+1,1);
x = linspace(0,1,Nx+1);


%Exact Solution at t = INf
for i= 1:Nx+1
    u_exact(i) =  ( 1 / (exp(mu) - 1) ) * (exp(mu * x(i)) - 1 );
end

 % Initialise forplot frames
 frame(total_iteration) =  struct('cdata',[],'colormap',[]);


alpha  = c * time_step / (2* h);
beta = b*time_step/(h*h);

%%% FTCS %%
aa = zeros(Nx-1,1);
bb = zeros(Nx-1,1);
cc = zeros(Nx-1,1);
dd = zeros(Nx-1,1);
% Matrix Formation
for i = 1:Nx-1
bb(i) = 1 + (2*beta);
cc(i) =  alpha -  beta ;
aa(i) = - ( alpha + beta ) ;
end

%Boundary Conditions
u(1) = 0;
u(Nx+1) = 1;
u_new(1) = 0;
u_new(Nx+1) = 1;

iteration_count = 0;
temp_time = 1;
  count  = 1;
  format long;
while (iteration_count < total_iteration)
    
  time = time + time_step;
  
dd = appendBoundary(u,aa,cc,Nx);
u_new = thomasSolver(aa,bb,cc,dd,u_new,Nx-2);


  
if(iteration_count == temp_time) 
    figure(count);
    plot(x,u_exact);
    hold on
    plot(x,u_new);
    grid on;
    ylim([-1.0,1.5]);
    legend('Steady State Solution ' , ' Computed Solution');
    title([' BTCS : h = ',  num2str(h), '  \Delta (t)  : ',num2str(time_step),'-  Time = ' , num2str(time) , ' sec '  ]) 
    temp_time = temp_time *10;
    drawnow;
    error = 0.0;
for i = 1:Nx+1
    error = error + (u_exact(i) - u_new(i))^2 ;
end

error_norm = 0.0;
error_norm =  sqrt(error)/(Nx+1)
    frame(count) = getframe(gcf);
    count = count+1;
    hold off
end



 u = u_new;
 iteration_count = iteration_count +1 ;

    
    
end


%Ploting Exact Solution
plot(x,u_exact);
hold on
plot(x,u_new);
    grid on;
legend('Steady State Solution ' , ' Computed Solution');
    title([' BTCS : h = ',  num2str(h), ',  \Delta(t)  : ',num2str(time_step),'sec ,  Time : ' , num2str(time) , ' sec '  ]) 
ylim([-1.0,1.5]);
     %%% Plot Animation 
% fig = figure;
% movie(fig,frame,10)

%%% Error
error = 0.0;
for i = 1:Nx+1
    error = error + (u_exact(i) - u_new(i))^2 ;
end

error_norm = 0.0;
error_norm =  sqrt(error)/(Nx+1)


function dd = appendBoundary(u,aa,cc,Nx)
 dd = zeros(Nx-1,1); 
    for i = 1:Nx-1
        dd(i) = u(i+1);
    end
    
  dd(1) = dd(1)  - (aa(1)* 0.0);
  dd(Nx-1) = dd(Nx-1)- (cc(Nx-1)*1.0);
    
  
end




function u_new =  thomasSolver(a,b,c,d,u_new,Nx)

% b - Leadin Diagonal )
% a Sub Diagonal
% c Super Diagonal 
% d Right hand Side force vector
% a(i)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i)
p = zeros(Nx+1,1);
q = zeros(Nx+1,1);


p(1) = c(1)/b(1);
q(1) = d(1)/b(1);

for i = 2:Nx+1  
    p(i) = c(i) / (b(i) - (a(i) * p(i-1))) ; 
    q(i) = (d(i) -  (a(i) * q(i-1)) ) / ( b(i) - (a(i) * p(i-1)) ) ;
 end
 
 u_new(Nx+1+1) = q(Nx+1) ;
 
 for i = Nx:-1:1
    u_new(i+1) = -(p(i)*u_new(i+1+1)) + q(i) ;
 end
end
