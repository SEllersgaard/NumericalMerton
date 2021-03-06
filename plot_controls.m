%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




global beta;
beta = 0.02;
b    = 0.1;
r    = 0.05;
vol  = 0.3;
global g;
g    = 0.5;
global I;
 T    = 1;
 N    = T*30;
global dt;
dt   = T/N;
A    = (beta-r*(1-g))/(g) - 0.5*(1-g)*(b-r)^2/((g)^2*vol^2);
epsilon = 0.0001;
global K;
K = 1.5;
maxrun = 30;







    I    = 50;

    xmax = 100;
    h    = xmax/I;
   
   

    con   = zeros(I+1,N+1);
    theta = zeros(I+1,N+1);
    V     = zeros(I+1,N+1);
    contrue   = zeros(I+1,N+1);
    thetatrue = zeros(I+1,N+1);
    Vtrue     = zeros(I+1,N+1);
    iteration = zeros(N,1);

    M         = zeros(I+1,I+1);
    d         = zeros(I+1,1);
    a         = zeros(I+1,1);

    for n = 1:1:N+1
        for i = 1:1:I+1
             gfun = (1+(A-1)*exp(-A*(T-(n-1)*dt)))/A;
             contrue(i,n) = (i-1)*h/gfun;
             thetatrue(i,n) = (b-r)*(i-1)*h/(vol^2*(g));
             Vtrue(i,n) = gfun^g*((i-1)*h)^(1-g)/(1-g);

        end
    end
    
    
    
    
    figure()
    time = T:-dt:0;
    wealth = 0:h:xmax;
    
    surf(time,wealth,contrue)
    colormap(white)
    title('Consumption   ','FontSize',16)
    ylabel('Wealth   ','FontSize',16)
    xlabel('Time to expiry   ','FontSize',16)
    
    
    figure()
    
    surf(time,wealth,thetatrue)
    colormap(white)
    title('Theta   ','FontSize',16)
    ylabel('Wealth   ','FontSize',16)
    xlabel('Time to expiry  ','FontSize',16)
    