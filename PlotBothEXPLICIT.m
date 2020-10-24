%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we calculate the optimal controls for the Merton problem using a
% full grid (not a trinomial tree). We can thus specify T as we please
% provided dt is chosen sufficiently small.
% We search for the optimum over a compact grid of possible control values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic %start timer

K = 1.5;
beta = 0.02;
b    = 0.1;
r    = 0.05;
vol  = 0.3;
g    = 0.5;

I    = 100; %even number
xmax = 100;
h    = xmax/I;

T    = 1;
N    = ceil(K^2*vol^2*I.^2+(r+(b-r)*K+K)*I+beta);   %T=1, I=100, N=5000 takes 36 mins
dt   = T/N;
time = N + 1;

V     = zeros(I+1,time); 
theta = zeros(I+1,time);
con   = zeros(I+1,time);
indextheta = zeros(I+1,time);
indexcon   = zeros(I+1,time);

Vtrue     = zeros(I+1,time); 
thetatrue = zeros(I+1,time);
contrue   = zeros(I+1,time);


Cmax = 150;
thetamax = 150;
Ic = 1500;
It = 1500;
search = 10;
dc = Cmax/Ic;
dth = thetamax/It;
MAXgrid = zeros(It+1,Ic+1);
Disc = exp(-beta*dt)/(1-beta*dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Easy (rough) convergence test


MAXtest = 1/(K^2*vol^2*I^2 + (r+(b-r)*K+K)*I + beta);
[MAXtest,dt]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Fill in the true values
A    = (beta-r*(1-g))/(g) - 0.5*(1-g)*(b-r)^2/((g)^2*vol^2);
for t = time-1:-1:1
    
    ttm = (time-t)*dt;
    gfun = (1+(A-1)*exp(-A*ttm))/A;
    
    for i = 1:1:I+1
         
         contrue(i,t) = (i-1)*h/gfun;
         thetatrue(i,t) = (b-r)*(i-1)*h/(vol^2*(g));
         Vtrue(i,t) = gfun^g*((i-1)*h)^(1-g)/(1-g);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%Terminal
 for i = 2:1:I+1
    V(i,time) = ((i-1)*h)^(1-g)/(1-g); 
 end
 %Side boundary
 BRHS = (I*h)^(1-g)/(1-g);
 for t = time-1:-1:1
     ttm = (time-t)*dt;
     V(I+1,t) = BRHS;
     %exp(-beta*ttm)*BRHS;  
 end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Solve optimisation problem iteratively backwards in time
 for t = time-1:-1:1
     
 
     
     for i = 2:1:I
     
         if t == time-1
             
             for jc = 0:1:Ic
                for jt = 0:1:It

                   Ctry = jc*dc;
                   Thtry = jt*dth;

                   MAXgrid(jt+1,jc+1) = Disc*( ... 
                        V(i-1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*Ctry) +...
                        V(i,t+1)*(1-beta*dt-(dt/h)*((i-1)*h*r+Thtry*(b-r))- (dt/h)*Ctry - (dt/h^2)*Thtry^2*vol^2) + ...
                        V(i+1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*((i-1)*h*r+Thtry*(b-r))) ) + ... 
                        dt*Ctry^(1-g)/(1-g);

                end
             end
             
         else   %Collapse the control search region based on controls from previous time step
             
            for jc = max(indexcon(i,t+1)-search,0):1:min(indexcon(i,t+1)+search,Ic)
                for jt = max(indextheta(i,t+1)-search,0):1:min(indextheta(i,t+1)+search,It)

                   Ctry = jc*dc;
                   Thtry = jt*dth;

                   MAXgrid(jt+1,jc+1) = Disc*( ... 
                        V(i-1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*Ctry) +...
                        V(i,t+1)*(1-beta*dt-(dt/h)*((i-1)*h*r+Thtry*(b-r))- (dt/h)*Ctry - (dt/h^2)*Thtry^2*vol^2) + ...
                        V(i+1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*((i-1)*h*r+Thtry*(b-r))) ) + ... 
                        dt*Ctry^(1-g)/(1-g);

                end
            end
            
         end
         
         
         
         [maxA,ind] = max(MAXgrid(:)); %compute max value
         [m,n] = ind2sub(size(MAXgrid),ind); %compute indices for max value
         
         V(i,t) = maxA;
         theta(i,t) = (m-1)*dth;
         indextheta(i,t) = m-1;
         con(i,t) = (n-1)*dc;
         indexcon(i,t) = n-1;
         
         
         %Convergence check 
         if (1-beta*dt-(dt/h)*((i-1)*h*r+theta(i,t)*(b-r))- (dt/h)*con(i,t) - (dt/h^2)*theta(i,t)^2*vol^2) < 0
            disp('Less Than Zero Funny Business') 
            return %ends program
         end
         
         
         MAXgrid = zeros(It+1,Ic+1);
     
     end
 end
 
 thetaerror = 100.*(theta-thetatrue)./thetatrue;
 thetaerror(isnan(thetaerror)) = 0 ;
 
 conerror = 100.*(con-contrue)./contrue;
 conerror(isnan(conerror)) = 0 ;
 
 toc %end timer
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
 
V     = zeros(I+1,time); 
theta = zeros(I+1,time);
con   = zeros(I+1,time);


Vtrue     = zeros(I+1,time); 
thetatrue = zeros(I+1,time);
contrue   = zeros(I+1,time);


Disc = exp(-beta*dt)/(1-beta*dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Easy (rough) convergence test

K = 1.5;
MAXtest = 1/(K^2*vol^2*I^2 + (r+(b-r)*K+K)*I + beta);
[MAXtest,dt]
if MAXtest < dt
    disp('HEY')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Fill in the true values
A    = (beta-r*(1-g))/(g) - 0.5*(1-g)*(b-r)^2/((g)^2*vol^2);
for t = time-1:-1:1
    
    ttm = (time-t)*dt;
    gfun = (1+(A-1)*exp(-A*ttm))/A;
    
    for i = 1:1:I+1
         
         contrue(i,t) = (i-1)*h/gfun;
         thetatrue(i,t) = (b-r)*(i-1)*h/(vol^2*(g));
         Vtrue(i,t) = gfun^g*((i-1)*h)^(1-g)/(1-g);

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
%Terminal
 for i = 2:1:I+1
    V(i,time) = ((i-1)*h)^(1-g)/(1-g); 
 end
 %Side boundary
 BRHS = (I*h)^(1-g)/(1-g);
 for t = time-1:-1:1
     ttm = (time-t)*dt;
     V(I+1,t) = exp(-beta*ttm)*BRHS;  
 end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Solve optimisation problem iteratively backwards in time
 for t = time-1:-1:1
     
 
     
     for i = 2:1:I
     
         control1 = -(b-r)*h*(V(i+1,t+1)-V(i,t+1))/...
               (vol^2*(V(i+1,t+1)-2*V(i,t+1)+V(i-1,t+1)));
         theta(i,t) = min(max(control1,0),K*(i-1)*h);
         
         control2 = (Disc*(V(i,t+1)-V(i-1,t+1))/(h))^(-1/g);
         con(i,t) = min(max(control2,0),K*(i-1)*h);
         
        Thtry = theta(i,t);
        Ctry  = con(i,t);
         
         
         V(i,t) = Disc*( ... 
                        V(i-1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*Ctry) +...
                        V(i,t+1)*(1-beta*dt-(dt/h)*((i-1)*h*r+Thtry*(b-r))- (dt/h)*Ctry - (dt/h^2)*Thtry^2*vol^2) + ...
                        V(i+1,t+1)*(0.5*(dt/h^2)*Thtry^2*vol^2 + (dt/h)*((i-1)*h*r+Thtry*(b-r))) ) + ... 
                        dt*Ctry^(1-g)/(1-g);
         
         
         
 
         
         
         %Convergence check 
         if (1-beta*dt-(dt/h)*((i-1)*h*r+theta(i,t)*(b-r))- (dt/h)*con(i,t) - (dt/h^2)*theta(i,t)^2*vol^2) < 0
            disp('Less Than Zero Funny Business') 
            theta(i,t)
            con(i,t)
            t
            return %ends program
         end
         
         
        
     
     end
 end
 
 thetaerror2 = 100.*(theta-thetatrue)./thetatrue;
 thetaerror2(isnan(thetaerror2)) = 0 ;
 
 conerror2 = 100.*(con-contrue)./contrue;
 conerror2(isnan(conerror2)) = 0 ;
 
 
 toc
 
 
 
 
 
 
 
 
 
 
 
 
 wealth = h:h:xmax;
 
figure(1)
plot(wealth,conerror2(2:end,1),'k')
hold on
plot(wealth,conerror(2:end,1),'-.k')
title('Consumption   ','FontSize',16)
xlabel('Wealth   ','FontSize',16)
ylabel('Percentage Error   ','FontSize',16)
xlim([0,100])
ylim([-5,2])
h1=legend({'With FOCs  ','Without FOCs  '},'location','northwest');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(h1,'FontSize',14)


figure(2)
plot(wealth,thetaerror2(2:end,1),'k')
hold on
plot(wealth,thetaerror(2:end,1),'-.k')
title('Theta   ','FontSize',16)
xlabel('Wealth   ','FontSize',16)
ylabel('Percentage Error   ','FontSize',16)
xlim([0,100])
ylim([-5,2])
h2=legend({'With FOCs  ','Without FOCs  '},'location','northwest');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(h2,'FontSize',14)
 
 
 
 
 
 
 
 
 
 
 
 
 