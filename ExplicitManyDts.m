tic

beta = 0.02;
b    = 0.1;
r    = 0.05;
vol  = 0.3;
g    = 0.5;
K = 1.5;


Isize = [50;100;200;400];

for q = 1:1:4

I    = Isize(q); %even number
xmax = 100;
h    = xmax/I;

T    = 1;
N    = ceil(K^2*vol^2*I.^2+(r+(b-r)*K+K)*I+beta);   %T=1, I=100, N=5000 takes 36 mins
dt   = T/N;
time = N + 1;

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
     V(I+1,t) = BRHS;
     %exp(-beta*ttm)*BRHS;  
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
 
 thetaerror = 100.*(theta-thetatrue)./thetatrue;
 thetaerror(isnan(thetaerror)) = 0 ;
 
 conerror = 100.*(con-contrue)./contrue;
 conerror(isnan(conerror)) = 0 ;
 
 toc %end timer
 
 wealth = h:h:xmax;
 
figure(2)
if q == 1
    plot(wealth,thetaerror(2:end,1),'k')
elseif q == 2
    plot(wealth,thetaerror(2:end,1),'--k')
elseif q == 3
    plot(wealth,thetaerror(2:end,1),'-.k')
else
    plot(wealth,thetaerror(2:end,1),':k')
end
title('Theta   ','FontSize',16)
xlabel('Wealth   ','FontSize',16)
ylabel('Percentage Error   ','FontSize',16)
xlim([0,100])
ylim([-5,2]) 
hold on
 
 
% figure(1)
% if q == 1
%     plot(wealth,conerror(2:end,1),'k')
% elseif q ==2
%     plot(wealth,conerror(2:end,1),'--k')
% elseif q == 3
%     plot(wealth,conerror(2:end,1),'-.k')
% else
%     plot(wealth,conerror(2:end,1),':k')
% end
% title('Consumption   ','FontSize',16)
% xlabel('Wealth   ','FontSize',16)
% ylabel('Percentage Error   ','FontSize',16)
% xlim([0,100])
% ylim([-5,2])
% hold on



end
h1=legend({'I = 50  ','I = 100  ','I = 200  ','I = 400  '},'location','northwest');
set(h1,'FontSize',14)

set(gca, 'YGrid', 'on', 'XGrid', 'off')


% figure(2)
% plot(wealth,thetaerror(2:end,1))
% title('Theta   ','FontSize',16)
% xlabel('Wealth   ','FontSize',16)
% ylabel('Percentage Error   ','FontSize',16)
% xlim([0,100])
% ylim([-10,5])
 
 