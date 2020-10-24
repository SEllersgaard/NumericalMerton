%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As in MunkExplicit.m only here we use the FOCs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

K = 1.5;
beta = 0.02;
b    = 0.1;
r    = 0.05;
vol  = 0.3;
g    = 0.5;

I    = 16; %even number
xmax = 100;
h    = xmax/I;

T    = 0.1;
dt   = T/(I/2);
time = I/2+1;

V     = zeros(I+1,time); 
theta = zeros(I+1,time);
con   = zeros(I+1,time);


Vtrue     = zeros(I+1,time); 
thetatrue = NaN(I+1,time);
contrue   = NaN(I+1,time);






A    = (beta-r*(1-g))/(g) - 0.5*(1-g)*(b-r)^2/((g)^2*vol^2);
start = 1;
for t = time-1:-1:1
    
    start = start + 1;
    ttm = (time-t)*dt;
    
    gfun = (1+(A-1)*exp(-A*ttm))/A;
    
    for i = start:1:(I+2-start)
         
         contrue(i,t) = (i-1)*h/gfun;
         thetatrue(i,t) = (b-r)*(i-1)*h/(vol^2*(g));
         Vtrue(i,t) = gfun^g*((i-1)*h)^(1-g)/(1-g);

    end
end




 
%Terminal
 for i = 2:1:I+1
    V(i,time) = ((i-1)*h)^(1-g)/(1-g); 
 end
 
 Disc = exp(-beta*dt)/(1-beta*dt);
 strt = 1;
 
 for t = time-1:-1:1
     
     strt = strt + 1;
     
     for i = strt:1:(I+2-strt)
     
         control1 = -(b-r)*h*(V(i+1,t+1)-V(i,t+1))/...
               (vol^2*(V(i+1,t+1)-2*V(i,t+1)+V(i-1,t+1)));
         theta(i,t) = min(control1,K*(i-1)*h);
         
         control2 = (Disc*(V(i,t+1)-V(i-1,t+1))/(h))^(-1/g);
         con(i,t) = min(control2,K*(i-1)*h);
         
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
         end
         
         
     
     end
 end
 
 thetaerror = 100.*(theta-thetatrue)./thetatrue;
 thetaerror(isnan(thetaerror)) = 0 ;
 
 conerror = 100.*(con-contrue)./contrue;
 conerror(isnan(conerror)) = 0 ;
 
 conerror(conerror == 0) = NaN;
 thetaerror(thetaerror == 0) = NaN;
 
 
 toc
 