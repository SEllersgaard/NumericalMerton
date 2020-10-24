clc
clear



beta = 0.02;
r    = 0.05;
g    = 0.5;
lambda = 0.3;
kappa = 3;
thelr = 0.3^2;
eta   = 1;


 T    = 1;
 N    = T*10;

dt   = T/N;
epsilon = 0.0001;
K = 1.5;
maxrun = 30; %30;


%%

khat = kappa-((g-1)/g)*eta*lambda;
qcond = khat^2 + ((g-1)/g^2)*lambda^2*eta^2;
if qcond < 0
   disp('Quadratic constraint not satisfied')
   return
end
nu   = sqrt(qcond);

% A1 function
A1FACT = (lambda^2/g);
A1 = @(tau) A1FACT*(exp(nu*tau)-1)./((khat+nu)*(exp(nu*tau)-1)+2*nu);

% A0 function
A0FACT = kappa*thelr*g/((g-1)*eta^2);
A0 = @(tau) r*tau - A0FACT*((nu+khat)*tau + 2*log(abs( 2*nu./((khat+nu)*(exp(nu*tau)-1)+2*nu) )));

% gtilde function
gtilde = @(v,t,s) exp(-beta*(s-t)/g-((g-1)/g)*( A0(s-t)+v.*A1(s-t) )); 

% zeta function
zeta = @(v,t) ( integral( @(s) exp(-beta*(s-t)/g-((g-1)/g)*( A0(s-t)+v*A1(s-t) )),t,T ) + ...
                 gtilde(v,t,T))^(-1);
             
%pi function
pifun = @(v,t) lambda/g + ((g-1)/g)*eta*zeta(v,t)...
    *( integral( @(s) A1(s-t).*exp(-beta*(s-t)/g-((g-1)/g)*( A0(s-t)+v*A1(s-t) )),t,T ) + ...
                 A1(T-t).*gtilde(v,t,T));
             



             


% figure()
% tau = 0:0.01:1;
% subplot(2,1,1)
% plot(tau,A1(tau),'k')
% subplot(2,1,2)
% plot(tau,A0(tau),'k')
% subplot(3,1,3)
% plot(t,gtilde(sqrt(0.3),t))

%%

Isize = [50;100;200;400];

%Isize = [50;400];

%Isize = 100;


for q = 1:1:length(Isize)

    I    = Isize(q);

    xmax = 1;
    h    = xmax/I;
   
    % Q function
    Qvol = @(x) h^2*(beta -(1-g)*(r+K*lambda*x-K) + 0.5*(1-g)*g*K^2*x + (1/dt) + ...
                (1/h)*abs(kappa*(thelr-x)-K*eta*x*(1-g)) + (1/h^2)*eta^2*x );
    
    Qvol2 = @(x,p,c) h^2*(beta -(1-g)*(r+p*lambda*x-c) + 0.5*(1-g)*g*p^2*x + (1/dt) + ...
                (1/h)*abs(kappa*(thelr-x)-p*eta*x*(1-g)) + (1/h^2)*eta^2*x );        
   

    con       = zeros(I+1,N+1);
    theta     = zeros(I+1,N+1);
    V         = zeros(I+1,N+1);
    contrue   = zeros(I+1,N+1);
    thetatrue = zeros(I+1,N+1);
    Vtrue     = zeros(I+1,N+1);
    iteration = zeros(N,1);
    betatilde = zeros(N,1);

    M         = zeros(I+1,I+1);
    d         = zeros(I+1,1);
    a         = zeros(I+1,1);

    for n = 1:1:N+1
        for i = 1:1:I+1
            
            t = (n-1)*dt;
            v = (i-1)*h;

             contrue(i,n) = zeta(v,t);
             thetatrue(i,n) = pifun(v,t);
             Vtrue(i,n) = (1/contrue(i,n))^g;

        end
    end

    %Terminal
    for i = 1:1:I+1
       V(i,N+1) = 1; 
    end
    for i = 2:1:I+1
       con(i,N+1) = (V(i,N+1))^(-1/g);
       theta(i,N+1) = lambda/g - (eta/g)*(1/V(i,N+1))*(V(i,N+1)-V(i-1,N+1))/h;
    end
    % extrapolate the end controls
    con(1,N+1) = (V(1,N+1))^(-1/g);
    theta(1,N+1) = theta(2,N+1);
    P = zeros(I+1,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do iterations backwards in time
    
    
    
    for n=N:-1:1

        con(:,n) = con(:,n+1); %con(:,n+1);
        theta(:,n) = theta(:,n+1); %theta(:,n+1);


        dist = 1;
        k = 0;
        while dist > epsilon && (k<maxrun)
            if k == maxrun-1
               dist 
            end

            k = k+1;

            for i = 1:1:I+1
                
                v = (i-1)*h;
                pi = theta(i,n);
                xi = con(i,n);
                
                betatilde(i) = beta - (1-g)*(r+pi*lambda*v - xi) + 0.5*(1-g)*g*pi^2*v;
                
                % Calculate M(i,i-1)
                if i>1
                   M(i,i-1) = h*(kappa*v + pi*eta*v*(1-g)) + 0.5*eta^2*v;
                   Mim = M(i,i-1);
                else
                   Mim = 0;
                end
                
                % Calculate M(i,i+1)
                if i < I+1
                   M(i,i+1) = h*kappa*thelr + 0.5*eta^2*v;
                   Mip = M(i,i+1);
                else
                   Mip = 0; 
                end
                
                M(i,i) = (1-exp(betatilde(i)*h^2/Qvol2(v,pi,xi)))*Qvol2(v,pi,xi) -  Mim -  Mip - h^2/dt;
               

               d(i,1) = -exp(betatilde(i)*h^2/Qvol2(v,pi,xi))*h^2*xi^(1-g);
               a(i,1) = -h^2/dt;
            
            end

            %dist = epsilon/2;

            if k > 1
                Vc = V(:,n);
            end



            V(:,n) = M\(a.*V(:,n+1)+d);

            %V(:,n) = TDMAsolver([0;diag(M,-1)],diag(M),[diag(M,1);0],a.*V(:,n+1)+d);

            if k > 1
                dist = max(abs(V(:,n)-Vc(:)));
            end

            %%%%%Update controls***************************************************

            for i = 2:1:I+1
               control1 = lambda/g - (eta/g)*(1/V(i,n))*(V(i,n)-V(i-1,n))/(h);
               theta(i,n) = min(max(control1,0),K);
            end
            theta(1,n) = theta(2,n);
            theta(I+1,n) = theta(I+1,n);
            
            for i = 1:1:I+1
               control2 = (V(i,n))^(-1/g);
               con(i,n) = min(max(control2,0),K);
            end





        end
        iteration(N+1-n,1) = k; 

    end

    
    thetaerror = 100.*(theta-thetatrue)./thetatrue;
     thetaerror(isnan(thetaerror)) = 0 ;

     conerror = 100.*(con-contrue)./contrue;
     conerror(isnan(conerror)) = 0 ;


     Verror = 100.*(V-Vtrue)./Vtrue;
     Verror(isnan(Verror)) = 0 ;
     
     
     
     
     

     wealth = h:h:xmax;

   
     
    if q == 1
        figure(1)
        plot(wealth,thetaerror(2:end,1),'k')
        hold on
    elseif q == 2
        figure(1)
        plot(wealth,thetaerror(2:end,1),'--k')
        hold on
    elseif q == 3
        figure(1)
        plot(wealth,thetaerror(2:end,1),'-.k')
        hold on
    else
        figure(1)
        plot(wealth,thetaerror(2:end,1),':k')
        h1=legend({'I = 50  ','I = 100  ','I = 200  ','I = 400  '},'location','northwest');
        set(gca, 'YGrid', 'on', 'XGrid', 'off')
        set(h1,'FontSize',14)
    end
    title('Pi   ','FontSize',16)
    xlabel('Variance   ','FontSize',16)
    ylabel('Percentage Error   ','FontSize',16)
    xlim([0,xmax])
    ylim([0,1.6]) 


    if q == 1
         figure(2)
        plot(wealth,conerror(2:end,1),'k')
        hold on
    elseif q ==2
        figure(2)
        plot(wealth,conerror(2:end,1),'--k')
        hold on
    elseif q == 3
        figure(2)
        plot(wealth,conerror(2:end,1),'-.k')
        hold on
    else
        figure(2)
        plot(wealth,conerror(2:end,1),':k')
        h2=legend({'I = 50  ','I = 100  ','I = 200  ','I = 400  '},'location','northwest');
        set(gca, 'YGrid', 'on', 'XGrid', 'off')
        set(h2,'FontSize',14)
    end
    title('Zeta   ','FontSize',16)
    xlabel('Variance   ','FontSize',16)
    ylabel('Percentage Error   ','FontSize',16)
    xlim([0,xmax])
    ylim([0,1.6])


    
    if q == 1
         figure(3)
        plot(wealth,Verror(2:end,1),'k')
        hold on
    elseif q ==2
        figure(3)
        plot(wealth,Verror(2:end,1),'--k')
        hold on
    elseif q == 3
        figure(3)
        plot(wealth,Verror(2:end,1),'-.k')
        hold on
    else
        figure(3)
        plot(wealth,Verror(2:end,1),':k')
        h2=legend({'I = 50  ','I = 100  ','I = 200  ','I = 400  '},'location','northwest');
        set(gca, 'YGrid', 'on', 'XGrid', 'off')
        set(h2,'FontSize',14)
    end
    title('Value Function   ','FontSize',16)
    xlabel('Variance   ','FontSize',16)
    ylabel('Percentage Error   ','FontSize',16)
    xlim([0,xmax])
    
    

     hold on
    
    
    
end



% figure()
% subplot(2,1,1)
% surf(contrue)
% subplot(2,1,2)
% surf(thetatrue)

% figure()
% Verror = 100.*(V-Vtrue)./Vtrue;
% plot(Verror(:,1))
% 
% 
% figure()
% plot(theta(:,1))
% hold on
% plot(thetatrue(:,1))
% 
% 
% figure()
% plot(con(:,1))
% hold on
% plot(contrue(:,1))
% 
% 
% figure()
% plot(V(:,1))
% hold on
% plot(Vtrue(:,1))




