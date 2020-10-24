%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the optimal controls for the Merton problem [0,T] using an
% implicit discretisation. Lower boundary taken as 0. Upper boundary 
% handled as though there is zero probability of exiting the grid.
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
 N    = T*20;
global dt;
dt   = T/N;
A    = (beta-r*(1-g))/(g) - 0.5*(1-g)*(b-r)^2/((g)^2*vol^2);
epsilon = 0.0001;
global K;
K = 1.5;
maxrun = 30;

Isize = 50; %[50;100;200;400];
%Isize = [1000;10000;100000;1000000];

for q = 1:1:length(Isize)

    tic
    I    = Isize(q);

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

    %Terminal
    for i = 2:1:I+1
       V(i,N+1) = ((i-1)*h)^(1-g)/(1-g); 
    end
    for i = 2:1:I
       con(i,N+1) = ((V(i,N+1)-V(i-1,N+1))/(h))^(-1/g);
       theta(i,N+1) = -(V(i+1,N+1)-V(i-1,N+1))*(b-r)*(i-1)*h/...
           (2*vol^2*i*(V(i+1,N+1)-2*V(i,N+1)+V(i-1,N+1)));
    end
    % extrapolate the end controls
    con(1,N+1) = 0;
    con(I+1,N+1) = 2*con(I,N+1)-con(I-1,N+1);
    theta(1,N+1) = 0;
    theta(I+1,N+1) = 2*theta(I,N+1)-theta(I-1,N+1);
    P = zeros(I+1,1);




    for n=N:-1:1

        con(:,n) = con(:,n+1);
        theta(:,n) = theta(:,n+1);

        %for i = 2:1:I+1
        %con(i,n) = 10;
        %theta(i,n) = 50;
        %end

        dist = 1;
        k = 0;
        while dist > epsilon && (k<maxrun)
            if k == maxrun-1
               dist 
            end

            k = k+1;

            M(1,1) = 1;
            for i = 2:1:I+1
               M(i,i-1) = 0.5*vol^2*theta(i,n)^2 + h*con(i,n);
               if i < I+1
                   M(i,i+1) = 0.5*vol^2*theta(i,n)^2 + h*(r*(i-1)*h+theta(i,n)*(b-r));
                   M(i,i) = (1-exp(beta*h^2/Q3((i-1)*h)))*Q3((i-1)*h) -  M(i,i-1) -  M(i,i+1) - h^2/dt;
                   P(i,1) = 1-(M(i,i+1)+M(i,i-1)+h^2/dt)/Q3((i-1)*h);
               end
               if P(i,1) < 0
                   disp('Less Than Zero Funny Business') 
               end
            end

            M(I+1,I+1) = (1-exp(beta*h^2/Q3(I*h)))*Q3(I*h) - M(I+1,I) -h^2/dt;

            if 1-(M(I+1,I)+h^2/dt)/Q3(I*h) < 0
                disp('Less Than Zero Funny Business')
            end


            for i = 2:1:I+1
               d(i,1) = -exp(beta*h^2/Q3((i-1)*h))*h^2*con(i,n)^(1-g)/(1-g);
               a(i,1) = -h^2/dt;
            end

            %dist = epsilon/2;

            if k > 1
                Vc = V(:,n);
            end



            %V(:,n) = M\(a.*V(:,n+1)+d);

            V(:,n) = TDMAsolver([0;diag(M,-1)],diag(M),[diag(M,1);0],a.*V(:,n+1)+d);

            if k > 1
                dist = max(abs(V(:,n)-Vc(:)));
            end

            %%%%%Update controls***************************************************

            for i = 2:1:I
               control1 = -(b-r)*h*(V(i+1,n)-V(i,n))/...
                   (vol^2*(V(i+1,n)-2*V(i,n)+V(i-1,n)));
               theta(i,n) = min(max(control1,0),K*(i-1)*h);
            end
            theta(I+1,n) = 0;

            for i = 2:1:I+1
               control2 = (exp(-beta*h^2/Q3((i-1)*h))*(V(i,n)-V(i-1,n))/(h))^(-1/g);
               con(i,n) = min(max(control2,0),K*(i-1)*h);
            end

    %        for i = 2:1:I
    %            control1 = -(b-r)*h*(V(i+1,n)-V(i,n))/...
    %                (vol^2*(V(i+1,n)-2*V(i,n)+V(i-1,n)));
    %            theta(i,n) = min(control1,K*(i-1)*h);
    %        end
    %        theta(1,n) = 2*theta(2,n)-theta(3,n);
    %        theta(I+1,n) = 2*theta(I,n)-theta(I-1,n);
    % 
    %         for i = 2:1:I
    %            %control2 = (exp(-beta*h^2/Q((i-1)*h))*(V(i,1)-V(i-1,1))/(h*g))^(1/(g-1));
    %            control2 = (exp(-beta*h^2/Q3((i-1)*h))*(V(i,n)-V(i-1,n))/(h))^(-1/g);
    %            con(i,n) = min(control2,K*(i-1)*h);
    %         end
    %         con(1,n) = 2*con(2,n)-con(3,n);
    %         con(I+1,n) = 2*con(I,n)-con(I-1,n);



        end
        iteration(N+1-n,1) = k; 

    end

    
    

    thetaerror = 100.*(theta-thetatrue)./thetatrue;
     thetaerror(isnan(thetaerror)) = 0 ;

     conerror = 100.*(con-contrue)./contrue;
     conerror(isnan(conerror)) = 0 ;



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
    title('Theta   ','FontSize',16)
    xlabel('Wealth   ','FontSize',16)
    ylabel('Percentage Error   ','FontSize',16)
    xlim([0,100])
    ylim([-5,2]) 


    figure(1)
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
    title('Consumption   ','FontSize',16)
    xlabel('Wealth   ','FontSize',16)
    ylabel('Percentage Error   ','FontSize',16)
    xlim([0,100])
    ylim([-5,2])



     hold on
    
    toc

end


figure()
Verror = 100.*(V-Vtrue)./Vtrue;
plot(Verror(:,1))



figure()
    time = T:-dt:0;
    wealth = 0:h:xmax;
    
    surf(time,wealth,contrue)
    hold on
    surf(time,wealth,con)
    colormap(white)
    title('Consumption   ','FontSize',16)
    ylabel('Wealth   ','FontSize',16)
    xlabel('Time to expiry   ','FontSize',16)
    
    
    figure()
    
    surf(time,wealth,thetatrue)
    hold on
    surf(time,wealth,theta)
    colormap(white)
    title('Theta   ','FontSize',16)
    ylabel('Wealth   ','FontSize',16)
    xlabel('Time to expiry  ','FontSize',16)
    
