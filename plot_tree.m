
% Run MunkExplicit4 first

% dat = conerror;
% mytitle = 'Consumption';

dat = thetaerror;
mytitle = 'Theta';

[D,N] = size(dat);

xdim = zeros(D,N);
ydim = zeros(D,N);

for n = 1:1:N
   xdim(:,n) = n-1; 
end

for d = 1:1:D
    
   ydim(D+1-d,:) = d;
end


wn = ~isnan(dat);


figure('pos',[10 10 900 500])
for n = N-2:-1:1
   
    for d = 2:1:D-1
        
        [n,d]
        
        if wn(d,n) %&& wn(d-1,n) && wn(d+1,n)
                plot([xdim(d,n),xdim(d,n+1)],[ydim(d,n),ydim(d,n+1)],'k');
                hold on
                plot([xdim(d,n),xdim(d,n+1)],[ydim(d,n),ydim(d+1,n+1)],'k');
                hold on
                plot([xdim(d,n),xdim(d,n+1)],[ydim(d,n),ydim(d-1,n+1)],'k');         
        end

        
        
    end
   
end


for n = N-1:-1:1
   
    for d = 2:1:D-1
        if wn(d,n)
            hold on
            s = scatter(xdim(d,n),ydim(d,n),'ok');
            s.MarkerFaceColor = 'k';
            
            strmax = [sprintf('%.2f',dat(d,n)),'%'];
            text(xdim(d,n)+0.1,ydim(d,n),strmax,'HorizontalAlignment','right','Color','r','BackgroundColor','w','EdgeColor','k');
            
            
        end
    end
end


q={'0\delta','1\delta','2\delta','3\delta','4\delta','5\delta','6\delta','7\delta'};
set(gca,'XTickLabel',q);

q2={'15h','14h','13h','12h','11h','10h','9h','8h','7h','6h','5h','4h','3h','2h','1h'};
ax = gca;
yax = 2:1:D-1;
ax.YTick = yax;
set(gca,'YTickLabel',q2);

xlabel('Time Grid','FontSize',16)
ylabel('Wealth Grid','FontSize',16)
title(mytitle,'FontSize',16)




% 
% ax = gca
% ax.Visible = 'off'