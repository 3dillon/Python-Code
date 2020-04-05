
Radius = 0.2;
FreeStreamVelocity = 5.0;
LiftPerUnitSpan = 00;
Density = 1.18;
EPS = 1*10^-6

Lambda = 0;
Gamma = LiftPerUnitSpan/ (Density*FreeStreamVelocity);
Cl = LiftPerUnitSpan/ (0.5*Density*FreeStreamVelocity^2*Radius);

Xmin = -1;
Xmax = 1;
numX = 101;

Ymin = -1;
Ymax = 1;
numY = 101;

dx = (Xmax - Xmin) / (numX - 1);
dy = (Ymax - Ymin) / (numY - 1);


for ii = 1:numX
  for jj = 1:numY
    Xvals(ii, jj) = Xmin + dx * (ii - 1);
    Yvals(ii, jj) = Ymin + dy * (jj - 1);
    r = (Xvals(ii, jj)^2 + Yvals(ii, jj)^2);
    r = max(r, EPS);
    r = sqrt(r);
    theta = atan2(Yvals(ii,jj), Xvals(ii,jj));
    
    Phi(ii,jj) = FreeStreamVelocity*r*cos(theta)*(1+(Radius/r)^2);
    Psi(ii,jj) = FreeStreamVelocity*r*sin(theta)*(1-(Radius/r)^2);
    Vr = (1-(Radius/r)^2)*FreeStreamVelocity*cos(theta);
    Vtheta = -(1+(Radius/r)^2)*FreeStreamVelocity*sin(theta);
    
    Phi(ii,jj) = Phi(ii,jj) - Gamma/ (2.*pi)*theta;
    Psi(ii,jj) = Psi(ii,jj) + Gamma/ (2.*pi)*log(r/Radius);
    Vr = Vr + 0;
    Vtheta = Vtheta - Gamma/ (2.*pi*r);
    
    Phi(ii,jj) = Phi(ii,jj) + Lambda/(2*pi)*log(r/Radius);
    Psi(ii,jj) = Psi(ii,jj) + Lambda/ (2 * pi) * theta;
    
    Vr = Vr + Lambda / (2*pi*r);
    Vtheta = Vtheta + 0;
    
    if(r<Radius)
      Vmag(ii,jj) = -1;
      Cp(ii,jj) = 10;
    else
      Vmag(ii,jj) = sqrt(Vr^2 + Vtheta^2);
      Cp(ii,jj) = 1 - (Vmag(ii,jj) / FreeStreamVelocity)^2;
    endif

  end
endfor


figure(1) 
clf(1)  
hold on  
[c1,h1]=contour(Xvals(:,:),Yvals (:,:), Psi(:,:),'k','LineWidth',4,[-10:1:10]); %plots the x vs t  
clabel (c1, h1, -10:2:10, "fontsize", 24); 
 
[c2,h2]=contour(Xvals(:,:),Yvals (:,:), Cp(:,:),'g','LineWidth',4,[-2:0.1:1]); %plots the x vs t  
clabel (c2, h2, -2:0.2:1, "fontsize", 24); 
 
hold off 
t1=title('Coplot of Stream Function (\psi) and C_p') 
x1=xlabel('Position, x (m)') %x-axis label  
y1=ylabel('Position, y (m)') %y-axis label  
set(t1,'FontSize',24) %change fontsize 
set(x1,'FontSize',24) %change fontsize 
set(y1,'FontSize',24) %change fontsize 
 
%Now only plot the Vel Mag, what does it look like 
figure(2) 
clf(2)  
hold on  
[c2,h2]=contour(Xvals(:,:),Yvals (:,:), Vmag(:,:),'k','LineWidth',0.5,[0:0.5:8]); %plots the x vs t  
clabel (c2, h2, 0:1:8, "fontsize", 24); 
 
t1=title('Velocity Mag(|V|)') 
x1=xlabel('Position, x (m)') %x-axis label  
y1=ylabel('Position, y (m)') %y-axis label  
set(t1,'FontSize',24) %change fontsize 
set(x1,'FontSize',24) %change fontsize 
set(y1,'FontSize',24) %change fontsize 
 
%Now only plot the Cp, what does it look like 
figure(3) 
clf(3)  
hold on  
[c2,h2]=contourf(Xvals(:,:),Yvals (:,:), Cp(:,:),'k','LineWidth',0.5,[-2:0.1:1]); %plots the x vs t  
clabel (c2, h2, -2:0.2:1, "fontsize", 24); 
 
t1=title('Pressure Coefficent(C_p)') 
x1=xlabel('Position, x (m)') %x-axis label  
y1=ylabel('Position, y (m)') %y-axis label  
set(t1,'FontSize',24) %change fontsize 
set(x1,'FontSize',24) %change fontsize 
set(y1,'FontSize',24) %change fontsize 
 
%Now one more, a coplot of Psi and Phi  
figure(4) 
clf(4)  
hold on  
[c1,h1]=contourf(Xvals(:,:),Yvals (:,:), Psi(:,:),'k','LineWidth',4,[-10:1:10]); %plots the x vs t  
clabel (c1, h1, -10:2:10, "fontsize", 24); 
 
[c2,h2]=contour(Xvals(:,:),Yvals (:,:), Phi(:,:),'r','LineWidth',4,[-10:1:10]); %plots the x vs t  
clabel (c2, h2, -10:2:10, "fontsize", 24); 
 
 
t1=title('Coplot of Stream Function (\psi) and Velocity Potential (\phi)') 
x1=xlabel('Position, x (m)') %x-axis label  
y1=ylabel('Position, y (m)') %y-axis label  
set(t1,'FontSize',24) %change fontsize 
set(x1,'FontSize',24) %change fontsize 
set(y1,'FontSize',24) %change fontsize 
 
display(r);
display(theta);
 