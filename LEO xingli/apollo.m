function xdot=apollo(t,x) 

U0=1/82.45; 
U1=1-U0; 
R1=sqrt((x(1)+U0)^2+x(3)^2); 
R2=sqrt((x(1)-U1)^2+x(3)^2); 

xdot=[x(2);
    2*x(4)+x(1)-U1*(x(1)+U0)/R1^3;
    x(4);
    -2*x(2)+x(3)-U1*x(3)/R1^3-U0*x(3)/R2^3];