function [nc,A1,B1]=nec_cond_at_ddp(ds,d12s,d21,r1,r2,a1,b2,b1,a2,Lx)

us=(r1*b2-r2*b1)/(a1*b2-a2*b1);
vs=(r2*a1-r1*a2)/(a1*b2-a2*b1);
%disp([us,vs])

d_s=ds;%(sqrt(105)-5)/(32*pi^2);
gamma_s=d12s;%(477+79*sqrt(105))/(2*(45-sqrt(105))*pi^2);

k=pi/Lx;
%M_m=[-(d+gamma*vs)*m^2*k^2-a1*us, -(gamma*m^2*k^2+b1)*us; -a2*vs, -(d*m^2*k^2+b2*vs)];
m=1;
M1=[-(d_s+gamma_s*vs)*m^2*k^2-a1*us, -(gamma_s*m^2*k^2+b1)*us; -d21*vs*m^2*k^2-a2*vs, -((d_s+d21*us)*m^2*k^2+b2*vs)];
%trM1=trace(M1);
m=2;
M2=[-(d_s+gamma_s*vs)*m^2*k^2-a1*us, -(gamma_s*m^2*k^2+b1)*us; -d21*vs*m^2*k^2-a2*vs, -((d_s+d21*us)*m^2*k^2+b2*vs)];
%trM2=trace(M2);

T1=[M1(1,2) M1(1,1); -M1(1,1) M1(2,1)];
T2=[M2(1,2) M2(1,1); -M2(1,1) M2(2,1)];
detT1=det(T1);
detT2=det(T2);

%p0=-b1;
p1=-(gamma_s*k^2+b1);
p2=-(gamma_s*4*k^2+b1);
%p3=-(gamma_s*9*k^2+b1);
%p4=-(gamma_s*16*k^2+b1);
p1v=-(d21*k^2+a2);
p2v=-(d21*4*k^2+a2);

f1=[-2*a1*T1(1,1)*T2(1,1)+p1*(T1(1,1)*T2(2,1)+T2(1,1)*T1(2,1)), -a1*T1(1,1)^2+p2*T1(1,1)*T1(2,1)];
g1=[-2*b2*T1(2,1)*T2(2,1)+p1v*(T1(1,1)*T2(2,1)+T2(1,1)*T1(2,1)), -b2*T1(2,1)^2+p2v*T1(1,1)*T1(2,1)];


A1=1/detT1*(T1(2,2)*f1(1)-T1(1,2)*g1(1));
%A1_loro=13*(324175+62301*sqrt(105))/150880;
B1=1/detT2*(T2(2,2)*f1(2)-T2(1,2)*g1(2));
%B1_loro=13*(4323445-424489*sqrt(105))/104960;
% necessary condition
nc=sign(A1*B1);

end