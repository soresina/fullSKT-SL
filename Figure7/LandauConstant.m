function [sigma,L]=LandauConstant(c1,c2,a1,a2,b2,mu1,mu2,gamma11,gamma22,gamma12,gamma21,Gamma,b_c,ks_c)
% calcola la costante L
wos=gamma11*gamma22-gamma12*gamma21;
%%
% homogeneous equilibrium
u0=(mu1*gamma22-mu2*gamma12)/wos;
v0=(mu2*gamma11-mu1*gamma21)/wos;
% Jacobian of the reaction at (u0,v0)
K=[-gamma11*u0, -gamma12*u0; -gamma21*v0, -gamma22*v0];
%trK=trace(K);
%detK=det(K);
% Jacobian of the diffusion at (u0,v0)
D=[c1+2*a1*u0+b_c*v0, b_c*u0; b2*v0, c2+2*a2*v0+b2*u0];
%trD=trace(D);
detD=det(D);
%% Dispersion relation
%q=gamma11*u0*(2*a2*v0+c2)+gamma22*v0*(2*a1*u0+c1)+b_c*v0*(gamma22*v0-gamma21*u0)+b2*u0*(gamma11*u0-gamma12*v0);
%h=@(ks) detD*ks^2+Gamma*q*ks+Gamma^2*detK;
%disp_rel=@(lambda,ks) lambda.^2+(ks*trD-Gamma*trK)*lambda+h(ks);
%ks_c=-0.5*Gamma*q/detD;
%k_c=sqrt(ks_c);
%alpha=v0*(gamma21*u0-gamma22*v0);
%beta=gamma11*u0*(2*a2*v0+c2)+gamma22*v0*(2*a1*u0+c1)+b2*u0*(gamma11*u0-gamma12*v0);
%q=-alpha*b+beta
% equation 2.9
%A=0.25*alpha^2/detK;
%B=v0*(2*a2*v0+c2);
%C=beta/alpha*v0*(2*a2*v0+c2)+(2*a1*u0+c1)*(2*a2*v0+b2*u0+c2);
%xi1=(B+sqrt(B^2+4*A*C))/(2*A);
%xi2=(B-sqrt(B^2+4*A*C))/(2*A);
%if xi1>0
%    xi_plus=xi1;
%else
%    xi_plus=xi2;
%end
%b_c=beta/alpha+xi_plus;
%disp(b_c);
% b_c=5.296980278451512
% and corresponding critical wavenumber
% k_c=3.206914937366216
%% Weakly Nonlinear Analysis
M=(-D(2,1)*ks_c+Gamma*K(2,1))/(D(2,2)*ks_c-Gamma*K(2,2));
rho=[1;M];
Ms=(-D(1,2)*ks_c+Gamma*K(1,2))/(D(2,2)*ks_c-Gamma*K(2,2));

QK=@(xu,xv,yu,yv) Gamma*[-2*gamma11*xu*yu-gamma12*(xu*yv+xv*yu); -2*gamma22*xv*yv-gamma21*(xu*yv+xv*yu)];
QD=@(xu,xv,yu,yv) [2*a1*xu*yu+b_c*(xu*yv+xv*yu); 2*a2*xv*yv+b2*(xu*yv+xv*yu)];

QK_rr=QK(rho(1),rho(2),rho(1),rho(2));
QD_rr=QD(rho(1),rho(2),rho(1),rho(2));

w20=(Gamma*K)\(-0.25*QK_rr);
w22=(Gamma*K-4*ks_c*D)\(-0.25*(QK_rr-4*ks_c*QD_rr));

QK_rw20=QK(rho(1),rho(2),w20(1),w20(2));
QD_rw20=QD(rho(1),rho(2),w20(1),w20(2));

QK_rw22=QK(rho(1),rho(2),w22(1),w22(2));
QD_rw22=QD(rho(1),rho(2),w22(1),w22(2));

G_1_3=-(QK_rw20-ks_c*QD_rw20)-0.5*(QK_rw22-ks_c*QD_rw22);

sigma=-c1*ks_c;
L=(G_1_3'*[1; Ms])/(1+M*Ms);
end