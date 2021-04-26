clear 
close all
clc
% signL - lambda_k fixed, (d12,d21)plane -> d_c _> L
% Parameters (usual parameter set)
% Domain
Lx=1;

% Reaction part
r1=5;
r2=2;
a1=3;
a2=3;
b1=1;
b2=1;

% homogeneous equilibrium state
wos=a1*a2-b1*b2;
us=(r1*a2-r2*b1)/wos;
vs=(r2*a2-r1*b2)/wos;

% Linearization of the reaction part at (us,vs)
K=[-a1*us, -b1*us; -b2*vs, -a2*vs];
detK=det(K);
trK=trace(K);

alpha_cross=(b2*us-a2*vs)*vs;
beta_cross=(b1*vs-a1*us)*us;


d21=0.025;
D=0:0.0001:0.04;
NSC=@(d,lambda_k) (d.^2*lambda_k^2+d21*us*lambda_k^2*d+d21*lambda_k*(a1*us-b1*vs)*us+(a1*us+a2*vs)*lambda_k*d+wos*us*vs)./(-d*vs*lambda_k^2+(b2*us-a2*vs)*vs*lambda_k);

figure()
hold on
box on

for il=6:-1:1
    % first eigenvalue
lambda_1=(il*pi/Lx)^2;
% Curves
D12=0:0.01:4;
LL=zeros(numel(D12),1);

for i21=1:numel(D12)
    d12=D12(i21);
    A1=lambda_1^2;
    B1=d12*vs*lambda_1^2+d21*us*lambda_1^2-trK*lambda_1;
    C1=detK-d12*alpha_cross*lambda_1-d21*beta_cross*lambda_1;
    Delta=B1^2-4*A1*C1;
    if Delta>=0
    d_c=(-B1+sqrt(B1^2-4*A1*C1))/(2*A1);   
    
    [~,L]=LandauConstant(d_c,d_c,0,0,d21,r1,r2,a1,a2,b1,b2,1,d12,lambda_1);
    LL(i21)=L;
    if il<3
    if L>0
        plot(d_c,d12,'.','Color',[0.91 0.33 0.5])
    else
        plot(d_c,d12,'.','Color',[0 0.81 0.82])
    end
    else 
        plot(d_c,d12,'.','Color',[0.86 0.86 0.86])
    end
    end
end
end

% d12_ref=1.7;
% plot([0 0.04],[d12_ref d12_ref],':k')
% d12_ref=2.5;
% plot([0 0.04],[d12_ref d12_ref],':k')
% d12_ref=3;
% plot([0 0.04],[d12_ref d12_ref],':k')

set(gcf,'color','w');
axis([0 0.04 0 4])
ax = gca;
ax.FontSize = 16; 