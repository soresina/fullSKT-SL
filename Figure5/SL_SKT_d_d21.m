clear 
close all
clc
% Figure - sign of L on the neutral stability curves (d,d21)-plane

% Parameters (usual parameter set)
% Domain
Lx=1;
% !!! Please note !!!
% Here u and v are switched wrt the paper
% Reaction part
r1=2;
r2=5;
a1=3;
a2=3;
b1=1;
b2=1;
% Diffusion part
d21=3; 

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

% Curves
NSC=@(d,lambda_k) (d.^2*lambda_k^2+d21*us*lambda_k^2*d+d21*lambda_k*(a1*us-b1*vs)*us+(a1*us+a2*vs)*lambda_k*d+wos*us*vs)./(-d*vs*lambda_k^2+(b2*us-a2*vs)*vs*lambda_k);
D=0:0.0001:0.04;
n=numel(D);


figure()
hold on
box on

for j=1:2
    % first eigenvalue
    lambda_1=(j*pi/Lx)^2;
    D12=NSC(D,lambda_1);

LL=zeros(n,1);

for i=1:n
    d=D(i);
    d12=D12(i);
    [~,L]=LandauConstant(d,d,0,0,d21,r1,r2,a1,a2,b1,b2,1,d12,lambda_1);
    LL(i)=L;
    %SS(i)=S;
    if L>0
        plot(d,d12,'.','Color',[0.91 0.33 0.5])
    else
        plot(d,d12,'.','Color',[0 0.81 0.82])
    end
end


end
set(gcf,'color','w');
axis([0 0.04 0 0.06])
ax = gca;
ax.FontSize = 16; 
xlabel('d')
ylabel('d_{21}')
