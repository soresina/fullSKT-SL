clear 
%close all
clc
% signL - lambda_k fixed, (d12,d21)plane -> d_c _> L
% Parameters (usual parameter set)
% Domain
Lx=1;
% first eigenvalue
lambda_1=(1*pi/Lx)^2;
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

D12=0:0.005:4;

D21_cambio1=zeros(numel(D12),1);
D21_cambio2=zeros(numel(D12),1);
dthreshold=zeros(numel(D12),1);

for i12=1:numel(D12)

% Diffusion part
d12=D12(i12);
d21_thr=(alpha_cross*d12-detK/lambda_1)/abs(beta_cross);
D21=0:0.00001:d21_thr;
LL=zeros(numel(D21),1);

for i21=1:numel(D21)
    d21=D21(i21);
    A1=lambda_1^2;
    B1=d12*vs*lambda_1^2+d21*us*lambda_1^2-trK*lambda_1;
    C1=detK-d12*alpha_cross*lambda_1-d21*beta_cross*lambda_1;
    Delta=B1^2-4*A1*C1;
    if Delta>=0
    d_c=(-B1+sqrt(B1^2-4*A1*C1))/(2*A1);   
    end
    [~,L]=LandauConstant(d_c,d_c,0,0,d21,r1,r2,a1,a2,b1,b2,1,d12,lambda_1);
    LL(i21)=L;
%     if L>0
%         plot(d12,d21,'.','Color',[0.91 0.33 0.5])
%     else
%         plot(d12,d21,'.','Color',[0 0.81 0.82])
%     end
end
if d21_thr>0
   change=find(diff(sign(LL)))+1;
    if numel(change)==1    
    D21_cambio2(i12)=D21(change);
    elseif numel(change)>1
    D21_cambio1(i12)=min(D21(change));
    D21_cambio2(i12)=max(D21(change));
    end
dthreshold(i12)=D21(end);
end
end
%%
% figure()
% hold on
% box on
% area(D12, dthreshold,0,'EdgeColor','none','FaceColor',[0.91 0.33 0.5],'FaceAlpha',.3)
% area(D12, D21_cambio2,0,'EdgeColor','none','FaceColor',[1 1 1])
% area(D12, D21_cambio2,0,'EdgeColor','none','FaceColor',[0 0.81 0.82],'FaceAlpha',.3)
% area(D12, D21_cambio1,0,'EdgeColor','none','FaceColor',[1 1 1])
% area(D12, D21_cambio1,0,'EdgeColor','none','FaceColor',[0.91 0.33 0.5],'FaceAlpha',.3)
plot(D12,D21_cambio2,'b')
plot(D12,dthreshold,'b')
plot(D12,D21_cambio1,'b')
axis([0 4 0 0.08])