clear 
close all
clc

%% Parameters
% Domain
Lx=1;
% Reaction part
r1=5;
r2=2;
a1=3;
a2=3;
b1=1;
b2=1;
% Diffusion part
d21=0.012;
d=0:0.0001:0.04;

% homogeneous equilibrium state
wos=a1*a2-b1*b2;
us=(r1*a2-r2*b1)/wos;
vs=(r2*a2-r1*b2)/wos;

%% Neutral Stability Curves
NSC=@(d,lambda_k) (d.^2*lambda_k^2+d21*us*lambda_k^2*d+d21*lambda_k*(a1*us-b1*vs)*us+(a1*us+a2*vs)*lambda_k*d+wos*us*vs)./(-d*vs*lambda_k^2+(b2*us-a2*vs)*vs*lambda_k);

figure(2)
hold on
box on

% Grey region
for i=1:5
    lambda_k=(pi*i)^2;
    curve_k=NSC(d,lambda_k);
    if i==1 && d21==0
        ick=numel(curve_k);
        patch([d d(end) 0], [curve_k 4 4], [0.86 0.86 0.86])
    else
        ick=find(curve_k>4,1);
        patch([d(1:ick) 0], [curve_k(1:ick) 4], [0.86 0.86 0.86])
    end
end

% plot the curves
for i=1:5
    lambda_k=(pi*i)^2;
    curve_k=NSC(d,lambda_k);
    if i==1 && d21==0
        ick=numel(curve_k);
        plot(d(1:ick),curve_k(1:ick),'k')
    else
        ick=find(curve_k>4,1);
        plot(d(1:ick),curve_k(1:ick),'k')
    end
    
end
% some d_12 of interest
plot([0 0.04],[3 3],':','Linewidth',2,'Color',[0 0.8 0.6])
plot([0 0.04],[4 4],'k','Linewidth',0.1)
% to obtain a nice figure
ax = gca;
ax.FontSize = 16; 
axis([0 0.04 0 4])
xlabel('d')  
ylabel('d_{12}')
