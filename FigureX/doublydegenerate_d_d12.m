clear
close all
clc

% Parameters (usual parameter set)
% Domain
Lx=1;
% first eigenvalue
lambda_1=(pi/Lx)^2;
% Reaction part
r1=5;
r2=2;
a1=3;
a2=3;
b1=1;
b2=1;
% Diffusion part
D21=0:0.001:0.05;

Ds=zeros(size(D21));
D12s=zeros(size(D21));
nc_ddp=zeros(size(D21));
A1=zeros(size(D21));
B1=zeros(size(D21));

% doubly degenerate point varying d21
for i=1:numel(D21)
    d21=D21(i);
    [Ds(i),D12s(i)]=intersection_nsc12(d21,r1,r2,a1,a2,b1,b2,Lx);
end

% necessary condition at ddp
for i=1:numel(D21)
    d21=D21(i);
    ds=Ds(i);
    d12s=D12s(i);
    % paper
    [nc_ddp(i),A1(i),B1(i)]=nec_cond_at_ddp(ds,d12s,d21,r1,r2,a1,a2,b1,b2,Lx);
    %scambio u,v
    %[nc_ddp(i),A1(i),B1(i)]=nec_cond_at_ddp_nuova(ds,d12s,d21,r1,r2,a1,a2,b1,b2,Lx);
end

A1_loro=13*(324175+62301*sqrt(105))/150880;
B1_loro=13*(4323445-424489*sqrt(105))/104960;

% disp([A1 A1_loro])
% disp([B1 B1_loro])

%%
figure()
hold on
box on
for i=1:numel(D21)
    if nc_ddp(i)>0
        plot(Ds(i),D12s(i),'.r')
    else
        plot(Ds(i),D12s(i),'.b')
    end
end

%%
figure()
hold on
box on
for i=1:numel(D21)
    if nc_ddp(i)>0
        plot(D21(i),D12s(i),'.r')
    else
        plot(D21(i),D12s(i),'.b')
    end
end
%%
figure()
plot(D21,A1,D21,B1)
legend('A1','B1')
