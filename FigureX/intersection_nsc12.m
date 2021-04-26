function [ds,d12s]=intersection_nsc12(d21,r1,r2,a1,a2,b1,b2,Lx)

syms d d12 lambda_k

% Lx=1;
% r1=5;
% r2=2;
% a1=3;
% a2=3;
% b1=1;
% b2=1;
% 
% d21=0;

lambda_1=(pi/Lx)^2;
lambda_2=(2*pi/Lx)^2;

wos=a1*a2-b1*b2;
us=(r1*a2-r2*b1)/wos;
vs=(r2*a2-r1*b2)/wos;

% Curves
NSC=(d.^2*lambda_k^2+d21*us*lambda_k^2*d+d21*lambda_k*(a1*us-b1*vs)*us+(a1*us+a2*vs)*lambda_k*d+wos*us*vs)./(-d*vs*lambda_k^2+(b2*us-a2*vs)*vs*lambda_k);

curve_1=subs(NSC,lambda_k,lambda_1);
curve_2=subs(NSC,lambda_k,lambda_2);

ds=solve(curve_1-curve_2);
%disp(ds)
if isempty(ds)
    ds=0;
    d12s=0;
else
Ds=[double(ds(1)), double(ds(2))];
if Ds(1)>0
    ds=Ds(1);
else
    ds=Ds(2);
end
d12s=double(subs(curve_1,ds));
end
%disp([ds, d12s])
end