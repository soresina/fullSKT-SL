# fullSKT-SL-H

Related to the paper "Hopf bifurcations in the full SKT model and where to find them" by C. Soresina

Each folder corresponds to a figure of the paper and contains the Matlab files to generate it. In detail:
- Figure 1: Neutral stability curves in the (d,d12)-plane for a fixed value of d21.
- Figure 5: Neutral stability curves with the sign of L, in the (d,d21)-plane for a fixed value of d12.
- Figure 6: Sign of L in the (d12,d21)-plane for a fixed value of d12.
- Figure 7: Neutral stability curves with the sign of L, in the (d,d12)-plane for a fixed value of d21.
- Figure X: Coefficient A1 and B1.

To obtain the bifurcation diagrams (Fig.2), we refer to the GitHub folder fullSKT
https://github.com/soresina/fullSKT

To compute the branches starting from a Hopf point, use:

%---------------------------------------------------------------------------------

%% Hopf bifs from HBPs
para=4; ds=0.1; dsmax=0.5; xi=1e-1; nsteps=50; figure(2); clf; 
aux=[]; aux.tl=30; 
p=hoswibra('bpt1_up','hpt1',ds,para,'1dh1',aux); 
p.hopf.xi=xi; p.hopf.jac=1; p.nc.dsmax=dsmax; p.sw.verb=2;  p.file.smod=1; 
p=setbel(p,2,1e-4,5,@lss); p.hopf.fltol=1e-3; 
t1=tic; p=cont(p,nsteps); toc(t1) 

%---------------------------------------------------------------------------------
