function [piA,piD,nu,nubar]=transport_scaling(n,c,alpha)
%TRANSPORT_SCALING Explicit scaling of the single-state model and its
% constants nu and nubar. Pi_D = diag(q)^(-1), Pi_A = diag(q) make the coefficient
% matrix symmetric, diag(Gamma,Delta) - [s;s][s;s]' with s = sqrt(q), and
% nu is the smaller modulus of the two roots of the secular equation
%   sum q_i/(Gamma_i - lambda) + sum q_i/(Delta_i + lambda) = 1
% next to the origin; nubar, (3.4), is half the distance between these roots.
[omega,w]=load_gauss_legendre_matlab(n); q=w./(2*omega);
Delta=1./(c*omega*(1+alpha)); Gamma=1./(c*omega*(1-alpha));
piA=q; piD=1./q;
f=@(lambda) sum(q./(Gamma-lambda))+sum(q./(Delta+lambda))-1;
rp=fzero(f,[0,min(Gamma)*(1-1e-13)]);
rn=fzero(f,[-min(Delta)*(1-1e-13),0]);
nu=min(rp,-rn); nubar=(rp-rn)/2;
end
