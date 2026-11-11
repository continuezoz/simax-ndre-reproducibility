function ps=scale_problem(p,piA,piD)
%SCALE_PROBLEM Diagonal scaling of Section 2.1: Xhat = PiA^(1/2) X PiD^(-1/2).
% piA (n-by-1) and piD (p-by-1) are the positive diagonals of Pi_A and Pi_D.
% The returned problem has the operators of the scaled equation
%   Ahat = PiA^(1/2) A PiA^(-1/2),  Dhat = PiD^(1/2) D PiD^(-1/2),
%   Shat = PiD^(1/2) S PiA^(-1/2),  Qhat = PiA^(1/2) Q PiD^(-1/2),
% so every solver of this package runs on it unchanged. A solution of the
% original equation is X = diag(dA)\Xhat*diag(dD) with ps.scaling.
dA=sqrt(piA(:)); dD=sqrt(piD(:));
ps=p; ps.F=dA.*p.F; ps.G=p.G./dD;
ps.normQ=sqrt(sum(sum((ps.F'*ps.F).*(ps.G'*ps.G))));
ps.Aapply=@(X)dA.*p.Aapply(X./dA);
ps.Dapply=@(X)dD.*p.Dapply(X./dD);
ps.DTapply=@(X)p.DTapply(dD.*X)./dD;
ps.Sapply=@(X)dD.*p.Sapply(X./dA);
ps.Asolve=@(X)dA.*p.Asolve(X./dA);
ps.Dsolve=@(X)dD.*p.Dsolve(X./dD);
ps.DTsolve=@(X)p.DTsolve(dD.*X)./dD;
ps.Ashift=@(X,s)dA.*p.Ashift(X./dA,s);
ps.DTshift=@(X,s)p.DTshift(dD.*X,s)./dD;
ps.rightBasisApply=ps.DTapply; ps.rightBasisSolve=ps.DTsolve;
ps.scaling=struct('dA',dA,'dD',dD);
end
