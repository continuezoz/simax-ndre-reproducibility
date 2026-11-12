function ps=ndre_setup(p)
%NDRE_SETUP Scaled problem of section 2.1 and its constant nu, (3.2).
% Returns the problem structure of the scaled equation (scale_problem) with
% the additional field nu. Supported families:
%   "transport"   explicit scaling, nu from the secular equation;
%   "multistate"  Lemma 2.1 with u = e and v = K\e, where
%                 K = diag([k^-;k^+]) - U*V' has rank 4r; nu from the pencil
%                 (G,J), G = Khat + Khat', by eigs with solves with G.
% A problem that is already scaled is returned unchanged.
if isfield(p,'scaling'), ps=p; return; end
switch string(p.family)
    case "transport"
        c=0.9999; alpha=1e-8;
        if numel(p.parameter)==2, c=p.parameter(1); alpha=p.parameter(2); end
        [piA,piD,nu]=transport_scaling(p.n,c,alpha);
        ps=scale_problem(p,piA,piD); ps.nu=nu;
    case "multistate"
        f=p.factors; n=p.n; r=size(f.F1,2); Z=zeros(n,r);
        dg=[f.km;f.kp];
        U=[f.F1,f.B1,Z,Z;Z,Z,f.B1,f.F1];
        V=[f.km.*f.F2,Z,f.km.*f.B2,Z;Z,f.kp.*f.B2,Z,f.kp.*f.F2];
        v=lowrank_solve(dg,U,V,ones(2*n,1));
        assert(all(v>0),'setup:scaling','Scaling vector is not positive.');
        piD=1./v(1:n); piA=1./v(n+1:end);
        ps=scale_problem(p,piA,piD);
        d=sqrt([piD;piA]); Uh=d.*U; Vh=V./d;      % Khat = diag(dg) - Uh*Vh'
        Gsolve=@(x)lowrank_solve(2*dg,[Uh,Vh],[Vh,Uh],x);
        ps.nu=pencil_nu(Gsolve,n,n);
    otherwise
        error('setup:family','No scaling implemented for %s.',p.family);
end
end

function x=lowrank_solve(d,U,V,b)
% (diag(d) - U*V')\b by the Sherman--Morrison--Woodbury formula.
z=b./d; DU=U./d; x=z+DU*((eye(size(U,2))-V'*DU)\(V'*z));
end

function nu=pencil_nu(Gsolve,p,n)
% 2*nu is the smallest modulus of the eigenvalues of G x = lambda J x.
J=[ones(p,1);-ones(n,1)];
opts=struct('tol',1e-12,'maxit',500,'issym',false);
theta=eigs(@(x)Gsolve(J.*x),p+n,2,'largestabs',opts);
nu=1/(2*max(abs(theta)));
end
