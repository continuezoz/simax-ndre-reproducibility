function sol=li2021_lowrank_exp(op,times,h,scheme,opts)
%LI2021_LOWRANK_EXP Low-rank exponential Rosenbrock integrators of Li,
% Zhang and Liu (J. Comput. Appl. Math. 389, 2021) for the NDRE
%   X' = A X + X D + Q - X G X,  Q = C1t*C2t',  G = B1*B2',
% written in their notation (Algorithm 5 = "expeuler", Algorithm 6 = "erow3").
% X_n = L_n*S_n*U_n'. phi_k(h S_n) is applied by Gauss-Legendre quadrature,
% eq. (28), with exponential actions computed by scaling and Taylor
% truncation; factors are compressed by their Algorithm 4 (truncated SVD
% of the core) with tolerance N*eps. The second-stage difference uses
% G_n(X2)-G_n(Xn) = -(X2-Xn) G (X2-Xn), which equals their (54).
% op fields: Aop(V), ATop(V), Dop(V), DTop(V), C1t, C2t, B1, B2, normA, normD.
if nargin<5, opts=struct; end
if ~isfield(opts,'nodes'), opts.nodes=6; end
if ~isfield(opts,'taylorTol'), opts.taylorTol=1e-14; end
[M,~]=size(op.C1t); N=size(op.C2t,1);
ctol=max(M,N)*eps;
[gx,gw]=gauss_legendre01(opts.nodes);
L=zeros(M,0); S=zeros(0); U=zeros(N,0);
steps=round(times(end)/h); assert(abs(steps*h-times(end))<1e-12*max(1,times(end)));
sample=round(times/h); assert(max(abs(sample*h-times))<1e-12,'Output nodes must be steps.');
store=cell(1,numel(times)); store{1}={L,S,U}; next=2; ranks=zeros(1,steps);
for n=1:steps
    % Linearization A_n = A - X_n G and D_n = D - G X_n, applied matrix-free.
    % G = B1*B2' is N x M, so X_n G X_n = L*(W2*W1)*U' with the cores below.
    W1=(op.B2'*L)*S;  W2=S*(U'*op.B1);                       % q x r, r x q
    An=@(V) op.Aop(V)-L*(W2*(op.B2'*V));                      % (A - X_n G) V
    DnT=@(V) op.DTop(V)-U*(S'*(L'*(op.B2*(op.B1'*V))));       % (D - G X_n)' V
    nA=op.normA+norm(L,1)*norm(W2*op.B2',1);
    nD=op.normD+norm(U*S',1)*norm(L'*op.B2*op.B1',1);
    % F(X_n) = Lt*St*Ut' (their (51)).
    r=size(S,1); l=size(op.C1t,2);
    Lt=[op.C1t,op.Aop(L),L]; Ut=[op.C2t,op.DTop(U),U];
    St=blkdiag(eye(l),zeros(2*r)); St(l+1:l+r,l+r+1:end)=S;
    St(l+r+1:end,l+1:l+r)=S; St(l+r+1:end,l+r+1:end)=-W2*W1;
    [Lt,St,Ut]=compress(Lt,St,Ut,ctol);
    [Y,T,Z]=phi_lowrank(An,DnT,nA,nD,h,Lt,St,Ut,1,gx,gw,opts.taylorTol);
    [L2,S2,U2]=compress([L,Y],blkdiag(S,h*T),[U,Z],ctol);
    if scheme=="erow3"
        % Delta = X_{n,2}-X_n = h*Y*T*Z'; G_n(X2)-G_n(Xn) = -Delta*G*Delta.
        Dl=Y; Ds=h*T; Du=Z;
        Sb=-Ds*((Du'*op.B1)*(op.B2'*Dl))*Ds;
        [Lb,Sb,Ub]=compress(Dl,Sb,Du,ctol);
        [Y3,T3,Z3]=phi_lowrank(An,DnT,nA,nD,h,Lb,Sb,Ub,3,gx,gw,opts.taylorTol);
        [L,S,U]=compress([L2,Y3],blkdiag(S2,2*h*T3),[U2,Z3],ctol);
    else
        L=L2; S=S2; U=U2;
    end
    ranks(n)=size(S,1);
    while next<=numel(times) && sample(next)==n, store{next}={L,S,U}; next=next+1; end
end
sol=struct('factors',{store},'time',times,'ranks',ranks,'steps',steps);
end

function [Y,T,Z]=phi_lowrank(An,DnT,nA,nD,h,L,S,U,k,gx,gw,tol)
% phi_k(h S_n)(L S U') ~ sum_j w_j s_j^(k-1)/(k-1)! e^{(1-s_j)h A_n} L S
% (e^{(1-s_j)h D_n'} U)'. Nodes are processed by increasing (1-s_j) so each
% exponential action starts from the previous one.
[tau,order]=sort((1-gx)*h); w=gw(order).*gx(order).^(k-1)/factorial(k-1);
q=numel(tau); r=size(S,1); Y=zeros(size(L,1),q*r); Z=zeros(size(U,1),q*r);
EL=L; EU=U; prev=0;
for j=1:q
    EL=expm_action(An,nA,tau(j)-prev,EL,tol); EU=expm_action(DnT,nD,tau(j)-prev,EU,tol);
    prev=tau(j); Y(:,(j-1)*r+(1:r))=EL; Z(:,(j-1)*r+(1:r))=EU;
end
T=kron(diag(w),S);
end

function V=expm_action(Aop,normA,t,V,tol)
% exp(t A) V by scaling with s substeps of Taylor series truncated adaptively.
if t==0, return; end
s=max(1,ceil(t*normA)); dt=t/s;
for i=1:s
    term=V; F=V;
    for m=1:60
        term=(dt/m)*Aop(term); F=F+term;
        if norm(term,1)<=tol*norm(F,1), break; end
    end
    V=F;
end
end

function [L,S,U]=compress(L,S,U,tol)
% Their Algorithm 4: QR of both factors, truncated SVD of the core.
if isempty(S), return; end
[Q1,R1]=qr(L,0); [Q2,R2]=qr(U,0); [P,Sig,Qv]=svd(R1*S*R2');
sv=diag(Sig); k=sum(sv>tol*max(sv(1),realmin));
L=Q1*P(:,1:k); S=Sig(1:k,1:k); U=Q2*Qv(:,1:k);
end

function [x,w]=gauss_legendre01(m)
% Golub-Welsch nodes and weights on [0,1].
b=(1:m-1)./sqrt(4*(1:m-1).^2-1); [V,Dg]=eig(diag(b,1)+diag(b,-1));
x=(diag(Dg)+1)/2; w=V(1,:)'.^2; [x,i]=sort(x); w=w(i);
end
