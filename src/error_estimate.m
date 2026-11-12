function est=error_estimate(red,prep,nodes,RA,RD,nu)
%ERROR_ESTIMATE A posteriori error estimate (4.5) of the paper.
% nodes are increasing with nodes(1)=0. At each node the weighted residual
%   omega^2 = ||a*Delta*b||_F^2 + ||RA*Y*b||_F^2 + ||a*Y*RD'||_F^2,
%   a = (I-Y*Y')^(-1/2),  b = (I-Y'*Y)^(-1/2),
% is formed from the reduced state Y (composed increment evaluation), with
% Delta = -N(Ystar) the reduced defect, and RA, RD the triangular QR
% factors of the projection residuals. The bound
%   int_0^t exp(-2*nu*(t-s))*omega(s) ds
% is evaluated with the exponential weight carried exactly between nodes and
% the larger of the two weighted endpoint values of omega on each
% subinterval (an upper bound when the weighted integrand is monotone there).
% Absolute quantities (scaled Frobenius norm).
assert(nodes(1)==0 && all(diff(nodes)>0),'estimate:nodes','Nodes must start at 0 and increase.');
sol=incremental_centered_nodes(red,prep,nodes,"composed");
q=numel(nodes); omega=zeros(q,1); epsilon=zeros(q,1); normY=zeros(q,1); norm2Y=zeros(q,1);
N=-red.A*prep.Ystar-prep.Ystar*red.D+prep.Ystar*red.S*prep.Ystar+red.Q;
for j=1:q
    Y=sol.Yhistory(:,:,j); [U,sv,Z]=svd(Y,'econ'); sv=diag(sv);
    norm2Y(j)=max([sv;0]); normY(j)=norm(Y,'fro');
    assert(norm2Y(j)<1,'estimate:ball','Reduced state left the unit ball.');
    g=1./sqrt(1-sv.^2)-1;
    left=@(M) M+U*(g.*(U'*M)); right=@(M) M+((M*Z).*g')*Z';
    omega(j)=sqrt(norm(right(left(N)),'fro')^2+norm(right(RA*Y),'fro')^2 ...
        +norm(left(Y*RD'),'fro')^2);
    if j>1
        h=nodes(j)-nodes(j-1); w=exp(-2*nu*h);
        epsilon(j)=w*epsilon(j-1)+h*max(w*omega(j-1),omega(j));
    end
end
est=struct('nodes',nodes(:),'omega',omega,'epsilon',epsilon,'normY',normY, ...
    'norm2Y',norm2Y,'Yhistory',sol.Yhistory);
end
