function op=scale_ops(op,scaling)
%SCALE_OPS Scaled operators for li2021_lowrank_exp, in its notation
%   X' = A X + X D + C1t*C2t' - X B1*B2' X,
% under Xhat = diag(dA)*X/diag(dD) (scaling = ps.scaling of scale_problem).
dA=scaling.dA; dD=scaling.dD; o=op;
op.Aop=@(V)dA.*o.Aop(V./dA); op.ATop=@(V)o.ATop(dA.*V)./dA;
op.Dop=@(V)dD.*o.Dop(V./dD); op.DTop=@(V)o.DTop(dD.*V)./dD;
op.C1t=dA.*o.C1t; op.C2t=o.C2t./dD; op.B1=dD.*o.B1; op.B2=o.B2./dA;
end
