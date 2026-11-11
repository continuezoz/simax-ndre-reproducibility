function r=rksm_bound_controller(a,p,tauA,tau,o,shifts)
%RKSM_BOUND_CONTROLLER Algorithm 1: two-sided rational Krylov expansion of
% the scaled problem p (see ndre_setup), stopped by the a posteriori error
% estimate. The reduced equation is solved around its contractive
% stationary solution and the estimate epsilon_j of error_estimate is
% formed on geometric nodes with ratio at most 1.1 from 0.01/max|shift| to
% o.tf; the space is accepted when
%   epsilon_j <= tauA + tau*||Y(sigma_j)||_F   for all j.
% Candidate spaces are tested when the left dimension reaches a multiple
% of o.checkDimension (every block if o.checkDimension is 0). Dependent
% directions of a new block are deflated by pivoted QR. X0 = 0.
% r.basisSeconds and r.checkSeconds split the time between basis expansion
% and the tests of candidate spaces.
timer=tic; nu=p.nu; lo=.01/max(abs(shifts));
q=ceil(log(o.tf/lo)/log(1.1)); nodes=[0,logspace(log10(lo),log10(o.tf),q+1)];
[V,nv]=append(zeros(p.n,0),p.F); [W,nw]=append(zeros(p.p,0),p.G);
trace=table('Size',[0 6],'VariableTypes',repmat({'double'},1,6), ...
    'VariableNames',{'k','ell','ratio','margin','norm2Ystar','seconds'});
accepted=false; status="MAX_DIMENSION"; est=[]; next=o.checkDimension; tb=0; tc=0; poles=0;
for shift=shifts
    clock=tic; poles=poles+1;
    if nv>0, [V,nv]=append(V,p.Ashift(V(:,end-nv+1:end),shift)); end
    if nw>0, [W,nw]=append(W,p.DTshift(W(:,end-nw+1:end),shift)); end
    k=size(V,2); tb=tb+toc(clock);
    if o.checkDimension>0
        if k<next && nv+nw>0, continue; end
        next=next+o.checkDimension;
    end
    clock=tic; red=a.project(p,V,W); prep=a.prepare(red,false);
    [~,RA]=qr(p.Aapply(V)-V*red.A,0);
    [~,RD]=qr(p.DTapply(W)-W*red.D',0);
    est=error_estimate(red,prep,nodes,RA,RD,nu);
    ratio=max(est.epsilon./(tauA+tau*est.normY)); tc=tc+toc(clock);
    trace(end+1,:)={k,size(W,2),ratio,prep.minClosedLoopRealPart,norm(prep.Ystar),toc(timer)}; %#ok<AGROW>
    if ratio<=1 && ~(isfield(o,'noStop') && o.noStop), accepted=true; status="OK_ERROR_ESTIMATE"; break; end
    if nv+nw==0, status="INVARIANT"; break; end
end
r=struct('V',V,'W',W,'accepted',accepted,'status',status,'trace',trace, ...
    'estimate',est,'nu',nu,'nodes',nodes,'seconds',toc(timer), ...
    'basisSeconds',tb,'checkSeconds',tc,'poles',poles);
end

function [Q,added]=append(V,Z)
% Two-pass orthogonalization against V and pivoted QR with deflation.
scale=norm(Z,'fro'); Z=Z-V*(V'*Z); Z=Z-V*(V'*Z);
[q,R,~]=qr(Z,0); d=abs(diag(R)); added=sum(d>1e-10*max(scale,eps));
q=q(:,1:added); q=q-V*(V'*q); [q,~]=qr(q,0);
Q=[V,q];
end
