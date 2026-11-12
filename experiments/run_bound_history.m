function run_bound_history(n,tau)
%RUN_BOUND_HISTORY Figure 1 data: node-wise true error and error estimate
% of Algorithm 5.1 for the near-critical single-state model and the
% multistate model with gamma=0.9999; dense full-order reference as in
% run_bound_experiments("dense"). Output: data/bounds/history_*.csv.
if nargin<1, n=600; end
if nargin<2, tau=1e-7; end
out=fullfile(ndre_paths().data,'bounds'); if ~exist(out,'dir'), mkdir(out); end
previous=ndre_scaled; restore=onCleanup(@()ndre_scaled(previous)); %#ok<NASGU>
ndre_scaled(false); a=ndre_core;
p0=a.problem("transport",n,n,[],[0.9999 1e-8]);
cases(1)=struct('id',"transport",'p',ndre_setup(p0),'lo',p0.shiftRange(1), ...
    'hi',p0.shiftRange(2)*(1+1e-8)/(1-1e-8),'block',10);
p0=multistate_problem(n,10,0.9999,0.1);
cases(2)=struct('id',"multistate",'p',ndre_setup(p0),'lo',min(p0.factors.kp), ...
    'hi',max(p0.factors.km),'block',0);
for c=cases
    p=c.p; I=eye(p.n); full=a.project(p,I,I); fp=a.prepare(full,false);
    shifts=pole_order(-logspace(log10(c.lo),log10(c.hi),159),"leja");
    r=rksm_bound_controller(a,p,1e-12,tau,struct('tf',1,'checkDimension',c.block),shifts);
    e=r.estimate; q=numel(r.nodes);
    ref=incremental_centered_nodes(full,fp,r.nodes(2:end),"composed");
    % error: computed output (increment evaluation) against the dense
    % reference; projectionError: reduced flow evaluated independently by
    % scaled Radon composition, so that the evaluation error is excluded.
    red=a.project(p,r.V,r.W); err=zeros(q,1); refn=err; perr=err;
    for j=2:q
        X=ref.Yhistory(:,:,j-1); refn(j)=norm(X,'fro');
        err(j)=norm(X-r.V*e.Yhistory(:,:,j)*r.W','fro');
        mr=a.mobius(red,[0 r.nodes(j)],1,0); perr(j)=norm(X-r.V*mr.Y*r.W','fro');
    end
    T=table(r.nodes(:),err,perr,e.epsilon(:),refn,'VariableNames',{'t','error','projectionError','estimate','referenceNorm'});
    writetable(T,fullfile(out,char("history_"+c.id+".csv")));
    fprintf('HISTORY %s k=%d nodes=%d max err/est (est>1e-11)=%.3f\n',c.id,size(r.V,2),q, ...
        max(err(e.epsilon>1e-11)./e.epsilon(e.epsilon>1e-11)));
end
end
