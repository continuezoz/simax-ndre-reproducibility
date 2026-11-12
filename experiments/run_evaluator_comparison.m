function run_evaluator_comparison
%RUN_EVALUATOR_COMPARISON Table 5: evaluators of the reduced flow on one
% fixed space, the space accepted by Algorithm 5.1 for the single-state
% model with n=4000 and tau=1e-7. All evaluators use the same 121 output
% nodes and the reference of run_complete_solves for n=4000, which must be
% present. Timings include the stationary preparation.
% Output: data/evaluator_comparison/transport_evaluators.csv.
out=fullfile(ndre_paths().data,'evaluator_comparison');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'evaluator_comparison.log'));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
assert(ndre_scaled,'ndre:scaled','The evaluator comparison uses the scaled equation.');
a=ndre_core; p=a.problem('transport',4000);
loaded=load(fullfile(ndre_paths().data,'complete_solves','reference_4000.mat'),'ref'); ref=loaded.ref;
space="leja_bound_t7"; file=fullfile(out,char(space+".mat"));
if ~exist(file,'file')
    shifts=pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),319),"leja");
    r=rksm_bound_controller(a,p,1e-12,1e-7,struct('tf',1,'checkDimension',10),shifts);
    require_accepted(r.status,space); save(file,'r','-v7.3');
end
d=load(file,'r'); V=d.r.V; W=d.r.W; red=a.project(p,V,W);
methods=["centered","incremental-direct","incremental-composed","radon-composed"];
rows=struct([]);
for method=methods
    timer=tic; sol=evaluate(a,red,ref.time,method); seconds=toc(timer);
    sol.V=V; sol.W=W; e=a.errors(sol,ref);
    budget=e.absoluteError./(1e-12+1e-7*e.referenceNorm);
    [peak,j]=max(e.relativeError,[],'omitmissing');
    row=struct('space',space,'dimension',size(V,2),'method',method, ...
        'maxRelativeError',peak,'peakTime',e.time(j),'maxBudgetRatio',max(budget), ...
        'seconds',seconds);
    if isfield(sol,'sylvesterScaledResidual'), row.sylvesterResidual=sol.sylvesterScaledResidual; else, row.sylvesterResidual=NaN; end
    rows=[rows;row]; %#ok<AGROW>
    writetable(struct2table(rows),fullfile(out,'transport_evaluators.csv'));
    fprintf('EVAL %s m=%d %s err=%.3e@%.2e rho=%.3e time=%.4f\n',space,row.dimension, ...
        method,peak,row.peakTime,row.maxBudgetRatio,seconds);
end
end

function sol=evaluate(a,red,times,method)
if method=="radon-composed", sol=scaled_radon_nodes(red,times,"composed"); return; end
% Only the direct centered formula needs the Sylvester kernel L.
prep=a.prepare(red,~startsWith(method,"incremental"));
if startsWith(method,"incremental")
    sol=incremental_centered_nodes(red,prep,times,extractAfter(method,"incremental-")); return
end
[k,ell]=size(red.Q); history=zeros(k,ell,numel(times));
for j=1:numel(times)
    if times(j)>0, value=a.evaluate(prep,times(j),false); history(:,:,j)=value.Y; end
end
sol=struct('Yhistory',history,'Y',history(:,:,end),'time',times, ...
    'sylvesterScaledResidual',prep.sylvesterScaledResidual);
end
