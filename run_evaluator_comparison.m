function run_evaluator_comparison(mode)
%RUN_EVALUATOR_COMPARISON Reduced evaluators on a fixed space.
% Same space, same 357 output nodes and same reference for every method.
% Timings include stationary preparation for centered/increment methods.
% The space is the one accepted by Algorithm 5.1 (rksm_bound_controller,
% tau=1e-7, Leja poles, saved in the output folder) and the reference is
% that of run_complete_solves for n=4000, which must have been run first.
if nargin<1, mode="transport"; end
root=fileparts(mfilename('fullpath'));
out=fullfile(root,results_dir,'evaluator_comparison');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,char("run_"+string(mode)+".log")));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
assert(ndre_scaled,'simax:scaled','The evaluator comparison uses the scaled equation.');
a=nla_core; p=a.problem('transport',4000);
loaded=load(fullfile(root,results_dir,'complete_solves','reference_4000.mat'),'ref'); ref=loaded.ref;
spaces="leja_bound_t7";
for space=spaces
    file=fullfile(out,char(space+".mat")); if exist(file,'file'), continue; end
    shifts=pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),319),"leja");
    r=rksm_bound_controller(a,p,1e-12,1e-7,struct('tf',1,'checkDimension',10),shifts);
    require_accepted(r.status,space); save(file,'r','-v7.3');
end
spaceFolder=out;
methods=["centered","incremental-direct","incremental-composed","radon-composed"];
csv='transport_evaluators.csv'; repetitions=timing_reps(6)-1; rows=struct([]);
for space=spaces
    d=load(fullfile(spaceFolder,char(space+".mat")),'r');
    V=d.r.V; W=d.r.W; red=a.project(p,V,W);
    for method=methods
        seconds=zeros(1,repetitions+1);
        for rep=1:repetitions+1
            timer=tic; sol=evaluate(a,red,ref.time,method); seconds(rep)=toc(timer);
        end
        sol.V=V; sol.W=W; e=a.errors(sol,ref);
        budget=e.absoluteError./(1e-12+1e-7*e.referenceNorm);
        [peak,j]=max(e.relativeError,[],'omitmissing');
        row=struct('space',space,'dimension',size(V,2),'method',method, ...
            'maxRelativeError',peak,'peakTime',e.time(j),'maxBudgetRatio',max(budget), ...
            'medianSeconds',timing_pick(seconds),'minSeconds',timing_pick(seconds,"min"), ...
            'maxSeconds',timing_pick(seconds,"max"));
        if isfield(sol,'steps'), row.steps=sol.steps; else, row.steps=NaN; end
        if isfield(sol,'sylvesterScaledResidual'), row.sylvesterResidual=sol.sylvesterScaledResidual; else, row.sylvesterResidual=NaN; end
        rows=[rows;row]; %#ok<AGROW>
        writetable(struct2table(rows),fullfile(out,csv));
        fprintf('EVAL %s m=%d %s err=%.3e@%.2e rho=%.3e time=%.4f\n',space,row.dimension, ...
            method,peak,row.peakTime,row.maxBudgetRatio,row.medianSeconds);
    end
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
