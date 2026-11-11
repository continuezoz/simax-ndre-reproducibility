function run_complete_solves(mode)
%RUN_COMPLETE_SOLVES E6: complete adaptive solves on the near-critical
% transport benchmark for n = 4000, 10000, 20000, 40000.
% Timed work: basis construction, all candidate checks, projection, and the
% evaluation of all output states. References and error checks are untimed.
% The projected solver uses the transient check on geometric guard nodes
% between 0.01/lambda_max and t_f. Mode "rksm" reruns only the projected
% solver rows and replaces them in complete_solves.csv.
if nargin<1, mode="run"; end
root=fileparts(mfilename('fullpath'));
out=fullfile(root,results_dir,'complete_solves');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,char("run_"+string(mode)+".log")));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
a=nla_core; times=unique([0,logspace(-8,-2,21),0:.01:1]); target=1e-7;
rows=struct([]); csv=fullfile(out,'complete_solves.csv');
for n=[4000 10000 20000 40000]
    p=a.problem('transport',n);
    file=fullfile(out,sprintf('reference_%d.mat',n));
    if exist(file,'file'), s=load(file,'ref'); ref=s.ref;
    else
        % For n=40000 dimension 240 is not fine enough (run_reference_40000);
        % the scaled run uses 360 and 440 directly.
        dims=[240 280]; if ndre_scaled && n==40000, dims=[360 440]; end
        ref=eba_reference(a,p,dims(1),times); fine=eba_reference(a,p,dims(2),times);
        c=a.errors(ref,fine); ref.spaceRefinement=max(c.absoluteError./(1e-12+1e-7*c.referenceNorm));
        save(file,'ref','-v7.3');
    end
    fprintf('REFERENCE n=%d refinement=%.3e\n',n,ref.spaceRefinement);
    o=struct('tf',1,'checkDimension',10);
    shifts=pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),319),"leja");
    guard=logspace(log10(.01/p.shiftRange(2)),log10(o.tf),49);
    % 2026-11-11: for n=40000 the computed estimate has a floor at about
    % 2.35*1e-7 (rounding in A*V, ||A|| ~ 1.1e9), so tau=1e-6 there; the
    % rejected tau=1e-7 search is kept as the record of that floor.
    tau=target;
    if n>=40000
        tau=1e-6;
        floorRun=rksm_bound_controller(a,p,1e-12,target,o,shifts);
        writetable(floorRun.trace,fullfile(out,sprintf('floor_trace_n%d.csv',n)));
        fprintf('FLOOR n=%d tau=%.0e status=%s min ratio=%.3f\n',n,target,floorRun.status,min(floorRun.trace.ratio));
    end
    for method=["incremental-composed","radon-composed"]
        seconds=zeros(1,timing_reps(4));
        for rep=1:timing_reps(4)
            timer=tic;
            r=rksm_controller(a,p,"transient",tau,o,guard,shifts); require_accepted(r.status,"transport n="+n);
            red=a.project(p,r.V,r.W); sol=evaluate(a,red,times,method);
            seconds(rep)=toc(timer);
        end
        sol.V=r.V; sol.W=r.W; rows=record(a,rows,n,"RKSM-Leja",method,size(r.V,2),sol,ref,seconds,tau);
        if mode~="rksm", writetable(struct2table(rows),csv); end
    end
    if mode=="rksm", continue; end
    % Published EBA--BDF1 baseline with its own endpoint stopping rule.
    opts=a.options(p); opts.tol=1e-10; seconds=zeros(1,timing_reps(4));  % stopping tolerance of the 2020 paper
    for rep=1:timing_reps(4), timer=tic; r=a.adaptive(p,'EBA-BDF1',opts); seconds(rep)=toc(timer); end
    if isfield(r,'maxNewtonBackwardError'), fprintf('NEWTON max backward error %.3e\n',r.maxNewtonBackwardError); end
    sol=struct('V',r.V,'W',r.W,'Yhistory',r.Yhistory,'Y',r.Y,'time',r.time);
    rows=record(a,rows,n,"EBA (published)","BDF1, h=0.01",size(r.V,2),sol,ref,seconds);
    writetable(struct2table(rows),csv);
end
if mode=="rksm"
    old=readtable(csv,'TextType','string'); new=struct2table(rows);
    writetable([new;old(old.space~="RKSM-Leja",:)],csv);
end
end

function rows=record(a,rows,n,space,method,dim,sol,ref,seconds,tau)
if nargin<10, tau=NaN; end
if numel(sol.time)~=numel(ref.time) || max(abs(sol.time(:)-ref.time(:)))>1e-12
    sol=interpolate_to(sol,ref.time);
end
e=a.errors(sol,ref); budget=e.absoluteError./(1e-12+1e-7*e.referenceNorm);
[peak,j]=max(e.relativeError,[],'omitmissing');
row=struct('n',n,'space',space,'method',method,'tau',tau,'dimension',dim, ...
    'maxRelativeError',peak,'peakTime',e.time(j),'maxBudgetRatio',max(budget), ...
    'medianSeconds',timing_pick(seconds),'minSeconds',timing_pick(seconds,"min"), ...
    'maxSeconds',timing_pick(seconds,"max"));
rows=[rows;row];
fprintf('COMPLETE n=%d %s %s m=%d err=%.3e rho=%.3e time=%.3f\n',n,space,method,dim, ...
    peak,row.maxBudgetRatio,row.medianSeconds);
end

function sol=interpolate_to(sol,times)
% BDF1 returns its own grid; compare only on common nodes, NaN elsewhere.
[k,ell,~]=size(sol.Yhistory); history=NaN(k,ell,numel(times));
for j=1:numel(times)
    hit=find(abs(sol.time-times(j))<1e-12,1);
    if ~isempty(hit), history(:,:,j)=sol.Yhistory(:,:,hit); end
end
sol.Yhistory=history; sol.time=times;
end

function sol=evaluate(a,red,times,method)
if method=="radon-composed", sol=scaled_radon_nodes(red,times,"composed"); return; end
if method=="ode15s", sol=a.dense(red,times,1e-10); return; end
% Increment evaluation needs Ystar, Ac and Dc only, not the kernel L.
prep=a.prepare(red,false);
sol=incremental_centered_nodes(red,prep,times,"composed");
end

function ref=eba_reference(a,p,dim,times)
[V,W]=a.basis(p,'EBA',dim); red=a.project(p,V,W);
[k,ell]=size(red.Q); history=zeros(k,ell,numel(times));
for j=2:numel(times), v=a.mobius(red,[0 times(j)],1,0); history(:,:,j)=v.Y; end
ref=struct('V',V,'W',W,'Yhistory',history,'Y',history(:,:,end),'time',times,'dimension',dim);
end
