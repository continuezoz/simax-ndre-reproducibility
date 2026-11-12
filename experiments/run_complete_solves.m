function run_complete_solves(sizes)
%RUN_COMPLETE_SOLVES Table 3: complete solves of the near-critical
% single-state model for n = 4000, 10000, 20000, 40000.
% Timed work: basis construction, all candidate checks, projection, and the
% evaluation of all output states. References and error checks are untimed.
% The projected solver is Algorithm 5.1 (rksm_bound_controller) followed by
% composed increment or scaled Radon evaluation; the baseline is EBA--BDF1.
% Output: data/complete_solves/complete_solves.csv. References are stored
% as data/complete_solves/reference_<n>.mat and reused when present.
if nargin<1, sizes=[4000 10000 20000 40000]; end
out=fullfile(ndre_paths().data,'complete_solves');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'complete_solves.log'));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
a=ndre_core; times=unique([0,logspace(-8,-2,21),0:.01:1]);
rows=struct([]); csv=fullfile(out,'complete_solves.csv');
for n=sizes
    p=a.problem('transport',n);
    file=fullfile(out,sprintf('reference_%d.mat',n));
    if exist(file,'file'), s=load(file,'ref'); ref=s.ref;
    else
        % EBA spaces of dimension 240 (check: 280); 360 (check: 440) for n=40000.
        dims=[240 280]; if n==40000, dims=[360 440]; end
        ref=eba_reference(a,p,dims(1),times); fine=eba_reference(a,p,dims(2),times);
        c=a.errors(ref,fine); ref.spaceRefinement=max(c.absoluteError./(1e-12+1e-7*c.referenceNorm));
        save(file,'ref','-v7.3');
    end
    fprintf('REFERENCE n=%d refinement=%.3e\n',n,ref.spaceRefinement);
    o=struct('tf',1,'checkDimension',10);
    shifts=pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),319),"leja");
    % For n=40000 the computed estimate cannot fall below about 2.4e-7
    % (Figure 2), so tau=1e-6 is used there.
    tau=1e-7; if n>=40000, tau=1e-6; end
    for method=["incremental-composed","radon-composed"]
        timer=tic;
        r=rksm_bound_controller(a,p,1e-12,tau,o,shifts); require_accepted(r.status,"transport n="+n);
        red=a.project(p,r.V,r.W); sol=evaluate(a,red,times,method);
        seconds=toc(timer);
        sol.V=r.V; sol.W=r.W; rows=record(a,rows,n,"RKSM-Leja",method,size(r.V,2),sol,ref,seconds,tau,NaN);
        writetable(struct2table(rows),csv);
    end
    % EBA--BDF1 with the endpoint stopping rule and tolerance of Angelova et al.
    opts=a.options(p); opts.tol=1e-10;
    timer=tic; r=a.adaptive(p,'EBA-BDF1',opts); seconds=toc(timer);
    sol=struct('V',r.V,'W',r.W,'Yhistory',r.Yhistory,'Y',r.Y,'time',r.time);
    rows=record(a,rows,n,"EBA (published)","BDF1, h=0.01",size(r.V,2),sol,ref,seconds,NaN, ...
        r.maxNewtonBackwardError);
    writetable(struct2table(rows),csv);
end
end

function rows=record(a,rows,n,space,method,dim,sol,ref,seconds,tau,newton)
if numel(sol.time)~=numel(ref.time) || max(abs(sol.time(:)-ref.time(:)))>1e-12
    sol=interpolate_to(sol,ref.time);
end
e=a.errors(sol,ref); budget=e.absoluteError./(1e-12+1e-7*e.referenceNorm);
[peak,j]=max(e.relativeError,[],'omitmissing');
row=struct('n',n,'space',space,'method',method,'tau',tau,'dimension',dim, ...
    'maxRelativeError',peak,'peakTime',e.time(j),'maxBudgetRatio',max(budget), ...
    'seconds',seconds,'newtonBackwardError',newton);
rows=[rows;row];
fprintf('COMPLETE n=%d %s %s m=%d err=%.3e rho=%.3e time=%.3f\n',n,space,method,dim, ...
    peak,row.maxBudgetRatio,seconds);
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
