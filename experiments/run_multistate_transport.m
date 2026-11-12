function run_multistate_transport(sizes)
%RUN_MULTISTATE_TRANSPORT Table 4: the multistate transport model of Li,
% Chu, Juang and Lin (2011), R(0)=0, x in [0,1], rank-10 scattering
% (multistate_problem). Compares the projected solver
% (Algorithm 5.1 with Leja poles, every block tested, + composed increment /
% scaled Radon evaluation) and the published EBA--BDF1 solver against an
% independent rational reference.
% Timed work as in run_complete_solves; references and error checks untimed.
% Output: data/multistate_transport/multistate_run.csv. References are
% stored as data/multistate_transport/reference_n<n>_g<gamma>.mat.
if nargin<1, sizes=[10000 100000]; end
out=fullfile(ndre_paths().data,'multistate_transport');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'multistate_transport.log'));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
a=ndre_core; grid=0:.01:1; times=unique([0,logspace(-8,-2,21),grid]); target=1e-7;
onGrid=ismembertol(times,grid,1e-12);
rows=struct([]); csv=fullfile(out,'multistate_run.csv');
for n=sizes
    for gamma=[0.9999 0.5]
        p=multistate_problem(n,10,gamma,0.1);
        tag=sprintf('n%d_g%g',n,gamma); lo=p.shiftRange(1); hi=p.shiftRange(2);
        file=fullfile(out,['reference_' tag '.mat']);
        if exist(file,'file'), s=load(file,'ref'); ref=s.ref;
        else
            ref=rational_reference(a,p,times,[lo/20,50*hi],100);
            fine=rational_reference(a,p,times,[lo/50,100*hi],140);
            c=a.errors(ref,fine); ref.spaceRefinement=max(c.absoluteError./(1e-12+1e-7*c.referenceNorm));
            save(file,'ref','-v7.3');
        end
        fprintf('REFERENCE %s refinement=%.3e\n',tag,ref.spaceRefinement);
        shifts=pole_order(-logspace(log10(lo),log10(hi),79),"leja");
        for method=["incremental-composed","radon-composed"]
            timer=tic;
            r=rksm_bound_controller(a,p,1e-12,target,struct('tf',1,'checkDimension',0),shifts);
            require_accepted(r.status,"multistate n="+n);
            red=a.project(p,r.V,r.W); sol=evaluate(a,red,times,method);
            seconds=toc(timer);
            sol.V=r.V; sol.W=r.W;
            rows=record(a,rows,n,gamma,"RKSM-Leja",method,size(r.V,2),sol,ref,onGrid,seconds,r.status);
            writetable(struct2table(rows),csv);
        end
        opts=a.options(p); opts.tol=1e-10; opts.outputTimes=grid; ok=true;
        try
            timer=tic; r=a.adaptive(p,'EBA-BDF1',opts); seconds=toc(timer);
            if isfield(r,'maxNewtonBackwardError'), fprintf('NEWTON max backward error %.3e\n',r.maxNewtonBackwardError); end
        catch ME, ok=false; fprintf('EBA-BDF1 failed: %s\n',ME.message);
        end
        if ok && ~r.converged
            % The first extended block is nearly rank deficient, and the
            % published solver stops without meeting its own criterion.
            rows=add_row(rows,n,gamma,"EBA (published)","BDF1, h=0.01",size(r.V,2),NaN,NaN,seconds,"BREAKDOWN");
            writetable(struct2table(rows),csv);
        elseif ok
            sol=struct('V',r.V,'W',r.W,'Yhistory',r.Yhistory,'Y',r.Y,'time',r.time);
            rows=record(a,rows,n,gamma,"EBA (published)","BDF1, h=0.01",size(r.V,2),sol,ref,onGrid,seconds,"");
            writetable(struct2table(rows),csv);
        end
    end
end
end

function rows=record(a,rows,n,gamma,space,method,dim,sol,ref,onGrid,seconds,status)
if numel(sol.time)~=numel(ref.time) || max(abs(sol.time(:)-ref.time(:)))>1e-12
    sol=interpolate_to(sol,ref.time);
end
e=a.errors(sol,ref); budget=e.absoluteError./(1e-12+1e-7*e.referenceNorm);
rel=e.relativeError(:); budget=budget(:);
if space=="EBA (published)", rel=rel(onGrid); budget=budget(onGrid); end
rows=add_row(rows,n,gamma,space,method,dim,max(rel,[],'omitmissing'),max(budget,[],'omitmissing'),seconds,status);
end

function rows=add_row(rows,n,gamma,space,method,dim,err,rho,seconds,status)
row=struct('n',n,'gamma',gamma,'space',space,'method',method,'dimension',dim, ...
    'maxRelativeError',err,'maxBudgetRatio',rho,'seconds',seconds,'status',string(status));
rows=[rows;row];
fprintf('MULTISTATE n=%d gamma=%g %s %s k=%d err=%.3e rho=%.3e time=%.3f %s\n',n,gamma, ...
    space,method,dim,err,rho,seconds,status);
end

function sol=interpolate_to(sol,times)
[k,ell,~]=size(sol.Yhistory); history=NaN(k,ell,numel(times));
for j=1:numel(times)
    hit=find(abs(sol.time-times(j))<1e-12,1);
    if ~isempty(hit), history(:,:,j)=sol.Yhistory(:,:,hit); end
end
sol.Yhistory=history; sol.time=times;
end

function [V,Q]=append(V,Z)
% Two-pass block orthogonalization with deflation of dependent directions.
scale=max(norm(Z,'fro'),eps); Z=Z-V*(V'*Z); Z=Z-V*(V'*Z);
[q,R,~]=qr(Z,0); keep=abs(diag(R))>1e-10*scale; Q=q(:,keep); V=[V,Q];
end

function sol=evaluate(a,red,times,method)
if method=="radon-composed", sol=scaled_radon_nodes(red,times,"composed"); return; end
% Increment evaluation needs Ystar, Ac and Dc only, not the kernel L.
prep=a.prepare(red,false);
sol=incremental_centered_nodes(red,prep,times,"composed");
end

function ref=rational_reference(a,p,times,range,count)
shifts=pole_order(-logspace(log10(range(1)),log10(range(2)),count),"leja");
V=orth(p.F); W=orth(p.G); lastV=V; lastW=W;
for shift=shifts
    [V,lastV]=append(V,p.Ashift(lastV,shift)); [W,lastW]=append(W,p.DTshift(lastW,shift));
    if isempty(lastV), lastV=V(:,end); end
    if isempty(lastW), lastW=W(:,end); end
end
red=a.project(p,V,W); [k,ell]=size(red.Q); history=zeros(k,ell,numel(times));
for j=2:numel(times), v=a.mobius(red,[0 times(j)],1,0); history(:,:,j)=v.Y; end
ref=struct('V',V,'W',W,'Yhistory',history,'Y',history(:,:,end),'time',times,'dimension',[k,ell]);
fprintf('REFERENCE SPACE k=%d l=%d\n',k,ell);
end
