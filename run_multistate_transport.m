function run_multistate_transport(mode,sizes)
%RUN_MULTISTATE_TRANSPORT E7: NDRE (1.8) of the one-dimensional multistate
% transport model of Li, Chu, Juang and Lin (2011), R(0)=0, x in [0,1],
% rank-10 scattering (multistate_problem). Compares the projected solver
% (deflated block RKSM with Leja poles + composed increment / scaled Radon
% evaluation), the published EBA--BDF1 solver, and the low-rank ExpEuler
% and Erow3 schemes, all against an independent rational reference.
% Timed work as in run_complete_solves; references and error checks untimed.
% Mode "rksm" reruns only the projected solver rows and replaces them in
% multistate_run.csv.
if nargin<1, mode="run"; end
if nargin<2, sizes=[10000 100000]; end
root=fileparts(mfilename('fullpath'));
out=fullfile(root,results_dir,'multistate_transport');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,char("run_"+string(mode)+".log")));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
a=nla_core; grid=0:.01:1; times=unique([0,logspace(-8,-2,21),grid]); target=1e-7;
onGrid=ismembertol(times,grid,1e-12);
rows=struct([]); csv=fullfile(out,char("multistate_"+string(mode)+".csv"));
fullOrder=false; % 2026-11-11: full-order ExpEuler/Erow3 rows dropped from the paper
if mode=="rksm", csv=fullfile(out,'multistate_rksm.csv'); end
for n=sizes
    for gamma=[0.9999 0.5]
        [p,op]=multistate_problem(n,10,gamma,0.1);
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
        guard=logspace(log10(.01/hi),0,49);
        for method=["incremental-composed","radon-composed"]
            seconds=zeros(1,timing_reps(4));
            for rep=1:timing_reps(4)
                timer=tic;
                r=block_controller(a,p,target,1,guard,shifts); require_accepted(r.status,"multistate n="+n);
                red=a.project(p,r.V,r.W); sol=evaluate(a,red,times,method);
                seconds(rep)=toc(timer);
            end
            sol.V=r.V; sol.W=r.W;
            rows=record(a,rows,n,gamma,"RKSM-Leja",method,size(r.V,2),sol,ref,onGrid,timing_pick(seconds),r.status);
            writetable(struct2table(rows),csv);
        end
        if mode=="rksm", continue; end
        opts=a.options(p); opts.tol=1e-10; opts.outputTimes=grid; seconds=zeros(1,timing_reps(2)); ok=true;
        try
            for rep=1:timing_reps(2), timer=tic; r=a.adaptive(p,'EBA-BDF1',opts); seconds(rep)=toc(timer); end
            if isfield(r,'maxNewtonBackwardError'), fprintf('NEWTON max backward error %.3e\n',r.maxNewtonBackwardError); end
        catch ME, ok=false; fprintf('EBA-BDF1 failed: %s\n',ME.message);
        end
        if ok && ~r.converged
            % The first extended block deflates (dependent rank-r images), so
            % the published solver stops without meeting its own criterion.
            rows=add_row(rows,n,gamma,"EBA (published)","BDF1, h=0.01",size(r.V,2),NaN,NaN,timing_pick(seconds),"BREAKDOWN");
            writetable(struct2table(rows),csv);
        elseif ok
            sol=struct('V',r.V,'W',r.W,'Yhistory',r.Yhistory,'Y',r.Y,'time',r.time);
            rows=record(a,rows,n,gamma,"EBA (published)","BDF1, h=0.01",size(r.V,2),sol,ref,onGrid,timing_pick(seconds),"");
            writetable(struct2table(rows),csv);
        end
        if fullOrder && n==min(sizes)
            for scheme=["erow3","expeuler"]
                for j=0:2
                    h=0.01*2^-j; timer=tic; sol=li2021_lowrank_exp(op,grid,h,scheme); sec=toc(timer);
                    [rel,abserr,refnorm]=factor_errors(sol,ref,grid);
                    rows=add_row(rows,n,gamma,"full order",scheme+", h="+sprintf('%.4g',h), ...
                        max(sol.ranks),max(rel,[],'omitmissing'),max(abserr(:)./(1e-12+1e-7*refnorm(:))),sec,"");
                    writetable(struct2table(rows),csv);
                end
            end
        end
    end
end
if mode=="rksm"
    main=fullfile(out,'multistate_run.csv'); old=readtable(main,'TextType','string');
    new=struct2table(rows); keep=old.space~="RKSM-Leja" | ~ismember(old.n,sizes);
    writetable([new;old(keep,:)],main);
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

function [rel,abserr,refnorm]=factor_errors(sol,ref,grid)
% Errors of factored states L*S*U' against V*Y*W' without forming n x n matrices.
rel=NaN(numel(grid),1); abserr=zeros(numel(grid),1); refnorm=zeros(numel(grid),1);
for j=1:numel(grid)
    k=find(abs(ref.time-grid(j))<1e-12,1); Y=ref.Yhistory(:,:,k);
    f=sol.factors{j}; [~,RV]=qr([f{1},ref.V],0); [~,RW]=qr([f{3},ref.W],0);
    abserr(j)=norm(RV*blkdiag(f{2},-Y)*RW','fro'); refnorm(j)=norm(Y,'fro');
    if refnorm(j)>0, rel(j)=abserr(j)/refnorm(j); end
end
end

function r=block_controller(a,p,gate,tf,guard,shifts)
if isfield(p,'nu')
    % Scaled problem: Algorithm 1, candidate spaces tested after every block.
    r=rksm_bound_controller(a,p,1e-12,gate,struct('tf',tf,'checkDimension',0),shifts);
    r.indicator=Inf; if height(r.trace)>0, r.indicator=r.trace.ratio(end); end; return
end
% Transient check on the composed-increment lifted residual: at tf, and
% once that passes, also on the guard nodes.
timer=tic; V=orth(p.F); W=orth(p.G); lastV=V; lastW=W;
accepted=false; status="MAX_DIMENSION"; indicator=Inf;
for shift=shifts
    [V,lastV]=append(V,p.Ashift(lastV,shift));
    [W,lastW]=append(W,p.DTshift(lastW,shift));
    red=a.project(p,V,W);
    try, prep=a.prepare(red,false); catch ME, status=string(ME.identifier); continue; end
    [~,RA]=qr(p.Aapply(V)-V*red.A,0); [~,RD]=qr(p.DTapply(W)-W*red.D',0);
    z=increment_defects(red,prep,tf,RA,RD); indicator=z(1);
    if indicator<gate
        z=increment_defects(red,prep,unique([guard,tf]),RA,RD); indicator=max(z(:,1));
    end
    if indicator<gate, accepted=true; status="OK"; break; end
    if isempty(lastV) && isempty(lastW), status="INVARIANT"; break; end
    if isempty(lastV), lastV=V(:,end); end
    if isempty(lastW), lastW=W(:,end); end
end
r=struct('V',V,'W',W,'accepted',accepted,'status',status,'indicator',indicator,'seconds',toc(timer));
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
