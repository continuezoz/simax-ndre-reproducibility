function run_bound_experiments(stage,n)
%RUN_BOUND_EXPERIMENTS Tables 1 and 2 and the stationary bounds of section
% 6.2, on problems small enough for a dense full-order reference.
%   "dense"       Algorithm 5.1 for tau = 1e-4,...,1e-9: accepted dimensions,
%                 true error against the target and against the estimate.
%   "stationary"  ||Xstar - V*Ystar*W'||_F and the bound omega/(2*nubar) of
%                 Theorem 4.2 for growing spaces.
%   "wellposed"   every candidate space up to dimension 80: failures of the
%                 stationary solve and of the flow evaluation, closed-loop
%                 margin and ||Ystar||_2, for the scaled equation and for
%                 the unscaled equation with Euclidean orthonormal bases.
% Output: data/bounds/*.csv.
if nargin<1, stage="all"; end
if nargin<2, n=600; end
stage=string(stage);
out=fullfile(ndre_paths().data,'bounds'); if ~exist(out,'dir'), mkdir(out); end
previous=ndre_scaled; restore=onCleanup(@()ndre_scaled(previous)); %#ok<NASGU>
ndre_scaled(false);   % unscaled constructors; ndre_setup is called explicitly
a=ndre_core; cases=problem_cases(a,n);
if any(stage==["all","dense"]), dense_stage(a,cases,out); end
if any(stage==["all","stationary"]), stationary_stage(a,cases,out); end
if any(stage==["all","wellposed"]), wellposed_stage(a,cases,out); end
end

function cases=problem_cases(a,n)
cases=struct('label',{},'p0',{},'p',{},'lo',{},'hi',{},'block',{});
for par=[0.5 0.5;0.9999 1e-8]'
    p0=a.problem("transport",n,n,[],par');
    hi=p0.shiftRange(2)*(1+par(2))/(1-par(2));
    cases(end+1)=struct('label',sprintf("transport c=%g alpha=%g",par),'p0',p0, ...
        'p',ndre_setup(p0),'lo',p0.shiftRange(1),'hi',hi,'block',10); %#ok<AGROW>
end
for gamma=[0.5 0.9999]
    p0=multistate_problem(n,10,gamma,0.1);
    cases(end+1)=struct('label',sprintf("multistate gamma=%g",gamma),'p0',p0, ...
        'p',ndre_setup(p0),'lo',min(p0.factors.kp),'hi',max(p0.factors.km),'block',0); %#ok<AGROW>
end
end

function [full,fp]=dense_reference(a,p)
I=eye(p.n); full=a.project(p,I,I); fp=a.prepare(full,false);
end

function dense_stage(a,cases,out)
rows=struct([]);
for c=cases
    p=c.p; [full,fp]=dense_reference(a,p);
    shifts=pole_order(-logspace(log10(c.lo),log10(c.hi),159),"leja");
    o=struct('tf',1,'checkDimension',c.block); ref=[];
    for tau=10.^(-4:-1:-9)
        r=rksm_bound_controller(a,p,1e-12,tau,o,shifts); e=r.estimate; q=numel(r.nodes);
        if isempty(ref), ref=incremental_centered_nodes(full,fp,r.nodes(2:end),"composed"); end
        err=zeros(q,1); refn=err; perr=err; red=a.project(p,r.V,r.W);
        for j=2:q
            X=ref.Yhistory(:,:,j-1); refn(j)=norm(X,'fro');
            err(j)=norm(X-r.V*e.Yhistory(:,:,j)*r.W','fro');
            % projection error: reduced flow by scaled Radon composition
            mr=a.mobius(red,[0 r.nodes(j)],1,0); perr(j)=norm(X-r.V*mr.Y*r.W','fro');
        end
        pos=2:q; [pm,pj]=max(perr(pos)./e.epsilon(pos));
        % Compare only where the estimate is above the rounding level of
        % the reference and of the evaluator.
        big=find(e.epsilon>1e-11); ratio=err(big)./e.epsilon(big); [m,jm]=max(ratio);
        row=struct('case',c.label,'n',p.n,'nu',p.nu,'nubar',p.nubar,'fullMargin',fp.minClosedLoopRealPart, ...
            'tau',tau,'status',r.status,'k',size(r.V,2),'ell',size(r.W,2),'nodes',q-1, ...
            'errorOverTarget',max(err./(1e-12+tau*refn)),'errorOverEstimateMax',m, ...
            'timeAtMax',r.nodes(big(jm)),'estimateAtMax',e.epsilon(big(jm)), ...
            'errorOverEstimateMin',min(ratio),'errorOverEstimateFinal',err(end)/e.epsilon(end), ...
            'maxNorm2Y',max(e.norm2Y),'seconds',r.seconds, ...
            'projOverEstimateMax',pm,'timeAtProjMax',r.nodes(pos(pj)), ...
            'projOverEstimateFinal',perr(end)/e.epsilon(end), ...
            'evalErrorOverNormMax',max(abs(err(pos)-perr(pos))./refn(pos)), ...
            'projOverTarget',max(perr(pos)./(1e-12+tau*refn(pos))));
        rows=[rows;row]; %#ok<AGROW>
        fprintf('DENSE %s tau=%.0e %s k=%d err/target=%.2e err/est max %.3f (t=%.2e, est=%.1e) min %.2e final %.2e\n', ...
            c.label,tau,r.status,row.k,row.errorOverTarget,m,row.timeAtMax,row.estimateAtMax,row.errorOverEstimateMin,row.errorOverEstimateFinal);
        writetable(struct2table(rows),fullfile(out,'dense_bounds.csv'));
    end
end
end

function stationary_stage(a,cases,out)
rows=struct([]);
for c=cases
    p=c.p; [~,fp]=dense_reference(a,p); Xstar=fp.Ystar;
    shifts=pole_order(-logspace(log10(c.lo),log10(c.hi),159),"leja");
    V=orthonormal(zeros(p.n,0),p.F); W=orthonormal(zeros(p.p,0),p.G); b=size(V,2);
    for shift=shifts
        V=orthonormal(V,p.Ashift(V(:,end-b+1:end),shift)); W=orthonormal(W,p.DTshift(W(:,end-b+1:end),shift));
        k=size(V,2); if k>90, break; end
        if c.block>0 && mod(k,c.block)~=0, continue; end
        red=a.project(p,V,W); prep=a.prepare(red,false); Y=prep.Ystar;
        [~,RA]=qr(p.Aapply(V)-V*red.A,0); [~,RD]=qr(p.DTapply(W)-W*red.D',0);
        N=-red.A*Y-Y*red.D+Y*red.S*Y+red.Q;
        [U,sv,Z]=svd(Y,'econ'); sv=diag(sv); g=1./sqrt(1-sv.^2)-1;
        left=@(M) M+U*(g.*(U'*M)); right=@(M) M+((M*Z).*g')*Z';
        omega=sqrt(norm(right(left(N)),'fro')^2+norm(right(RA*Y),'fro')^2+norm(left(Y*RD'),'fro')^2);
        plain=sqrt(norm(N,'fro')^2+norm(RA*Y,'fro')^2+norm(Y*RD','fro')^2);
        err=norm(Xstar-V*Y*W','fro'); bound=omega/(2*p.nubar);
        row=struct('case',c.label,'n',p.n,'nu',p.nu,'nubar',p.nubar,'k',k,'norm2Ystar',max(sv), ...
            'norm2Xstar',norm(Xstar),'error',err,'bound',bound,'errorOverBound',err/bound, ...
            'weightedOverPlainResidual',omega/plain,'relativeError',err/norm(Xstar,'fro'));
        rows=[rows;row]; %#ok<AGROW>
        fprintf('STATIONARY %s k=%d err=%.3e bound=%.3e ratio=%.3e weight=%.2f\n',c.label,k,err,bound,err/bound,omega/plain);
    end
    writetable(struct2table(rows),fullfile(out,'stationary_bounds.csv'));
end
end

function wellposed_stage(a,cases,out)
rows=struct([]);
for c=cases
    for variant=["scaled","unscaled"]
        if variant=="scaled", p=c.p; else, p=c.p0; end
        for order=["leja","monotone"]
            shifts=pole_order(-logspace(log10(c.lo),log10(c.hi),159),order);
            V=orthonormal(zeros(p.n,0),p.F); W=orthonormal(zeros(p.p,0),p.G); b=size(V,2);
            tested=0; stationaryFailures=0; flowFailures=0; minMargin=Inf; maxNorm=0; reasons=strings(0,1);
            nodes=logspace(log10(.01/c.hi),0,60);
            for shift=shifts
                V=orthonormal(V,p.Ashift(V(:,end-b+1:end),shift)); W=orthonormal(W,p.DTshift(W(:,end-b+1:end),shift));
                if size(V,2)>80, break; end
                red=a.project(p,V,W); tested=tested+1;
                try
                    prep=a.prepare(red,false);
                catch ME
                    stationaryFailures=stationaryFailures+1; reasons(end+1,1)=string(ME.identifier); continue; %#ok<AGROW>
                end
                minMargin=min(minMargin,prep.minClosedLoopRealPart); maxNorm=max(maxNorm,norm(prep.Ystar));
                try
                    incremental_centered_nodes(red,prep,nodes,"composed");
                catch ME
                    flowFailures=flowFailures+1; reasons(end+1,1)=string(ME.identifier); %#ok<AGROW>
                end
            end
            row=struct('case',c.label,'variant',variant,'poleOrder',order,'spaces',tested, ...
                'stationaryFailures',stationaryFailures,'flowFailures',flowFailures, ...
                'minMargin',minMargin,'maxNorm2Ystar',maxNorm,'reasons',strjoin(unique(reasons),";"));
            rows=[rows;row]; %#ok<AGROW>
            fprintf('WELLPOSED %s %s %s spaces=%d stationary failures=%d flow failures=%d min margin=%.3e max||Ystar||=%.3f %s\n', ...
                c.label,variant,order,tested,stationaryFailures,flowFailures,minMargin,maxNorm,row.reasons);
        end
    end
    writetable(struct2table(rows),fullfile(out,'wellposed.csv'));
end
end

function Q=orthonormal(V,Z)
% Two-pass orthogonalization against V, thin QR of the new block.
Z=Z-V*(V'*Z); Z=Z-V*(V'*Z); [q,~]=qr(Z,0); q=q-V*(V'*q); [q,~]=qr(q,0); Q=[V,q];
end
