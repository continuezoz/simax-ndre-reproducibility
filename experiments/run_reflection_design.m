function run_reflection_design(n)
%RUN_REFLECTION_DESIGN Figure 3: design lengths of a multistate reflector.
% Multistate model (multistate_problem, rank 10, beta=0.1), R(0)=0. For an
% incident spectrum g (Gaussian in the energy variable, centre 0.5, width
% 0.05, sum 1) the albedo of a rod of length x is a(x)=e'R(x)g, the
% expected number of reflected particles per incident particle. It is
% nondecreasing in x. The design length x_p is the smallest x with
% a(x)=p*a(inf), p = 0.5, 0.9, 0.99.
%   Algorithm 5.1 (tau=1e-7 on [0,t_f]): a(inf) from the reduced stationary
%   solution, a(x) by direct increment evaluation, x_p by fzero. With the
%   bounds (4.4) and (4.5) on the scaled state,
%     |a(x)-a_m(x)| <= delta,  |a(inf)-a_m(inf)| <= deltaInf,
%   and x_p lies in [xLo,xHi] with a_m(xLo)+delta = p*(a_m(inf)-deltaInf)
%   and a_m(xHi)-delta = p*(a_m(inf)+deltaInf).
%   EBA--BDF1 with deflation (eba_deflated_bdf; the published method stops
%   at dimension 20 on this model, Table 4), h=0.01 as published and h=0.001
%   on the same space; x_p by linear interpolation of the grid values,
%   using the reference a(inf) (which the method does not provide) and,
%   alternatively, its own value at t_f.
%   Reference: deflated rational space from 100 Leja poles in
%   [min k^+/20, 50 max k^-], scaled Radon evaluation (as for Table 4).
% Output: data/reflection_design/reflection_design.csv, and for Figure 3
% albedo_g*.csv, albedo_bdf1_g*.csv (a(x)/a(inf)) and design_g*.csv (x_p
% for p in [0.3,0.995]: Algorithm 5.1, enclosure, EBA--BDF1 h=0.01, 0.001).
if nargin<1, n=10000; end
out=fullfile(ndre_paths().data,'reflection_design');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'reflection_design.log'));
cleanup=onCleanup(@()diary('off')); %#ok<NASGU>
a=ndre_core; levels=[0.5 0.9 0.99]; gammas=[0.5 0.9 0.9999]; tfs=[0.5 2 15];
u=linspace(0,1,n)'; g=exp(-(u-0.5).^2/(2*0.05^2)); g=g/sum(g); e=ones(n,1);
rows=struct([]); csv=fullfile(out,'reflection_design.csv');
for i=1:numel(gammas)
    gamma=gammas(i); tf=tfs(i);
    p=multistate_problem(n,10,gamma,0.1);
    wl=e./p.scaling.dA; wr=p.scaling.dD.*g;        % a = wl'*Xhat*wr
    lo=p.shiftRange(1); hi=p.shiftRange(2);
    % --- reference (untimed)
    rb=rational_basis(p,pole_order(-logspace(log10(lo/20),log10(50*hi),100),"leja"));
    redb=a.project(p,rb.V,rb.W); prepb=a.prepare(redb,false);
    lb=wl'*rb.V; rbv=rb.W'*wr; ainfRef=lb*prepb.Ystar*rbv;
    albRef=@(x) lb*radon_at(redb,x)*rbv;
    xRef=crossings(albRef,levels*ainfRef,tf);
    fprintf('gamma=%g tf=%g nu=%.4g nubar=%.4g ref dim=%d ainf=%.10f x=%s\n', ...
        gamma,tf,p.nu,p.nubar,size(rb.V,2),ainfRef,mat2str(xRef,8));
    % --- Algorithm 5.1
    timer=tic;
    shifts=pole_order(-logspace(log10(lo),log10(hi),79),"leja");
    r=rksm_bound_controller(a,p,1e-12,1e-7,struct('tf',tf,'checkDimension',0),shifts);
    require_accepted(r.status,"reflection gamma="+gamma);
    red=a.project(p,r.V,r.W); prep=a.prepare(red,false);
    lv=wl'*r.V; rv=r.W'*wr; ainf=lv*prep.Ystar*rv;
    alb=@(x) lv*increment_at(red,prep,x)*rv;
    xs=crossings(alb,levels*ainf,tf);
    seconds=toc(timer);
    % certified enclosure (untimed)
    [~,RA]=qr(p.Aapply(r.V)-r.V*red.A,0); [~,RD]=qr(p.DTapply(r.W)-r.W*red.D',0);
    weight=norm(wl)*norm(wr);
    delta=weight*max(r.estimate.epsilon);
    deltaInf=weight*stationary_omega(red,prep.Ystar,RA,RD)/(2*p.nubar);
    xLo=crossings(@(x)alb(x)+delta,levels*(ainf-deltaInf),tf);
    xHi=crossings(@(x)alb(x)-delta,levels*(ainf+deltaInf),tf);
    rows=add_row(rows,gamma,tf,"Algorithm 5.1",size(r.V,2),0,ainf,ainfRef,xs,xRef,xLo,xHi,seconds,delta,deltaInf);
    writetable(struct2table(rows),csv);
    % curves for Figure 3 (untimed): albedo on a fine grid and design
    % lengths for a range of fractions p
    xg=[0,logspace(-4,log10(tf),400)]; ag=[0,arrayfun(alb,xg(2:end))];
    pg=[0.3:0.005:0.99,0.995]; tag=sprintf('g%g',gamma);
    design=table(pg',pre_crossings(alb,xg,ag,pg'*ainf), ...
        pre_crossings(@(x)alb(x)+delta,xg,ag+delta,pg'*(ainf-deltaInf)), ...
        pre_crossings(@(x)alb(x)-delta,xg,ag-delta,pg'*(ainf+deltaInf)), ...
        'VariableNames',{'p','x','xLo','xHi'});
    writetable(table(xg',ag'/ainf,'VariableNames',{'x','albedoRatio'}), ...
        fullfile(out,['albedo_' tag '.csv']));
    % --- deflated EBA--BDF1
    opts=a.options(p); opts.tol=1e-10; opts.h=0.01; opts.tf=tf;
    opts.outputTimes=0:0.01:tf; opts.maxEBASteps=60; opts.maxDim=1200;
    timer=tic; rr=eba_deflated_bdf(a,p,opts); basisSeconds=toc(timer);
    redE=a.project(p,rr.V,rr.W); lvb=wl'*rr.V; rvb=rr.W'*wr;
    for h=[0.01 0.001]
        if h==0.01, hist=rr.Yhistory; tt=opts.outputTimes; secondsH=basisSeconds;
        else
            o3=opts; o3.h=h; o3.outputTimes=0:h:tf;
            timer=tic; s3=a.bdf(redE,o3); secondsH=basisSeconds+toc(timer);
            hist=s3.Yhistory; tt=o3.outputTimes;
        end
        ab=zeros(numel(tt),1);
        for j=1:numel(ab), ab(j)=lvb*hist(:,:,j)*rvb; end
        xb=grid_crossings(tt,ab,levels*ainfRef);
        xOwn=grid_crossings(tt,ab,levels*ab(end));
        design.("xBDF1_h"+round(1/h))=grid_crossings(tt,ab,pg'*ainfRef);
        if h==0.01
            writetable(table(tt',ab/ainfRef,'VariableNames',{'x','albedoRatio'}), ...
                fullfile(out,['albedo_bdf1_' tag '.csv']));
        end
        rows=add_row(rows,gamma,tf,"EBA-BDF1 deflated ("+rr.status+")",size(rr.V,2),h,ab(end),ainfRef,xb,xRef,xOwn,NaN(size(xb)),secondsH,NaN,NaN);
        writetable(struct2table(rows),csv);
    end
    writetable(design,fullfile(out,['design_' tag '.csv']));
end
end

function x=pre_crossings(f,grid,vals,targets)
% As crossings, with f already evaluated on grid (grid(1)=0).
x=NaN(size(targets));
for i=1:numel(targets)
    j=find(vals>=targets(i),1);
    if isempty(j) || j==1, continue; end
    x(i)=fzero(@(s)f(s)-targets(i),[grid(j-1) grid(j)],optimset('TolX',1e-15));
end
end

function x=crossings(f,targets,tf)
% Smallest x in (0,tf] with f(x)=target; f is nondecreasing up to rounding.
x=NaN(size(targets)); grid=[0,logspace(-6,log10(tf),200)];
vals=arrayfun(f,grid(2:end)); vals=[0,vals];
for i=1:numel(targets)
    j=find(vals>=targets(i),1);
    if isempty(j), continue; end
    if j==1, x(i)=0; continue; end
    x(i)=fzero(@(s)f(s)-targets(i),[grid(j-1) grid(j)],optimset('TolX',1e-15));
end
end

function x=grid_crossings(tt,ab,targets)
x=NaN(size(targets));
for i=1:numel(targets)
    j=find(ab>=targets(i),1);
    if isempty(j) || j==1, continue; end
    x(i)=tt(j-1)+(targets(i)-ab(j-1))/(ab(j)-ab(j-1))*(tt(j)-tt(j-1));
end
end

function w=stationary_omega(red,Y,RA,RD)
% Weighted residual (4.3) of the constant matrix Y, as in error_estimate.
N=-red.A*Y-Y*red.D+Y*red.S*Y+red.Q;
[U,sv,Z]=svd(Y,'econ'); sv=diag(sv); assert(max(sv)<1,'ndre:ball','Ystar not in the ball.');
gg=1./sqrt(1-sv.^2)-1;
left=@(M) M+U*(gg.*(U'*M)); right=@(M) M+((M*Z).*gg')*Z';
w=sqrt(norm(right(left(N)),'fro')^2+norm(right(RA*Y),'fro')^2+norm(left(Y*RD'),'fro')^2);
end

function Y=increment_at(red,prep,x)
s=incremental_centered_nodes(red,prep,[0 x],"direct"); Y=s.Yhistory(:,:,end);
end

function Y=radon_at(red,x)
s=scaled_radon_nodes(red,[0 x],"composed"); Y=s.Yhistory(:,:,end);
end

function r=rational_basis(p,shifts)
[V,nv]=app(zeros(p.n,0),p.F); [W,nw]=app(zeros(p.p,0),p.G);
for s=shifts
    if nv>0, [V,nv]=app(V,p.Ashift(V(:,end-nv+1:end),s)); end
    if nw>0, [W,nw]=app(W,p.DTshift(W(:,end-nw+1:end),s)); end
    if nv+nw==0, break; end
end
r=struct('V',V,'W',W);
end

function [Q,added]=app(V,Z)
scale=norm(Z,'fro'); Z=Z-V*(V'*Z); Z=Z-V*(V'*Z);
[q,R,~]=qr(Z,0); d=abs(diag(R)); added=sum(d>1e-10*max(scale,eps));
q=q(:,1:added); q=q-V*(V'*q); [q,~]=qr(q,0); Q=[V,q];
end

function rows=add_row(rows,gamma,tf,method,dim,h,ainf,ainfRef,x,xRef,x2,x3,seconds,delta,deltaInf)
rel=abs(x-xRef)./xRef;
row=struct('gamma',gamma,'tf',tf,'method',method,'dimension',dim,'h',h, ...
    'ainf',ainf,'ainfRef',ainfRef,'x50',x(1),'x90',x(2),'x99',x(3), ...
    'x50Ref',xRef(1),'x90Ref',xRef(2),'x99Ref',xRef(3), ...
    'rel50',rel(1),'rel90',rel(2),'rel99',rel(3), ...
    'aux50',x2(1),'aux90',x2(2),'aux99',x2(3),'hi50',x3(1),'hi90',x3(2),'hi99',x3(3), ...
    'seconds',seconds,'delta',delta,'deltaInf',deltaInf);
rows=[rows;row];
fprintf('  %-34s k=%3d h=%-6g a=%.10f x=%s relerr=%s aux=%s hi=%s t=%.2f delta=%.2e deltaInf=%.2e\n', ...
    method,dim,h,ainf,mat2str(x,8),mat2str(rel,2),mat2str(x2,8),mat2str(x3,8),seconds,delta,deltaInf);
end
