function run_floor_traces(sizes,tau)
%RUN_FLOOR_TRACES Figure 2: rounding floor of the error estimate. Algorithm
% 5.1 on the near-critical single-state model, continued to the last pole
% without stopping (o.noStop), for several n.
% Output in data/floor/: trace_n<n>.csv, where ratio is
% max_j {eps_j/(tau_a+tau*||Y(sigma_j)||_F)}, and floor_summary.csv with
% u*||A|| (largest diagonal entry of the scaled A).
if nargin<1, sizes=[4000 10000 20000 40000]; end
if nargin<2, tau=1e-7; end
out=fullfile(ndre_paths().data,'floor');
if ~exist(out,'dir'), mkdir(out); end
ndre_scaled(true); a=ndre_core; rows=struct([]);
for n=sizes
    p=a.problem('transport',n);
    shifts=pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),319),"leja");
    r=rksm_bound_controller(a,p,1e-12,tau,struct('tf',1,'checkDimension',10,'noStop',true),shifts);
    writetable(r.trace,fullfile(out,sprintf('trace_n%d.csv',n)));
    tail=r.trace.ratio(end-9:end);
    row=struct('n',n,'normA',p.shiftRange(2),'uNormA',eps/2*p.shiftRange(2), ...
        'plateauMedian',tau*median(tail),'plateauMin',tau*min(r.trace.ratio),'nu',p.nu);
    rows=[rows;row]; %#ok<AGROW>
    writetable(struct2table(rows),fullfile(out,'floor_summary.csv'));
    fprintf('FLOOR n=%d u||A||=%.3e plateau(median of last 10)*tau=%.3e min=%.3e\n',n,row.uNormA,row.plateauMedian,row.plateauMin);
end
end
