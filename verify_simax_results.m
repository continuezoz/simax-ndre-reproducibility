function summary = verify_simax_results(root)
%VERIFY_SIMAX_RESULTS Audit recorded outcomes; this mode does not run solvers.
if nargin < 1, root = fileparts(mfilename('fullpath')); end
res = fullfile(root,'results_scaled');
read = @(name) readtable(fullfile(res,name),'TextType','string','VariableNamingRule','preserve');
D = read('bounds/dense_bounds.csv');
assert(height(D)==24 && all(D.status=="OK_ERROR_ESTIMATE"),'simax:dense','Expected 24 accepted dense tests.');
assert(all(isfinite(D.projOverEstimateMax)) && max(D.projOverEstimateMax)<1, ...
    'simax:bound','Projection error exceeds the archived bound.');
assert(max(D.projOverTarget)<1,'simax:target','The dense target was not met.');
W = read('bounds/wellposed.csv');
S = W(W.variant=="scaled",:);
assert(height(W)==16 && all(S.stationaryFailures==0) && all(S.flowFailures==0), ...
    'simax:wellposed','Scaled stationary/flow checks failed.');
assert(all(S.maxNorm2Ystar<1) && all(S.minMargin>0),'simax:unitBall','Scaled states must lie in the unit ball.');
C = read('complete_solves/complete_solves.csv');
assert(height(C)==12 && numel(unique(C.n))==4,'simax:singleState','Incomplete single-state records.');
I = C(C.method=="incremental-composed",:);
assert(all(isfinite(I.maxRelativeError)) && all(I.maxRelativeError<1e-7), ...
    'simax:singleStateError','Single-state accuracy is outside the reported range.');
assert(all(I.tau(I.n<40000)==1e-7) && I.tau(I.n==40000)==1e-6, ...
    'simax:tolerance','Unexpected controller tolerances.');
M = read('multistate_transport/multistate_run.csv');
I = M(M.method=="incremental-composed",:);
B = M(M.space=="EBA (published)",:);
assert(height(M)==12 && height(I)==4 && all(I.status=="OK_ERROR_ESTIMATE"), ...
    'simax:multistate','Incomplete accepted multistate results.');
assert(all(I.maxRelativeError<1e-7) && height(B)==4 && all(B.status=="BREAKDOWN"), ...
    'simax:baseline','Unexpected multistate error/baseline outcomes.');
E = read('evaluator_comparison/transport_evaluators.csv');
assert(height(E)==4 && all(isfinite(E.maxRelativeError)),'simax:evaluators','Missing evaluator results.');
F = read('floor/floor_summary.csv');
assert(height(F)==4 && all(F.plateauMedian>0),'simax:floor','Missing accuracy-floor traces.');
for name = ["transport","multistate"]
    H = read("bounds/history_"+name+".csv");
    assert(all(isfinite(H.projectionError)) && all(H.projectionError<=H.estimate+1e-18), ...
        'simax:history','Projection-error history exceeds the estimate.');
end
summary = struct('denseCases',height(D),'wellposedRecords',height(W), ...
    'singleStateRows',height(C),'multistateRows',height(M), ...
    'evaluatorRows',height(E),'maxProjectionErrorOverBound',max(D.projOverEstimateMax));
disp(summary);
fprintf('Recorded-results audit passed (no solver rerun in this mode).\n');
end
