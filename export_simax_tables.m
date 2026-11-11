function export_simax_tables(root)
%EXPORT_SIMAX_TABLES Derive all five numerical tables from recorded CSVs.
if nargin < 1, root = fileparts(mfilename('fullpath')); end
res = fullfile(root,'results_scaled');
out = fullfile(root,'generated','tables');
if ~exist(out,'dir'), mkdir(out); end
read = @(name) readtable(fullfile(res,name),'TextType','string','VariableNamingRule','preserve');
W = read('bounds/wellposed.csv');
D = read('bounds/dense_bounds.csv');
models = unique(W.("case"),'stable');
rows = cell(0,8);
for model = models'
    for variant = ["scaled","unscaled"]
        R = W(W.("case") == model & W.variant == variant,:);
        nu = D.nu(find(D.("case") == model,1));
        rows(end+1,:) = {model,variant,nu,R.stationaryFailures(R.poleOrder=="leja"), ...
            R.stationaryFailures(R.poleOrder=="monotone"), ...
            min(R.minMargin),max(R.maxNorm2Ystar),max(R.flowFailures)}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Problem','Equation','Nu','FailuresLeja', ...
    'FailuresIncreasing','MinClosedLoopMargin','MaxStationaryNorm','FlowFailures'});
emit(T,'table2_wellposed',out);
keep = any(abs(log10(D.tau)-[-4 -7 -9]) < 1e-8,2);
T = D(keep,{'case','tau','k','projOverTarget','projOverEstimateMax', ...
    'timeAtProjMax','projOverEstimateFinal','evalErrorOverNormMax'});
T.Properties.VariableNames{1} = 'Problem';
emit(T,'table3_dense_bounds',out);
C = read('complete_solves/complete_solves.csv');
newton = parse_newton(fullfile(res,'complete_solves','run_run.log'));
rows = cell(0,11);
for n = unique(C.n,'stable')'
    I = C(C.n==n & C.method=="incremental-composed",:);
    R = C(C.n==n & C.method=="radon-composed",:);
    B = C(C.n==n & C.space=="EBA (published)",:);
    assert(height(I)==1 && height(R)==1 && height(B)==1,'simax:table4','Incomplete n=%d.',n);
    rows(end+1,:) = {n,I.tau,I.dimension,I.maxRelativeError,I.medianSeconds, ...
        R.maxRelativeError,R.medianSeconds,B.dimension,B.maxRelativeError, ...
        B.medianSeconds,newton(n)}; %#ok<AGROW>
end
T = cell2table(rows,'VariableNames',{'N','Tau','Dimension','IncrementError', ...
    'IncrementSeconds','RadonError','RadonSeconds','BDF1Dimension','BDF1Error', ...
    'BDF1Seconds','MaxNewtonBackwardError'});
emit(T,'table4_single_state',out);
M = read('multistate_transport/multistate_run.csv');
pairs = unique(M(:,{'n','gamma'}),'rows','stable');
rows = cell(0,9);
for j = 1:height(pairs)
    subset = M(M.n==pairs.n(j) & M.gamma==pairs.gamma(j),:);
    I = subset(subset.method=="incremental-composed",:);
    R = subset(subset.method=="radon-composed",:);
    B = subset(subset.space=="EBA (published)",:);
    rows(end+1,:) = {I.n,I.gamma,I.dimension,I.maxRelativeError, ...
        I.seconds,R.maxRelativeError,R.seconds,B.dimension,B.status}; %#ok<AGROW>
end
T = cell2table(rows,'VariableNames',{'N','Gamma','Dimension','IncrementError', ...
    'IncrementSeconds','RadonError','RadonSeconds','BDF1Dimension','BDF1Status'});
emit(T,'table5_multistate',out);
E = read('evaluator_comparison/transport_evaluators.csv');
T = E(:,{'method','maxRelativeError','peakTime','medianSeconds','sylvesterResidual'});
emit(T,'table6_evaluators',out);
fprintf('Exported Tables 2--6 to %s\n',out);
end

function values = parse_newton(filename)
values = containers.Map('KeyType','double','ValueType','double');
lines = splitlines(string(fileread(filename)));
current = NaN;
for line = lines'
    token = regexp(line,'REFERENCE n=(\d+)','tokens','once');
    if ~isempty(token), current = str2double(token{1}); end
    token = regexp(line,'NEWTON max backward error ([\d.eE+-]+)','tokens','once');
    if ~isempty(token) && isfinite(current), values(current) = str2double(token{1}); end
end
end

function emit(T,name,folder)
writetable(T,fullfile(folder,[name '.csv']));
fid = fopen(fullfile(folder,[name '.tex']),'w');
assert(fid>=0,'simax:write','Cannot write the table.');
closer = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Numerical entries derived from the archived CSV files.\n');
fprintf(fid,'\\begin{tabular}{%s}\n\\hline\n',repmat('l',1,width(T)));
headers = replace(string(T.Properties.VariableNames),'_','\_');
fprintf(fid,'%s \\\\\n\\hline\n',strjoin(headers,' & '));
for i = 1:height(T)
    cells = strings(1,width(T));
    for j = 1:width(T)
        value = T{i,j};
        if iscell(value), value = value{1}; end
        if isnumeric(value)
            if isnan(value), cells(j) = "--";
            else, cells(j) = string(sprintf('%.6g',value)); end
        else
            cells(j) = replace(string(value),["_","&","%"],["\_","\&","\%"]);
        end
    end
    fprintf(fid,'%s \\\\\n',strjoin(cells,' & '));
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
end
