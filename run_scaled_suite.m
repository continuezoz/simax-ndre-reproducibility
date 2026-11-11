function run_scaled_suite(first,boundSize)
%RUN_SCALED_SUITE All experiments of the revised paper on the scaled
% equation (ndre_scaled true), in one session on an idle machine; results
% in results_scaled/. Stages are run in order from stage FIRST (default 1),
% so an interrupted run can be resumed; references and calibrations already
% saved are reused by the stage scripts. Stages that read references of
% earlier stages come after them (EBA deflated after E5/E7, evaluators
% after E6).
if nargin<1, first=1; end
if nargin<2, boundSize=600; end
ndre_scaled(true);
root=fileparts(mfilename('fullpath')); res=fullfile(root,results_dir);
if ~exist(res,'dir'), mkdir(res); end
logFile=fullfile(res,'run_scaled_suite.log');
note=@(varargin) log_line(logFile,sprintf(varargin{:}));
% 2026-11-11: reduced to the experiments that support the contributions;
% every solve is run once (timing_reps).
stages={ ...
    "B bounds and well-posedness", @() run_bound_experiments("all",boundSize); ...
    "E6 complete solves",          @() run_complete_solves("run"); ...
    "E7 multistate",               @() run_multistate_transport("run"); ...
    "EV evaluators",               @() run_evaluator_comparison("transport")};
for j=first:size(stages,1)
    timer=tic;
    note('STAGE %d %s started %s',j,stages{j,1},string(datetime('now')));
    steps=stages{j,2}; if ~iscell(steps), steps={steps}; end
    for q=1:numel(steps), ndre_scaled(true); steps{q}(); end
    assert(ndre_scaled,'suite:mode','Stage %d left ndre_scaled false.',j);
    note('STAGE %d %s finished in %.0f s',j,stages{j,1},toc(timer));
end
note('ALL STAGES DONE %s',string(datetime('now')));
end

function log_line(file,text)
fprintf('%s\n',text);
fid=fopen(file,'a'); fprintf(fid,'%s\n',text); fclose(fid);
end
