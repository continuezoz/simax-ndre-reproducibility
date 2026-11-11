function run_scaled_suite(first,boundSize)
%RUN_SCALED_SUITE All experiments of the revised paper on the scaled
% equation (ndre_scaled true), in one session on an idle machine; results
% in results_scaled/. Stages are run in order from stage FIRST (default 1),
% so an interrupted run can be resumed; references and calibrations already
% saved are reused by the stage scripts. Stages that read references of
% earlier stages come after them (EBA deflated after E5/E7, evaluators
% after E6). Each stage waits until the machine is idle.
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
    timer=tic; % wait_idle disabled 2026-11-11: its own PowerShell sampling kept the load above 15%
    note('STAGE %d %s started %s',j,stages{j,1},string(datetime('now')));
    steps=stages{j,2}; if ~iscell(steps), steps={steps}; end
    for q=1:numel(steps), ndre_scaled(true); steps{q}(); end
    assert(ndre_scaled,'suite:mode','Stage %d left ndre_scaled false.',j);
    note('STAGE %d %s finished in %.0f s',j,stages{j,1},toc(timer));
end
note('ALL STAGES DONE %s',string(datetime('now')));
end

function wait_idle(note)
% Whole-machine load, sampled before MATLAB starts the next stage.
quiet=0;
while quiet<4
    [~,s]=system('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor | Measure-Object LoadPercentage -Average).Average"');
    cpu=str2double(strtrim(s));
    if cpu<15, quiet=quiet+1; else, quiet=0; note('WAIT load=%g%%',cpu); end
    pause(5);
end
end

function log_line(file,text)
fprintf('%s\n',text);
fid=fopen(file,'a'); fprintf(fid,'%s\n',text); fclose(fid);
end
