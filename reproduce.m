function folder = reproduce(mode)
%REPRODUCE Reproduce the numerical results of the paper.
%   reproduce("check")       check the recorded results in data/
%   reproduce("report")      write Tables 1-5 and Figures 1-2 from data/
%   reproduce("smoke")       run Algorithm 5.1 on two small problems
%   reproduce("dense")       rerun Tables 1-2 and Figure 1 (n = 600)
%   reproduce("evaluators")  rerun Table 5 (needs the transport reference)
%   reproduce("cached")      rerun everything, reusing the stored references
%   reproduce("all")         rerun everything, including the references
% Every mode except "check" and "smoke" writes into a new folder under
% reproduced/ and leaves data/ and figures/ unchanged; the folder is
% returned and printed.
if nargin < 1, mode = "report"; end
mode = string(mode);
modes = ["check","report","smoke","dense","evaluators","cached","all"];
assert(any(mode == modes),'ndre:mode','Use one of: %s.',strjoin(modes,', '));
root = fileparts(mfilename('fullpath'));
oldPath = path;
oldVisibility = get(groot,'defaultFigureVisible');
restore = onCleanup(@() restore_session(oldPath,oldVisibility)); %#ok<NASGU>
addpath(fullfile(root,'src'),fullfile(root,'experiments'),'-begin');
set(groot,'defaultFigureVisible','off');
ndre_scaled(true);
folder = root;
switch mode
    case "check"
        check_results;
        return
    case "smoke"
        run_smoke;
        return
end
folder = fullfile(root,'reproduced',char(mode+"_"+string(datetime('now','Format','yyyyMMdd_HHmmss'))));
mkdir(folder);
recorded = fullfile(root,'data');
if any(mode == ["report","dense","evaluators"])
    copy_files(recorded,fullfile(folder,'data'),'*.csv');
end
if any(mode == ["evaluators","cached"])
    references = dir(fullfile(recorded,'**','reference_*.mat'));
    if mode == "evaluators"
        references = references(strcmp({references.name},'reference_4000.mat'));
        assert(numel(references) == 1,'ndre:references', ...
            'Extract ndre-reference-transport.zip into the repository first.');
    else
        assert(numel(references) == 8,'ndre:references', ...
            'Extract the four ndre-reference-*.zip archives into the repository first.');
    end
    for j = 1:numel(references)
        target = fullfile(folder,'data',erase(references(j).folder,[recorded filesep]));
        if ~exist(target,'dir'), mkdir(target); end
        copyfile(fullfile(references(j).folder,references(j).name),target);
    end
end
ndre_paths(folder);
switch mode
    case "dense"
        run_bound_experiments("all",600);
        run_bound_history(600,1e-7);
        plot_bound_history;
    case "evaluators"
        run_evaluator_comparison;
    case {"cached","all"}
        run_all_experiments(600);
        plot_bound_history;
        plot_floor_traces;
    case "report"
        plot_bound_history;
        plot_floor_traces;
end
check_results;
export_tables;
fprintf('Output folder: %s\n',folder);
end

function copy_files(source,target,pattern)
files = dir(fullfile(source,'**',pattern));
for j = 1:numel(files)
    destination = fullfile(target,erase(files(j).folder,source));
    if ~exist(destination,'dir'), mkdir(destination); end
    copyfile(fullfile(files(j).folder,files(j).name),destination);
end
end

function restore_session(searchPath,visibility)
ndre_paths("");
path(searchPath);
set(groot,'defaultFigureVisible',visibility);
end
