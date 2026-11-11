function folder = reproduce_simax(mode)
%REPRODUCE_SIMAX Reproduce the numerical results of the SIMAX manuscript.
% verify: inspect recorded CSVs; report: regenerate tables and MATLAB figures.
% smoke: solve small generated problems. dense/evaluators: numerical reruns.
% all: full numerical rerun with fresh references. cached: reuse references.
% Every mode that writes numerical output uses an isolated reproduced/ folder.
if nargin < 1, mode = "report"; end
mode = string(mode);
root = fileparts(mfilename('fullpath'));
if mode == "verify"
    verify_simax_results(root);
    folder = root;
    return
end
assert(any(mode == ["report","smoke","dense","evaluators","all","cached"]), ...
    'simax:mode','Use verify, report, smoke, dense, evaluators, all, or cached.');
folder = fullfile(root,'reproduced',char(mode+"_"+string(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))));
mkdir(folder);
sources = dir(fullfile(root,'*.m'));
for j = 1:numel(sources)
    copyfile(fullfile(sources(j).folder,sources(j).name),folder);
end
if ~any(mode == ["all","cached","smoke"])
    recorded = dir(fullfile(root,'results_scaled','**','*'));
    for j = 1:numel(recorded)
        [~,~,extension] = fileparts(recorded(j).name);
        if recorded(j).isdir || ~any(string(extension) == [".csv",".log"]), continue; end
        source = fullfile(recorded(j).folder,recorded(j).name);
        relative = erase(string(source),string(root)+filesep);
        destination = fullfile(folder,relative);
        if ~exist(fileparts(destination),'dir'), mkdir(fileparts(destination)); end
        copyfile(source,destination);
    end
end
quadrature = fullfile(root,'results','matlab_quadrature_20261007');
if exist(quadrature,'dir')
    destination = fullfile(folder,'results','matlab_quadrature_20261007');
    mkdir(destination);
    copyfile(fullfile(quadrature,'*.mat'),destination);
end
if any(mode == ["cached","evaluators"])
    references = dir(fullfile(root,'results_scaled','**','reference_*.mat'));
    if mode == "cached"
        assert(numel(references) == 8,'simax:references', ...
            'Extract the four reference-data ZIPs into the repository first (8 MAT files).');
    else
        references = references(strcmp({references.name},'reference_4000.mat'));
        assert(numel(references) == 1,'simax:references', ...
            'Extract simax-reference-transport-v1.0.0.zip into the repository first.');
    end
    for j = 1:numel(references)
        source = fullfile(references(j).folder,references(j).name);
        relative = erase(string(source),string(root)+filesep);
        destination = fullfile(folder,relative);
        if ~exist(fileparts(destination),'dir'), mkdir(fileparts(destination)); end
        copyfile(source,destination);
    end
end
oldDirectory = pwd;
oldPath = path;
oldVisibility = get(groot,'defaultFigureVisible');
restore = onCleanup(@() restore_session(oldDirectory,oldPath,oldVisibility)); %#ok<NASGU>
cd(folder);
addpath(folder,'-begin');
set(groot,'defaultFigureVisible','off');
ndre_scaled(true);
switch mode
    case "smoke"
        run_simax_smoke;
    case {"all","cached"}
        run_scaled_suite(1,600);
        run_bound_history(600,1e-7);
        run_floor_traces;
        verify_simax_results(folder);
        export_simax_tables(folder);
        plot_bound_history;
        plot_floor_traces;
    case "dense"
        run_bound_experiments("all",600);
        run_bound_history(600,1e-7);
        verify_simax_results(folder);
        export_simax_tables(folder);
        plot_bound_history;
    case "evaluators"
        run_evaluator_comparison("transport");
        verify_simax_results(folder);
        export_simax_tables(folder);
    case "report"
        verify_simax_results(folder);
        export_simax_tables(folder);
        plot_bound_history;
        plot_floor_traces;
end
fprintf('Output directory: %s\n',folder);
end

function restore_session(directory,searchPath,visibility)
cd(directory);
path(searchPath);
set(groot,'defaultFigureVisible',visibility);
end
