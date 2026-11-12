function p=ndre_paths(output)
%NDRE_PATHS Folders used by the experiments.
%   p = ndre_paths returns the repository root, the data folder, the figure
%   folder, the table folder and the quadrature folder. By default the
%   experiments read and write data/ and figures/ of the repository.
%   ndre_paths(folder) redirects data/, figures/ and tables/ to FOLDER, so
%   that a rerun leaves the recorded results unchanged; ndre_paths("")
%   restores the default. Quadrature rules are always read from the
%   repository.
persistent out
root=fileparts(fileparts(mfilename('fullpath')));
if nargin>0, out=char(output); end
base=out; if isempty(base), base=root; end
p=struct('root',root,'data',fullfile(base,'data'),'figures',fullfile(base,'figures'), ...
    'tables',fullfile(base,'tables'),'quadrature',fullfile(root,'data','quadrature'));
end
