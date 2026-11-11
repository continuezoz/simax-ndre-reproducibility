function [omega,w] = load_gauss_legendre_matlab(n)
%LOAD_GAUSS_LEGENDRE_MATLAB Load only a MATLAB-generated quadrature cache.
root=fileparts(mfilename('fullpath'));
folder=fullfile(root,'results','matlab_quadrature_20261007');
if ~exist(folder,'dir'), mkdir(folder); end
cache=fullfile(folder,sprintf('gauss_legendre_%d.mat',n));
if exist(cache,'file')
    data=load(cache,'omega','w','info');
    assert(isfield(data,'info') && strcmp(data.info.generator,'gauss_legendre_matlab'), ...
        'quadrature:provenance','Quadrature cache lacks MATLAB provenance.');
    omega=data.omega; w=data.w;
else
    [omega,w,info]=gauss_legendre_matlab(n);
    save(cache,'omega','w','info');
end
end
