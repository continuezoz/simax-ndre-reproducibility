function plot_bound_history(tau)
%PLOT_BOUND_HISTORY Figure 1: true error, estimate and target on the nodes
% of the accepted space (data of run_bound_history).
if nargin<1, tau=1e-7; end
paths=ndre_paths(); in=fullfile(paths.data,'bounds');
out=paths.figures; if ~exist(out,'dir'), mkdir(out); end
ids=["transport","multistate"]; titles=["(a) single-state, c = 0.9999","(b) multistate, \gamma = 0.9999"];
f=figure('Units','centimeters','Position',[2 2 16 6.8],'Color','w');
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
for i=1:2
    T=readtable(fullfile(in,char("history_"+ids(i)+".csv")));
    T=T(T.t>0,:); ax=nexttile(tl); hold(ax,'on');
    loglog(ax,T.t,T.estimate,'-','Color',[0 0.32 0.62],'LineWidth',1.2);
    loglog(ax,T.t,T.projectionError,'-','Color',[0.80 0.25 0.10],'LineWidth',1.2);
    loglog(ax,T.t,T.error,':','Color',[0.55 0.55 0.55],'LineWidth',1.1);
    loglog(ax,T.t,1e-12+tau*T.referenceNorm,'--','Color',[0.35 0.35 0.35],'LineWidth',0.9);
    set(ax,'XScale','log','YScale','log','FontSize',8,'Box','on','TickLabelInterpreter','tex');
    xlim(ax,[min(T.t) 1]); grid(ax,'on'); ax.GridAlpha=0.15; ax.MinorGridAlpha=0.05;
    xlabel(ax,'t'); title(ax,titles(i),'FontWeight','normal','FontSize',8.5);
    if i==1, ylabel(ax,'Frobenius norm'); end
    if i==1, lg=legend(ax,{'estimate \epsilon_j','projection error','computed output error','target \tau_a + \tau||X||_F'}, ...
            'FontSize',7.5,'Box','off','NumColumns',4);
        lg.Layout.Tile='south'; end
end
exportgraphics(f,fullfile(out,'bound_history.pdf'),'ContentType','vector');
close(f);
end
