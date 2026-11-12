function plot_floor_traces(tau)
%PLOT_FLOOR_TRACES Figure 2: estimated relative accuracy tau*max_j
% {eps_j/(tau_a+tau*||Y(sigma_j)||_F)} of Algorithm 5.1 against the dimension
% for several n (data of run_floor_traces), with the levels u*||A||.
if nargin<1, tau=1e-7; end
paths=ndre_paths(); in=fullfile(paths.data,'floor');
out=paths.figures; if ~exist(out,'dir'), mkdir(out); end
S=readtable(fullfile(in,'floor_summary.csv'));
colors=[0.00 0.32 0.62; 0.10 0.55 0.35; 0.85 0.50 0.05; 0.80 0.25 0.10];
f=figure('Units','centimeters','Position',[2 2 9.5 7.2],'Color','w');
ax=axes(f); hold(ax,'on'); h=gobjects(0); names=strings(0);
for i=1:height(S)
    T=readtable(fullfile(in,sprintf('trace_n%d.csv',S.n(i))));
    c=colors(min(i,size(colors,1)),:);
    h(end+1)=semilogy(ax,T.k,tau*T.ratio,'-o','Color',c,'MarkerSize',3,'MarkerFaceColor',c,'LineWidth',1.1); %#ok<AGROW>
    yline(ax,S.uNormA(i),':','Color',c,'LineWidth',1.1);
    names(end+1)=sprintf('n = %d',S.n(i)); %#ok<AGROW>
end
h(end+1)=yline(ax,tau,'--','Color',[0.35 0.35 0.35],'LineWidth',0.9); names(end+1)="\tau = 10^{-7}";
set(ax,'YScale','log','FontSize',8,'Box','on'); grid(ax,'on'); ax.GridAlpha=0.15;
xlabel(ax,'dimension k_m'); ylabel(ax,'\tau \cdot max_j \{\epsilon_j / (\tau_a + \tau||Y(\sigma_j)||_F)\}');
legend(ax,h,names,'Location','northeast','FontSize',7.5,'Box','off');
exportgraphics(f,fullfile(out,'floor_traces.pdf'),'ContentType','vector');
close(f);
end
