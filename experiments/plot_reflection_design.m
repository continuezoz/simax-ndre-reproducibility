function plot_reflection_design
%PLOT_REFLECTION_DESIGN Figure 3: (a) albedo ratio a(x)/a(inf) of the
% multistate reflector from Algorithm 5.1 and from EBA--BDF1 (h=0.01);
% (b) relative half-width of the enclosure of the design length x_p and
% relative errors of EBA--BDF1, measured against x_p of Algorithm 5.1
% (data of run_reflection_design).
paths=ndre_paths(); in=fullfile(paths.data,'reflection_design');
out=paths.figures; if ~exist(out,'dir'), mkdir(out); end
gammas=[0.5 0.9 0.9999]; tfs=[0.5 2 15];
colors=[0 0.32 0.62; 0.80 0.25 0.10; 0.20 0.55 0.25];
f=figure('Units','centimeters','Position',[2 2 16 6.8],'Color','w');
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
% (a) albedo curves
ax=nexttile(tl); hold(ax,'on'); axa=ax;
for i=1:3
    tag=sprintf('g%g',gammas(i));
    A=readtable(fullfile(in,['albedo_' tag '.csv'])); A=A(A.x>0,:);
    B=readtable(fullfile(in,['albedo_bdf1_' tag '.csv'])); B=B(B.x>0,:);
    plot(ax,A.x,A.albedoRatio,'-','Color',colors(i,:),'LineWidth',1.2);
    plot(ax,B.x,B.albedoRatio,'--','Color',colors(i,:),'LineWidth',1.0);
    early=B.x<=0.1+1e-12;
    plot(ax,B.x(early),B.albedoRatio(early),'o','Color',colors(i,:), ...
        'MarkerSize',3,'MarkerFaceColor','w','LineWidth',0.8);
end
for level=[0.5 0.9 0.99]
    plot(ax,[1e-4 max(tfs)],[level level],':','Color',[0.45 0.45 0.45],'LineWidth',0.8);
end
set(ax,'XScale','log','FontSize',8,'Box','on');
xlim(ax,[1e-4 max(tfs)]); ylim(ax,[0 1.02]); grid(ax,'on'); ax.GridAlpha=0.15; ax.MinorGridAlpha=0.05;
xlabel(ax,'rod length x'); ylabel(ax,'a(x) / a(\infty)');
title(ax,'(a) albedo','FontWeight','normal','FontSize',8.5);
% (b) design lengths
ax=nexttile(tl); hold(ax,'on');
for i=1:3
    D=readtable(fullfile(in,sprintf('design_g%g.csv',gammas(i))));
    plot(ax,D.p,(D.xHi-D.xLo)./(2*D.x),'-','Color',colors(i,:),'LineWidth',1.2);
    plot(ax,D.p,abs(D.xBDF1_h100-D.x)./D.x,'--','Color',colors(i,:),'LineWidth',1.0);
    plot(ax,D.p,abs(D.xBDF1_h1000-D.x)./D.x,':','Color',colors(i,:),'LineWidth',1.2);
end
set(ax,'YScale','log','FontSize',8,'Box','on');
xlim(ax,[0.3 1]); grid(ax,'on'); ax.GridAlpha=0.15; ax.MinorGridAlpha=0.05;
xlabel(ax,'fraction p'); ylabel(ax,'relative error / half-width of x_p');
title(ax,'(b) design length x_p','FontWeight','normal','FontSize',8.5);
% legends: colours (gamma) in panel (a), line styles (methods) in panel (b)
gray=[0.4 0.4 0.4];
h=gobjects(3,1);
for i=1:3, h(i)=plot(ax,NaN,NaN,'-','Color',gray,'LineWidth',1.2); end
set(h(2),'LineStyle','--','LineWidth',1.0); set(h(3),'LineStyle',':');
legend(ax,h,{'Algorithm 5.1 (enclosure)','EBA-BDF1, h = 0.01','EBA-BDF1, h = 0.001'}, ...
    'FontSize',7,'Box','off','Location','west');
g=gobjects(3,1);
for i=1:3, g(i)=plot(axa,NaN,NaN,'-','Color',colors(i,:),'LineWidth',2.5); end
legend(axa,g,{'\gamma = 0.5','\gamma = 0.9','\gamma = 0.9999'}, ...
    'FontSize',7,'Box','off','Location','northwest');
exportgraphics(f,fullfile(out,'reflection_design.pdf'),'ContentType','vector');
close(f);
end
