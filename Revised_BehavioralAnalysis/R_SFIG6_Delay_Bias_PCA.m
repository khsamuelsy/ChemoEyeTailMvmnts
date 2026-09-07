close all;
%% R_SFIG6_Delay_Bias_PCA.m
% Run AFTER:
%   R_FIG3G_Delay_biasFIG4H_PC1_PC2_MixedModel_Angle.m
%
% Expected variables already in workspace:
%   all_x        : N x 5 matrix of temporal predictors
%   all_group3   : N x 1 group labels (1=blank, 2=app, 3=av)
%   coeff        : 5 x 5 PCA loading matrix from pca(all_x)
%   score        : N x 5 PCA scores
%   explained    : 5 x 1 percent variance explained
%
% Optional variables:
%   PC1, PC2     : if not present, they are taken from score(:,1:2)

%% -------------------- checks --------------------
if ~exist('all_x','var') || isempty(all_x)
    error('all_x not found. Run R_FIG3G.m first.');
end
if ~exist('all_group3','var') || isempty(all_group3)
    error('all_group3 not found. R_FIG3G.m first.');
end
if ~exist('coeff','var') || isempty(coeff) || ...
   ~exist('score','var') || isempty(score) || ...
   ~exist('explained','var') || isempty(explained)
    error('PCA outputs coeff/score/explained not found. Run R_FIG3G.m first.');
end


%% -------------------- labels and colors --------------------
varLabels = {'x1','x2','x3','x4','x5'};
blankColor = [0.55 0.55 0.55];
appeColor  = [119 136 172]./255;
averColor  = [238 129 114]./255;

g1 = all_group3 == 1;
g2 = all_group3 == 2;
g3 = all_group3 == 3;

%% -------------------- pooled correlation matrix --------------------
showSigOnly = true;     % true = display only significant entries
alpha       = 0.05 / 10;     % uncorrected threshold for display

[R, P] = corr(all_x, 'rows', 'pairwise');

if showSigOnly
    Rplot = R;
    Rplot(P >= alpha) = NaN;   % non-significant entries become NaN
    figName  = 'Supp: Temporal predictor correlation (sig only)';
else
    Rplot = R;
    figName  = 'Supp: Temporal predictor correlation';
end

f1 = figure('Name', figName, ...
            'Position',[100 100 440 380], 'Color','w');

hImg = imagesc(Rplot, [-1 1]);
axis square;

% Make NaNs transparent so white axes background shows through
alphaMask = ones(size(Rplot));
alphaMask(isnan(Rplot)) = 0;
set(hImg, 'AlphaData', alphaMask);

% Harmonic diverging colormap
n = 256;
negColor = [76 129 139]./255;
midColor = [246 242 235]./255;
posColor = [214 110 89]./255;

n1 = floor(n/2);
n2 = n - n1;

cmap1 = [linspace(negColor(1), midColor(1), n1)', ...
         linspace(negColor(2), midColor(2), n1)', ...
         linspace(negColor(3), midColor(3), n1)'];

cmap2 = [linspace(midColor(1), posColor(1), n2)', ...
         linspace(midColor(2), posColor(2), n2)', ...
         linspace(midColor(3), posColor(3), n2)'];

cmap = [cmap1; cmap2];
colormap(gca, cmap);

ax = gca;
ax.Color = [1 1 1];   % white background for transparent NaN cells
set(ax, 'Layer', 'top');

cb = colorbar;
cb.Label.String = 'Pearson r';

xticks(1:5);
yticks(1:5);
xticklabels(varLabels);
yticklabels(varLabels);
set(gca, 'TickDir','out', 'Box','off', 'FontSize',10, 'LineWidth',1);

for i = 1:5
    for j = 1:5
        if ~isnan(Rplot(i,j))
            if abs(Rplot(i,j)) > 0.5
                txtcol = [1 1 1];
            else
                txtcol = [0.15 0.15 0.15];
            end
            text(j, i, sprintf('%.2f', Rplot(i,j)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color', txtcol, ...
                'FontSize', 8);
        end
    end
end

%% -------------------- PC1 and PC2 loading bars --------------------
f2 = figure('Name','Supp: PC loadings', ...
            'Position',[580 100 520 360], 'Color','w');

load12 = coeff(:,1:2);
b = bar(load12, 'grouped', 'LineStyle','none');
b(1).FaceColor = [0.28 0.28 0.28];
b(2).FaceColor = [0.70 0.70 0.70];

yline(0, 'k-');
xticks(1:5);
xticklabels(varLabels);
ylabel('Loading');
title(sprintf('PC1 (%.1f%%) and PC2 (%.1f%%)', explained(1), explained(2)), ...
      'FontWeight','normal');
legend({'PC1','PC2'}, 'Location','northoutside', 'Orientation','horizontal');
set(gca,'TickDir','out','Box','off','FontSize',10);

%% -------------------- scree plot --------------------
f3 = figure('Name','Supp: Scree plot', ...
            'Position',[1140 100 380 360], 'Color','w');

plot(1:numel(explained), explained, '-o', ...
    'Color', [0.22 0.22 0.22], ...
    'MarkerFaceColor', [0.22 0.22 0.22], ...
    'LineWidth', 1.5, ...
    'MarkerSize', 6);
hold on;

plot(1:numel(explained), cumsum(explained), '--s', ...
    'Color', [0.55 0.55 0.55], ...
    'MarkerFaceColor', [0.55 0.55 0.55], ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5);

xlim([0.8 numel(explained)+0.2]);
ylim([0 100]);
xticks(1:numel(explained));
xlabel('Principal component');
ylabel('Variance explained (%)');
legend({'Individual','Cumulative'}, 'Location','eastoutside');
set(gca,'TickDir','out','Box','off','FontSize',10);