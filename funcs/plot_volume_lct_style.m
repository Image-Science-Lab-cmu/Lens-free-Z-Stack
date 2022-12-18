function plot_volume_lct_style(vol, type, xtl, ytl, ztl)
% plot volume (z, x, y)
D = mat2gray(vol);
colormap gray;
h = vol3d('cdata', D, 'Alpha', D, 'texture', '3D');
alphamap('rampup'); alphamap(.06 .* alphamap);

% set background to black
set(gca,'color','k');

% add red box
box_width = 2;
hold on,
% lines in x
plot3([0 size(vol,2)],                [0 0],                [0 0], 'r-', 'LineWidth', box_width)
plot3([0 size(vol,2)], size(vol, 1) * [1 1],                [0 0], 'r-', 'LineWidth', box_width)
plot3([0 size(vol,2)],                [0 0], size(vol, 3) * [1 1], 'r-', 'LineWidth', box_width)
plot3([0 size(vol,2)], size(vol, 1) * [1 1], size(vol, 3) * [1 1], 'r-', 'LineWidth', box_width)
% lines in y
plot3(               [0 0],                [0 0],[0 size(vol, 3)], 'r-', 'LineWidth', box_width)
plot3(               [0 0], size(vol, 1) * [1 1],[0 size(vol, 3)], 'r-', 'LineWidth', box_width)
plot3(size(vol, 2) * [1 1],                [0 0],[0 size(vol, 3)], 'r-', 'LineWidth', box_width)
plot3(size(vol, 2) * [1 1], size(vol, 1) * [1 1],[0 size(vol, 3)], 'r-', 'LineWidth', box_width)
% lines in z
plot3(               [0 0], [0 size(vol, 1)],                [0 0], 'r-', 'LineWidth', box_width)
plot3(size(vol, 2) * [1 1], [0 size(vol, 1)],                [0 0], 'r-', 'LineWidth', box_width)
plot3(               [0 0], [0 size(vol, 1)], size(vol, 3) * [1 1], 'r-', 'LineWidth', box_width)
plot3(size(vol, 2) * [1 1], [0 size(vol, 1)], size(vol, 3) * [1 1], 'r-', 'LineWidth', box_width)
if nargin < 3
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'ZTick', []);
else
    nxtl = length(xtl);
    nytl = length(ytl);
    nztl = length(ztl);
    xticks(linspace(0, size(vol, 2), nxtl));
    xticklabels(xtl);
    yticks(linspace(0, size(vol, 1), nztl));
    yticklabels(ztl);
    zticks(linspace(0, size(vol, 3), nytl));
    zticklabels(ytl);
end
axis tight; %axis off;
if nargin > 1 && strcmp(type, 'angle_diopter') == 1
    xlabel('x/z'); ylabel('1/z'); zlabel('y/z');
    set(gca, 'YDir', 'reverse');
else
    xlabel('x'); ylabel('z'); zlabel('y');
end

az = 20;
el = 37.5;
view([az el]);

end