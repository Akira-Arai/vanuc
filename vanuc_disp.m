function fim = vanuc_disp(Y, c, Vx, Vy, Vz)
% Display axial & sagittal slice of 3D image
% 
% Input:
% Y - (char): NIfTI file of 3D image
%     (3D double): 3x3x3 matrix of image
% c - 1: Color
%     0 (default): Grayscale
% Vx (double): Width of pixel
% Vy (double): Height of pixel
% Vz (double): Thickness of slice
% 
% Return:
% fim (char): Figure window name
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Reading file
% ----------------------------------------------------------------
if ischar(Y) && strcmp(Y(numel(Y) - 3 : numel(Y)), '.nii')
	V = spm_vol(Y);
	Vec = spm_imatrix(V.mat);
	Y = spm_read_vols(V);
	Dim = V.dim;
	Size = Dim .* abs(Vec(7 : 9));
elseif isnumeric(Y) && numel(size(Y))==3 && (exist('Vx', 'var') && exist('Vy', 'var') && exist('Vz', 'var'))
	Dim = size(Y);
	Size = Dim .* [Vx Vy Vz];
else
	return;
end
Y(find(isinf(Y))) = max(Y(find(isfinite(Y))));
Y(find(isnan(Y))) = 0;

% Dispray data (size, window, color)
% ----------------------------------------------------------------
DispSize = ceil(Size / max(Size) * 350);
DispSize(find(DispSize<1)) = 1;
DispSize(find(DispSize>350)) = 350;
Yposition = floor(175 - DispSize * 0.5) + 1;
YA = zeros(350, 350);
YA(Yposition(1) : Yposition(1) + DispSize(1) - 1, Yposition(2) : Yposition(2) + DispSize(2) - 1) = imresize(Y(:, :, ceil(0.5 * Dim(3))), DispSize(1 : 2));
YS = zeros(350, 350);
YS(Yposition(2) : Yposition(2) + DispSize(2) - 1, Yposition(3) : Yposition(3) + DispSize(3) - 1) = imresize(reshape(Y(ceil(0.5 * Dim(1)), :, :), Dim(2 : 3)), DispSize(2 : 3));
YAS = [YA(350:-1:1, 350:-1:1)' YS(350:-1:1, 350:-1:1)'];
sortY = sort(YAS(:));
Low = sortY(1225);
High = sortY(243775);
if exist('c', 'var') && isnumeric(c) && c
	YAS = (YAS - Low) / (High - Low) * 255;
	rYAS = zeros(350, 700);
	rYAS(find(YAS < 36)) = 3 + YAS(find(YAS < 36)) * 69 / 36;
	rYAS(find((YAS >= 36) .* (YAS < 72))) = 72 - (YAS(find((YAS >= 36) .* (YAS < 72))) - 36) * 72 / 36;
	rYAS(find((YAS >= 136) .* (YAS < 144))) = (YAS(find((YAS >= 136) .* (YAS < 144))) - 136) * 128 / 8;
	rYAS(find((YAS >= 144) .* (YAS < 160))) = 128 + (YAS(find((YAS >= 144) .* (YAS < 160))) - 144) * 127 / 16;
	rYAS(find(YAS >= 160)) = 255;
	gYAS = zeros(350, 700);
	gYAS(find((YAS >= 36) .* (YAS < 72))) = (YAS(find((YAS >= 36) .* (YAS < 72))) - 36) * 72 / 36;
	gYAS(find((YAS >= 72) .* (YAS < 96))) = 72 + (YAS(find((YAS >= 72) .* (YAS < 96))) - 72) * 64 / 24;
	gYAS(find((YAS >= 96) .* (YAS < 136))) = 136 + (YAS(find((YAS >= 96) .* (YAS < 136))) - 96) * 119 / 40;
	gYAS(find((YAS >= 136) .* (YAS < 160))) = 255;
	gYAS(find((YAS >= 160) .* (YAS < 192))) = 255 - (YAS(find((YAS >= 160) .* (YAS < 192))) - 160) * 63 / 32;
	gYAS(find(YAS >= 192)) = 192 - (YAS(find(YAS >= 192)) - 192) * 189 / 63;
	bYAS = zeros(350, 700);
	bYAS(find(YAS < 36)) = YAS(find(YAS < 36)) * 108 / 36;
	bYAS(find((YAS >= 36) .* (YAS < 72))) = 108 + (YAS(find((YAS >= 36) .* (YAS < 72))) - 36) * 84 / 36;
	bYAS(find((YAS >= 72) .* (YAS < 96))) = 192 - (YAS(find((YAS >= 72) .* (YAS < 96))) - 72) * 22 / 24;
	bYAS(find((YAS >= 96) .* (YAS < 136))) = 170 - (YAS(find((YAS >= 96) .* (YAS < 136))) - 96) * 170 / 40;
	YAS = rYAS / 255;
	YAS(:, :, 2) = gYAS / 255;
	YAS(:, :, 3) = bYAS / 255;
	fim = figure;
	imshow(YAS, 'Border', 'tight');
	fim.Position=[300 200 700 350];
	fim.Name = 'Axial/Sagittal';
	fim.NumberTitle = 'off';
	fim.MenuBar='none';
	fim.ToolBar='none';
else
	fim = figure;
	imshow(YAS, [Low High], 'Border', 'tight');
	fim.Position=[300 200 700 350];
	fim.Name = 'Axial/Sagittal';
	fim.NumberTitle = 'off';
	fim.MenuBar='none';
	fim.ToolBar='none';
end
text(171, 11, 'A', 'Color', [1 1 0], 'FontSize', 16);
text(171, 341, 'P', 'Color', [1 1 0], 'FontSize', 16);
text(3, 171, 'R', 'Color', [1 1 0], 'FontSize', 16);
text(333, 171, 'L', 'Color', [1 1 0], 'FontSize', 16);
text(521, 11, 'H', 'Color', [1 1 0], 'FontSize', 16);
text(521, 341, 'F', 'Color', [1 1 0], 'FontSize', 16);
text(353, 171, 'A', 'Color', [1 1 0], 'FontSize', 16);
text(683, 171, 'P', 'Color', [1 1 0], 'FontSize', 16);
end