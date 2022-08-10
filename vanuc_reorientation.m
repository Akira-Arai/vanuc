function fileout = vanuc_reorientation(filein)
% Reorientation
% 
% Input:
% filein (char): Input file name
% 
% Return:
% fileout (char): Output file name
% 
% Output data:
% 'ro*.nii': Reorientated image
% 'fro*.nii' (option): Reorientated and flipped image
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Keep or change
% ----------------------------------------------------------------
V = spm_vol(filein);
kmat = zeros(4);
kmat(4, 4) = 1;
[a, b] = max(abs(V.mat(1:3, 1:3)));
kmat(b(1), 1) = 1;
kmat(b(2), 2) = 1;
kmat(b(3), 3) = 1;
kmat = kmat .* sign(V.mat);
kmat(4, 1:3) = -1 * (V.dim + 1) .* min(kmat(1:3, 1:3));
kmat = kmat';
VV(1:2) = V;
VV(1).mat = V.mat * kmat;
VV(1).dim = V.dim * abs(kmat(1:3, 1:3));
spm_reslice(VV, struct('interp', 1, 'mask', 0, 'mean', 0, 'which', 1, 'prefix', 'ro'));
file = ['ro' filein];
fim = vanuc_disp(file);
figure(fim);
c1 = uicontrol;
c1.Position = [561 326 140 25];
c1.String = 'Confirm directions';
c1.FontSize = 12;
c1.Callback = @Confirm;
c2 = uicontrol;
c2.Position = [561 306 140 20];
c2.String = 'Rotate/Flip';
c2.FontSize = 12;
c2.Callback = @RotateFlip;
	function Confirm(src, event)
		fileout = file;
		uiresume;
		close(fim);
	end
	function RotateFlip(src, event)
		close(fim);
		fileout = rotateHF2(file);
	end
uiwait;
end

function [fileout, fileorigin] = rotateHF2(filein, fileorigin)
% Rotate HF
% ----------------------------------------------------------------

if ~exist('fileorigin', 'var')
	fileorigin = filein;
end
V = spm_vol(filein);
Y = spm_read_vols(V);
fileout = ['f' fileorigin];
V.fname = fileout;
fim = vanuc_disp(filein);
figure(fim);
ca1 = uicontrol;
ca1.Position = [136 331 80 20];
ca1.String = 'This is H';
ca1.FontSize = 14;
ca1.Callback = @AH;
ca2 = uicontrol;
ca2.Position = [136 1 80 20];
ca2.String = 'This is H';
ca2.FontSize = 14;
ca2.Callback = @PH;
ca3 = uicontrol;
ca3.Position = [1 166 80 20];
ca3.String = 'This is H';
ca3.FontSize = 14;
ca3.Callback = @RH;
ca4 = uicontrol;
ca4.Position = [271 166 80 20];
ca4.String = 'This is H';
ca4.FontSize = 14;
ca4.Callback = @LH;
cs1 = uicontrol;
cs1.Position = [486 331 80 20];
cs1.String = 'This is H';
cs1.FontSize = 14;
cs1.Callback = @HH;
cs2 = uicontrol;
cs2.Position = [486 1 80 20];
cs2.String = 'This is H';
cs2.FontSize = 14;
cs2.Callback = @FH;
cs3 = uicontrol;
cs3.Position = [351 166 80 20];
cs3.String = 'This is H';
cs3.FontSize = 14;
cs3.Callback = @AH;
cs4 = uicontrol;
cs4.Position = [621 166 80 20];
cs4.String = 'This is H';
cs4.FontSize = 14;
cs4.Callback = @PH;
text(361, 41, 'Which direction is HEAD (upward)?', 'Color', 'w', 'FontSize', 14);
	function AH(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([1 3 2 4 6 5 7 9 8 10 12 11]));
		V.dim = V.dim([1 3 2]);
		Y = permute(Y, [1 3 2]);
		Y = Y(V.dim(1):-1:1, :, :);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function PH(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([1 3 2 4 6 5 7 9 8 10 12 11]));
		V.dim = V.dim([1 3 2]);
		Y = permute(Y, [1 3 2]);
		Y = Y(:, :, V.dim(3):-1:1);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function RH(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([3 2 1 6 5 4 9 8 7 12 11 10]));
		V.dim = V.dim([3 2 1]);
		Y = permute(Y, [3 2 1]);
		Y = Y(V.dim(1):-1:1, :, :);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function LH(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([3 2 1 6 5 4 9 8 7 12 11 10]));
		V.dim = V.dim([3 2 1]);
		Y = permute(Y, [3 2 1]);
		Y = Y(:, :, V.dim(3):-1:1);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function HH(src, event)
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function FH(src, event)
		Y = Y(V.dim(1):-1:1, :, V.dim(3):-1:1);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = rotateAP2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
uiwait;
end

function [fileout, fileorigin] = rotateAP2(filein, fileorigin)
% Rotate AP
% ----------------------------------------------------------------

if ~exist('fileorigin', 'var')
	fileorigin = filein;
end
V = spm_vol(filein);
Y = spm_read_vols(V);
fileout = ['f' fileorigin];
V.fname = fileout;
fim = vanuc_disp(filein);
figure(fim);
ca1 = uicontrol;
ca1.Position = [136 331 80 20];
ca1.String = 'This is A';
ca1.FontSize = 14;
ca1.Callback = @AA;
ca2 = uicontrol;
ca2.Position = [136 1 80 20];
ca2.String = 'This is A';
ca2.FontSize = 14;
ca2.Callback = @PA;
ca3 = uicontrol;
ca3.Position = [1 166 80 20];
ca3.String = 'This is A';
ca3.FontSize = 14;
ca3.Callback = @RA;
ca4 = uicontrol;
ca4.Position = [271 166 80 20];
ca4.String = 'This is A';
ca4.FontSize = 14;
ca4.Callback = @LA;
cs3 = uicontrol;
cs3.Position = [351 166 80 20];
cs3.String = 'This is A';
cs3.FontSize = 14;
cs3.Callback = @AA;
cs4 = uicontrol;
cs4.Position = [621 166 80 20];
cs4.String = 'This is A';
cs4.FontSize = 14;
cs4.Callback = @PA;
text(361, 41, 'Which direction is Anterior (forward)?', 'Color', 'w', 'FontSize', 14);
	function AA(src, event)
		spm_write_vol(V, Y);
		[fileout, fileorigin] = flipRL2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function PA(src, event)
		Y = Y(V.dim(1):-1:1, V.dim(2):-1:1, :);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = flipRL2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function RA(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([2 1 3 5 4 6 8 7 9 11 10 12]));
		V.dim = V.dim([2 1 3]);
		Y = permute(Y, [2 1 3]);
		Y = Y(V.dim(1):-1:1, :, :);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = flipRL2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
	function LA(src, event)
		Vec = spm_imatrix(V.mat);
		V.mat = spm_matrix(Vec([2 1 3 5 4 6 8 7 9 11 10 12]));
		V.dim = V.dim([2 1 3]);
		Y = permute(Y, [2 1 3]);
		Y = Y(:, V.dim(2):-1:1, :);
		spm_write_vol(V, Y);
		[fileout, fileorigin] = flipRL2(fileout, fileorigin);
		uiresume;
		close(fim);
	end
uiwait;
end

function [fileout, fileorigin] = flipRL2(filein, fileorigin)
% Flip
% ----------------------------------------------------------------

if ~exist('fileorigin', 'var')
	fileorigin = filein;
end
V = spm_vol(filein);
Y = spm_read_vols(V);
fileout = ['f' fileorigin];
V.fname = fileout;
fim = vanuc_disp(filein);
figure(fim);
c1 = uicontrol;
c1.Position = [211 326 140 25];
c1.String = 'Confirm directions';
c1.FontSize = 12;
c1.Callback = @RR;
c2 = uicontrol;
c2.Position = [211 301 140 25];
c2.String = 'Flip R/L';
c2.FontSize = 14;
c2.Callback = @LR;
text(11, 41, 'Flip Right & Left?', 'Color', 'w', 'FontSize', 14);
	function RR(src, event)
		spm_write_vol(V, Y);
		uiresume;
		close(fim);
	end
	function LR(src, event)
		Y = Y(V.dim(1):-1:1, :, :);
		spm_write_vol(V, Y);
		uiresume;
		close(fim);
	end
uiwait;
end