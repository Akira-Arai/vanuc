function [Sx, Sy, Sz] = vanuc_phantom(xyLlimit, xyUlimit, xyprecision, zLlimit, zUlimit, zprecision, savepara)
% Estimation of sigma of PSF from phantom data
% 
% (x-sigma equals y-sigma)
% 
% Input:
% xyLlimit (double): Lower limit of xy-sigma
% xyUlimit (double): Upper limit of xy-sigma
% xyprecision (double): Precision of xy-sigma
% zLlimit (double): Lower limit of z-sigma
% zUlimit (double): Upper limit of z-sigma
% zprecision (double): Precision of z-sigma
% savepara - 1(default): Save parameters
% 
% Return:
% Sx (double): x-sigma
% Sy (double): y-sigma
% Sz (double): z-sigma
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

if ~exist('xyLlimit');
	xyLlimit = 0;
	xyUlimit = 6;
	xyprecision = 0.001;
	zLlimit = 0;
	zUlimit = 6;
	zprecision = 0.001;
end
Sx = 0;
Sy = 0;
Sz = 0;

% Read phantom image
% ----------------------------------------------------------------
file = vanuc_dcmconvert('single');
if numel(file)==0
	return
end
[file, Vh, Vw, Vt] = vanuc_crop(file{1});

% Create digital phantom
% ----------------------------------------------------------------
PathH = pwd;
PathV = which('vanuc.m');
PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
cd([PathV 'phantom']);
[fileP, pathP] = uigetfile('*.nii', 'Select Phantom Model');
if ~ischar(fileP)
	return
end
cd(PathH);
if ~strcmp([pathP, fileP], [pwd '\DigitalPhantom.nii'])
	copyfile([pathP, fileP], 'DigitalPhantom.nii');
end
V = spm_vol('DigitalPhantom.nii');
Y = spm_read_vols(V);
ValueList = flip(unique(Y));
if numel(ValueList)>256
	Y = 255 * (Y - min(Y)) / (max(Y) - min(Y));
	Y = double(uint8(Y));
end
V.dt = [16 1];
spm_write_vol(V, Y);
for k = 1 : numel(ValueList)
	Y2 = zeros(size(Y));
	Y2(find(Y == ValueList(k))) = 1;
	V2 = V;
	V2.fname = ['c' num2str(k) 'Phantom.nii'];
	spm_write_vol(V2, Y2);
	FileList{k, 1} = V2.fname;
end
clear Y Y2 V2

% Resize
% ----------------------------------------------------------------
dimat = spm_imatrix(V.mat);
px = min(abs(dimat(7:9)));
copyfile(file, ['copy_' file]);
Vpet = spm_vol(file);
G = spm_read_vols(Vpet);
dimat = spm_imatrix(Vpet.mat);
Pratio = ceil(abs(dimat(7:9)/px));
dimat(7:9) = dimat(7:9)./Pratio;
Vpet.mat = spm_matrix(dimat);
Size = size(G);
G2 = zeros(Size.* Pratio);
Vpet.dim = size(G2);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			G2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1)) = G;
		end
	end
end
spm_write_vol(Vpet, G2);
clear G2

% Coregister digital phantom to the observed image
% ----------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {file};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'DigitalPhantom.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = FileList;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear matlabbatch

% Create true distribution image
% ----------------------------------------------------------------
delete(file);
copyfile(['copy_' file], file);
delete(['copy_' file]);
V = spm_vol(file);
Size = V.dim;
R = zeros(Size);
M = zeros([Size numel(FileList)]);
for m = 1 : numel(FileList);
	delete(FileList{m});
	FileList{m} = ['r' FileList{m}];
	Vr = spm_vol(FileList{m});
	Pratio = Vr.dim./V.dim;
	Y2 = spm_read_vols(Vr);
	Y2(find(isfinite(Y2)==0)) = 0;
	for i=1:Pratio(1)
		for j=1:Pratio(2)
			for k=1:Pratio(3)
				M(:, :, :, m) = M(:, :, :, m) + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
			end
		end
	end
	M(:, :, :, m) = M(:, :, :, m) / Pratio(1) / Pratio(2) / Pratio(3);
	R = R + M(:, :, :, m) * ValueList(m);
end
delete('rDigitalPhantom.nii');

% Count Normalise
% ----------------------------------------------------------------
G = G * mean(R(:)) / mean(G(:));
save('Phantom.mat', 'G', 'R', 'M');

% Estimate sigma
% ----------------------------------------------------------------
[Sx, Sz] = vanuc_measureresol(G, R, xyLlimit, xyUlimit, xyprecision, zLlimit, zUlimit, zprecision);
Sy = Sx;
if ~exist('savepara')
	savepara = 1;
end
if savepara
	Fh = Sx * Vh * 2 * sqrt(2 * log(2));
	Fw = Sy * Vw * 2 * sqrt(2 * log(2));
	Ft = Sz * Vt * 2 * sqrt(2 * log(2));
	Input = inputdlg({'Protocol name', 'FWHM_x of PSF (mm)   (0<fx<100)', 'FWHM_y of PSF(mm)   (0<fy<100)', 'FWHM_z of PSF(mm)   (0<fz<100)', 'Distance_x from center of local region (vx)   (dx<15)', 'Distance_y from center of local region (vx)   (dy<15)', 'Distance_z from center of local region (vx)   (dz<15)', 'Regularization parameter'}, 'Input Parameters', [1 30], {'(no name)', num2str(Fh), num2str(Fw), num2str(Ft), num2str(ceil(20/Vh)), num2str(ceil(20/Vw)), num2str(ceil(20/Vt)), '1'});
	if numel(Input)==0
		return
	end
	Pname = Input{1};
	Fh = str2num(Input{2});
	Fw = str2num(Input{3});
	Ft = str2num(Input{4});
	Lh = str2num(Input{5});
	Lw = str2num(Input{6});
	Lt = str2num(Input{7});
	Lam = str2num(Input{8});
	PathV = which('vanuc.m');
	PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
	if ~isdir([PathV 'parameters'])
		mkdir([PathV 'parameters']);
	end
	save([PathV 'parameters\' Pname '.mat'], 'Pname', 'Fh', 'Fw', 'Ft', 'Vh', 'Vw', 'Vt', 'Lh', 'Lw', 'Lt', 'Lam');
end
end