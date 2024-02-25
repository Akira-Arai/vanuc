function Success = vanuc_pre(varargin)
% Prepare for partial volume effect correction
% 
% This function consists of the following steps:
%    1. Selection of input images
%    2. Parameters setting
%    3. Reorientation and cropping in preparation
%    4. Segmentation and normalization of 3D T1WI
%    5. Coregistration of MRI to PET/SPECT
% 
% Return:
% Success - 0: Program did not execute completely.
%           1: Program executed successfully.
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Image file reading and parameters setting
% ----------------------------------------------------------------
Success = 0;
[Vx, Vy, Vz, pathH] = vanuc_imagepre('trim_PET.nii', 'PET or SPECT');
if Vx==0 || Vy==0 || Vz==0
	return
end
V = spm_vol('trim_PET.nii');
Y = spm_read_vols(V);
V.dt = [16 1];
spm_write_vol(V, Y);
clear V Y
fim = vanuc_disp('trim_PET.nii', 1);
fim.Name = 'Cropped PET or SPECT image';
Apara = questdlg({'Which do you use for PET parameters?', 'Preset (recommended) or Image-based estimation'}, 'Option for parameters setting', 'Preset', 'Image-based', 'Preset');
if strcmp(Apara, 'Preset')
% ----------------------------------------------------------------
% Read the last set of parameters
	PathV = which('vanuc.m');
	PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
	while 1
		if isfile([PathV 'temp\parameters.mat'])
			load([PathV 'temp\parameters.mat']);
			if exist('Pname') && ischar(Pname) && exist('Fh') && isnumeric(Fh) && Fh>0 && Fh<100 && exist('Fw') && isnumeric(Fw) && Fw>0 && Fw<100 && exist('Ft') && isnumeric(Ft) && Ft>0 && Ft<100 && exist('Vh') && isnumeric(Vh) && Vh>0 && Vh<100 && exist('Vw') && isnumeric(Vw) && Vw>0 && Vw<100 && exist('Vt') && isnumeric(Vt) && Vt>0 && Vt<100 && exist('Lh') && isnumeric(Lh) && Lh>=0 && Lh<15 && ~mod(Lh, 1) && exist('Lw') && isnumeric(Lw) && Lw>=0 && Lw<15 && ~mod(Lw, 1) && exist('Lt') && isnumeric(Lt) && Lt>=0 && Lt<15 && ~mod(Lt, 1) && exist('Lam') && isnumeric(Lam) && Lam>=0
				break;
			end
		end
		Pname = '(default)'
		Fh = 8;
		Fw = 8;
		Ft = 8;
		Vh = Vx;
		Vw = Vy;
		Vt = Vz;
		Lh = ceil(20 / Vh);
		Lw = ceil(20 / Vw);
		Lt = ceil(20 / Vt);
		Lam = 1;
		if ~isdir([PathV 'temp'])
			mkdir([PathV 'temp']);
		end
		save([PathV 'temp\parameters.mat'], 'Pname', 'Fh', 'Fw', 'Ft', 'Vh', 'Vw', 'Vt', 'Lh', 'Lw', 'Lt', 'Lam');
	end
	
% Set parameters
	while 1
		if round(Vh - Vx, 5)==0 && round(Vw - Vy, 5)==0 && round(Vt - Vz, 5)==0
			Sh = Fh / Vh * 0.5 / sqrt(2 * log(2));
			Sw = Fw / Vw * 0.5 / sqrt(2 * log(2));
			St = Ft / Vt * 0.5 / sqrt(2 * log(2));
			answer = questdlg({['Protocol name:                  ' Pname], ['Voxel size(x,y,z):                ' num2str(Vh) 'mm x ' num2str(Vw) 'mm x ' num2str(Vt) 'mm'], ['FWHM(x,y,z) of PSF:           ' num2str(Fh) 'mm x ' num2str(Fw) 'mm x ' num2str(Ft) 'mm'], ['Sigma(x,y,z) of PSF:           ' num2str(0.001 * round(1000 * Sh)) 'vx x ' num2str(0.001 * round(1000 * Sw)) 'vx x ' num2str(0.001 * round(1000 * St)) 'vx'], ['Local region size(x,y,z):       ' num2str(Lh * 2 + 1) 'vx x ' num2str(Lw * 2 + 1) 'vx x ' num2str(Lt * 2 + 1) 'vx'], ['Regularization parameter:    ' num2str(Lam)]}, 'Confirmation of PET parameters', 'OK', 'Other Protocol', 'Dilect Input', 'OK');
			if strcmp(answer, 'OK')
				break
			elseif strcmp(answer, 'Other Protocol')
				Plist = dir([PathV 'parameters\*.mat']);
				if numel(Plist)
					for k = 1 : numel(Plist)
						Plist(k).names = Plist(k).name(1 : numel(Plist(k).name) - 4);
					end
					Pselect = listdlg('SelectionMode', 'single', 'ListString', {Plist.names});
					if numel(Pselect)==0
						return
					end
					load([PathV 'parameters\' Plist(Pselect).name]);
					clear Plist
					if exist('Pname') && ischar(Pname) && exist('Fh') && isnumeric(Fh) && Fh>0 && Fh<100 && exist('Fw') && isnumeric(Fw) && Fw>0 && Fw<100 && exist('Ft') && isnumeric(Ft) && Ft>0 && Ft<100 && exist('Vh') && isnumeric(Vh) && Vh>0 && Vh<100 && exist('Vw') && isnumeric(Vw) && Vw>0 && Vw<100 && exist('Vt') && isnumeric(Vt) && Vt>0 && Vt<100 && exist('Lh') && isnumeric(Lh) && Lh>=0 && Lh<20 && ~mod(Lh, 1) && exist('Lw') && isnumeric(Lw) && Lw>=0 && Lw<20 && ~mod(Lw, 1) && exist('Lt') && isnumeric(Lt) && Lt>=0 && Lt<20 && ~mod(Lt, 1) && exist('Lam') && isnumeric(Lam) && Lam>=0
						if round(Vh - Vx, 5)==0 && round(Vw - Vy, 5)==0 && round(Vt - Vz, 5)==0
							continue
						end
						warning('The size does not match.');
					else
						warning('Some parameters are incorrect.');
					end
				else
					warning('No protocol registered.');
					clear Plist
				end
			elseif ~strcmp(answer, 'Dilect Input')
				return
			end
		end
		Input = inputdlg({'FWHM_x of PSF (mm)   (0<fx<100)', 'FWHM_y of PSF(mm)   (0<fy<100)', 'FWHM_z of PSF(mm)   (0<fz<100)', 'Distance_x from center of local region (vx)   (dx=0,1,2,3,4)', 'Distance_y from center of local region (vx)   (dy=0,1,2,3,4)', 'Distance_z from center of local region (vx)   (dz=0,1,2,3,4)', 'Regularization parameter'}, 'Input Parameters', [1 30], {num2str(Fh), num2str(Fw), num2str(Ft), num2str(Lh), num2str(Lw), num2str(Lt), num2str(Lam)});
		if numel(Input)==0
			return
		end
		if ~exist('Pname', 'var') || Fh ~= str2num(Input{1}) || Fw ~= str2num(Input{2}) || Ft ~= str2num(Input{3}) || Lh ~= str2num(Input{4}) || Lw ~= str2num(Input{5}) || Lt ~= str2num(Input{6}) || Lam ~= str2num(Input{7})
			Pname = '(no name)';
			Fh = str2num(Input{1});
			Fw = str2num(Input{2});
			Ft = str2num(Input{3});
			Lh = str2num(Input{4});
			Lw = str2num(Input{5});
			Lt = str2num(Input{6});
			Lam = str2num(Input{7});
		end
		clear Input
		Vh = Vx;
		Vw = Vy;
		Vt = Vz;
		answer = questdlg('Save parameters as new protocol?', 'Question', 'Yes', 'No', 'No');
		if strcmp(answer, 'Yes')
			Input = inputdlg('Protocol name:', 'Input', [1 50], {'(no name)'});
			if numel(Input)==0
				return
			end
			Pname = Input{1};
			clear Input
			if ~isdir([PathV 'parameters'])
				mkdir([PathV 'parameters']);
			end
			save([PathV 'parameters\' Pname '.mat'], 'Pname', 'Fh', 'Fw', 'Ft', 'Vh', 'Vw', 'Vt', 'Lh', 'Lw', 'Lt', 'Lam');
		elseif ~strcmp(answer, 'No')
			return
		end
	end
	
	if ~isdir([PathV 'temp'])
		mkdir([PathV 'temp']);
	end
	save([PathV 'temp\parameters.mat'], 'Pname', 'Fh', 'Fw', 'Ft', 'Vh', 'Vw', 'Vt', 'Lh', 'Lw', 'Lt', 'Lam');
	Sh = Fh / Vh * 0.5 / sqrt(2 * log(2));
	Sw = Fw / Vw * 0.5 / sqrt(2 * log(2));
	St = Ft / Vt * 0.5 / sqrt(2 * log(2));
	Lsize = (2 * Lh + 1) * (2 * Lw + 1) * (2 * Lt + 1);
	if ~isdir('temp')
		mkdir('temp');
	end
	save('temp\Ppsf.mat', 'Fh', 'Vh', 'Sh', 'Fw', 'Vw', 'Sw', 'Ft', 'Vt', 'St');
	save('temp\Ploc.mat', 'Lh', 'Lw', 'Lt', 'Lsize');
	save('temp\Plam.mat', 'Lam');
% ----------------------------------------------------------------
elseif strcmp(Apara, 'Image-based')
% ----------------------------------------------------------------
% Read the last set of parameters
	PathV = which('vanuc.m');
	PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
	while 1
		if isfile([PathV 'temp\parameters.mat'])
			load([PathV 'temp\parameters.mat']);
			if exist('Pname') && ischar(Pname) && exist('Vh') && isnumeric(Vh) && Vh>0 && Vh<100 && exist('Vw') && isnumeric(Vw) && Vw>0 && Vw<100 && exist('Vt') && isnumeric(Vt) && Vt>0 && Vt<100 && exist('Lh') && isnumeric(Lh) && Lh>=0 && Lh<15 && ~mod(Lh, 1) && exist('Lw') && isnumeric(Lw) && Lw>=0 && Lw<15 && ~mod(Lw, 1) && exist('Lt') && isnumeric(Lt) && Lt>=0 && Lt<15 && ~mod(Lt, 1) && exist('Lam') && isnumeric(Lam) && Lam>=0
				break;
			end
		end
		Pname = '(default)'
		Fh = 8;
		Fw = 8;
		Ft = 8;
		Vh = Vx;
		Vw = Vy;
		Vt = Vz;
		Lh = ceil(20 / Vh);
		Lw = ceil(20 / Vw);
		Lt = ceil(20 / Vt);
		Lam = 1;
		if ~isdir([PathV 'temp'])
			mkdir([PathV 'temp']);
		end
		save([PathV 'temp\parameters.mat'], 'Pname', 'Fh', 'Fw', 'Ft', 'Vh', 'Vw', 'Vt', 'Lh', 'Lw', 'Lt', 'Lam');
	end
	
% Set parameters
	while 1
		if round(Vh - Vx, 5)==0 && round(Vw - Vy, 5)==0 && round(Vt - Vz, 5)==0
			answer = questdlg({['Protocol name:                  ' Pname], ['Voxel size(x,y,z):                ' num2str(Vh) 'mm x ' num2str(Vw) 'mm x ' num2str(Vt) 'mm'], ['Local region size(x,y,z):       ' num2str(Lh * 2 + 1) 'vx x ' num2str(Lw * 2 + 1) 'vx x ' num2str(Lt * 2 + 1) 'vx'], ['Regularization parameter:    ' num2str(Lam)]}, 'Confirmation of PET parameters', 'OK', 'Other Protocol', 'Dilect Input', 'OK');
			if strcmp(answer, 'OK')
				break
			elseif strcmp(answer, 'Other Protocol')
				Plist = dir([PathV 'parameters\*.mat']);
				if numel(Plist)
					Pselect = listdlg('SelectionMode', 'single', 'ListString', {Plist.name});
					if numel(Pselect)==0
						return
					end
					load([PathV 'parameters\' Plist(Pselect).name]);
					clear Plist
					if exist('Pname') && ischar(Pname) && exist('Vh') && isnumeric(Vh) && Vh>0 && Vh<100 && exist('Vw') && isnumeric(Vw) && Vw>0 && Vw<100 && exist('Vt') && isnumeric(Vt) && Vt>0 && Vt<100 && exist('Lh') && isnumeric(Lh) && Lh>=0 && Lh<20 && ~mod(Lh, 1) && exist('Lw') && isnumeric(Lw) && Lw>=0 && Lw<20 && ~mod(Lw, 1) && exist('Lt') && isnumeric(Lt) && Lt>=0 && Lt<20 && ~mod(Lt, 1) && exist('Lam') && isnumeric(Lam) && Lam>=0
						if round(Vh - Vx, 5)==0 && round(Vw - Vy, 5)==0 && round(Vt - Vz, 5)==0
							continue
						end
						warning('The size does not match.');
					else
						warning('Some parameters are incorrect.');
					end
				else
					warning('No protocol registered.');
					clear Plist
				end
			elseif ~strcmp(answer, 'Dilect Input')
				return
			end
		end
		Input = inputdlg({'Distance_x from center of local region (vx)   (dx=0,1,2,3,4)', 'Distance_y from center of local region (vx)   (dy=0,1,2,3,4)', 'Distance_z from center of local region (vx)   (dz=0,1,2,3,4)', 'Regularization parameter'}, 'Input Parameters', [1 30], {num2str(Lh), num2str(Lw), num2str(Lt), num2str(Lam)});
		if numel(Input)==0
			return
		end
		Lh = str2num(Input{1});
		Lw = str2num(Input{2});
		Lt = str2num(Input{3});
		Lam = str2num(Input{4});
		clear Input
		Pname = '(no name)';
		Vh = Vx;
		Vw = Vy;
		Vt = Vz;
	end
	Lsize = (2 * Lh + 1) * (2 * Lw + 1) * (2 * Lt + 1);
	if ~isdir('temp')
		mkdir('temp');
	end
	save('temp\Ploc.mat', 'Lh', 'Lw', 'Lt', 'Lsize');
	save('temp\Plam.mat', 'Lam');
% ----------------------------------------------------------------
else
	return
end
[VxM, VyM, VzM] = vanuc_imagepre('T1.nii', 'T1 volume', pathH);
if VxM==0 || VyM==0 || VzM==0
	return
end
close(fim);
fim = vanuc_disp('T1.nii');
fim.Name = 'Cropped MRI image';
pause(0.1);
VolPET = spm_vol('trim_PET.nii');
G = spm_read_vols(VolPET);
mkdir('temp');
cd temp
save G.mat G
cd ..
clear VolPET G

% Create isovoxel volume image
% ----------------------------------------------------------------
if (exist('T1.nii', 'file') ~= 2)
	Filenames = dir('co*.nii');
	movefile(Filenames(1).name, 'T1.nii');
	clear Filenames
end
voxsiz = [0.8 0.8 0.8];
VolMRI = spm_vol('T1.nii');
VV(1:2) = VolMRI;
Imv = spm_imatrix(VV(1).mat);
VV(1).dim = ceil(VolMRI.dim .* Imv(7:9) ./ voxsiz);
Imv(7:9) = voxsiz;
VV(1).mat = spm_matrix(Imv);
spm_reslice(VV, struct('interp', 1, 'mask', 0, 'mean', 0, 'which', 1, 'prefix', 'co'));
clear voxsiz VolMRI VV Imv
close(fim);
fim = vanuc_disp('coT1.nii');
fim.Name = 'Isovoxel image';
pause(0.1);

% Resize
% ----------------------------------------------------------------
V = spm_vol('coT1.nii');
dimat = spm_imatrix(V.mat);
px = min(abs(dimat(7:9)));
copyfile trim_PET.nii trim_PETcopy.nii
V = spm_vol('trim_PET.nii');
Y = spm_read_vols(V);
dimat = spm_imatrix(V.mat);
Pratio = ceil(abs(dimat(7:9)/px));
dimat(7:9) = dimat(7:9)./Pratio;
V.mat = spm_matrix(dimat);
Size = size(Y);
Y2 = zeros(Size.* Pratio);
V.dim = size(Y2);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1)) = Y;
		end
	end
end
spm_write_vol(V, Y2);
clear V dimat px Y Pratio Size Y2 i j k
disp(' preparing images D O N E')

% Segmentation & normalisation
% ----------------------------------------------------------------
DirH = pwd;
PathSPM = which('spm.m');
PathSPM = PathSPM(1 : strfind(PathSPM, 'spm.m') - 1);
PathTPM = strcat(PathSPM, '\tpm\TPM.nii');
if ~exist(PathTPM)
	cd(PathSPM);
	[fileTPM, pathTPM] = uigetfile('*.nii', ['Select Tissue Probability Map for normalization']);
	cd(DirH);
end
PathV = which('vanuc.m');
PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
matlabbatch{1}.spm.spatial.preproc.channel.vols = {'coT1.nii,1'};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {strcat(PathTPM,',1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {strcat(PathTPM,',2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {strcat(PathTPM,',3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {strcat(PathTPM,',4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {strcat(PathTPM,',5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {strcat(PathTPM,',6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear matlabbatch
close(fim);
fim = vanuc_disp('c1coT1.nii');
fim.Name = 'Gray matter image';
pause(0.1);
matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {'coT1.nii,1'};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {strcat(PathV,'\standard\TPM.nii')};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;
matlabbatch{2}.spm.util.defs.comp{1}.def = {'y_coT1.nii'};
matlabbatch{2}.spm.util.defs.out{1}.push.fnames = {
                                                   strcat(PathV,'\standard\MNI_VOI_CBL.nii')
                                                   strcat(PathV,'\standard\MNI_VOI_STR.nii')
                                                   };
matlabbatch{2}.spm.util.defs.out{1}.push.weight = {''};
matlabbatch{2}.spm.util.defs.out{1}.push.savedir.saveusr = {''};
matlabbatch{2}.spm.util.defs.out{1}.push.fov.file = {'coT1.nii'};
matlabbatch{2}.spm.util.defs.out{1}.push.preserve = 0;
matlabbatch{2}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear matlabbatch PathV

% Segmentation of cerebellar GM
% ----------------------------------------------------------------
Volc1 = spm_vol('c1coT1.nii');
Yc1 = spm_read_vols(Volc1);
Volcbl = spm_vol('wMNI_VOI_CBL.nii');
Ycbl = spm_read_vols(Volcbl);
Yc1a = Yc1 .* Ycbl;
clear Volcbl Ycbl
Volc1a = Volc1;
Volc1a.fname = 'c1acoT1.nii';
spm_write_vol(Volc1a, Yc1a);
Yc1b = Yc1 - Yc1a;
clear Yc1 Yc1a Volc1a
Volc1b = Volc1;
Volc1b.fname = 'c1bcoT1.nii';
spm_write_vol(Volc1b, Yc1b);
clear Yc1b Volc1 Volc1b
disp(' CBL segmentation D O N E')
delete('wMNI_VOI_CBL.nii');
delete('wMNI_VOI_STR.nii');
close(fim);
fim = vanuc_disp('c1bcoT1.nii');
fim.Name = 'Gray matter without Cerebellum';
pause(0.1);

% Coregister
% ----------------------------------------------------------------
copyfile coT1.nii coT1copy.nii
copyfile c1acoT1.nii c1acoT1copy.nii
copyfile c1bcoT1.nii c1bcoT1copy.nii
copyfile c2coT1.nii c2coT1copy.nii
copyfile c3coT1.nii c3coT1copy.nii
copyfile c4coT1.nii c4coT1copy.nii
copyfile c5coT1.nii c5coT1copy.nii
copyfile c6coT1.nii c6coT1copy.nii

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'trim_PET.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'coT1.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                   'c1acoT1.nii,1'
                                                   'c1bcoT1.nii,1'
                                                   'c2coT1.nii,1'
                                                   'c3coT1.nii,1'
                                                   'c4coT1.nii,1'
                                                   'c5coT1.nii,1'
                                                   'c6coT1.nii,1'
                                                   };
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

delete coT1.nii
delete c1acoT1.nii
delete c1bcoT1.nii
delete c2coT1.nii
delete c3coT1.nii
delete c4coT1.nii
delete c5coT1.nii
delete c6coT1.nii
movefile coT1copy.nii coT1.nii
movefile c1acoT1copy.nii c1acoT1.nii
movefile c1bcoT1copy.nii c1bcoT1.nii
movefile c2coT1copy.nii c2coT1.nii
movefile c3coT1copy.nii c3coT1.nii
movefile c4coT1copy.nii c4coT1.nii
movefile c5coT1copy.nii c5coT1.nii
movefile c6coT1copy.nii c6coT1.nii

% Resize
% ----------------------------------------------------------------
delete trim_PET.nii
movefile trim_PETcopy.nii trim_PET.nii
V = spm_vol('trim_PET.nii');
Vr = spm_vol('rc1acoT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc1bcoT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc2coT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc3coT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc4coT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc5coT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
Vr = spm_vol('rc6coT1.nii');
Pratio = Vr.dim./V.dim;
Y2 = spm_read_vols(Vr);
Size = V.dim;
Y = zeros(Size);
for i=1:Pratio(1)
	for j=1:Pratio(2)
		for k=1:Pratio(3)
			Y = Y + Y2(i:Pratio(1):i+Pratio(1)*(Size(1)-1), j:Pratio(2):j+Pratio(2)*(Size(2)-1), k:Pratio(3):k+Pratio(3)*(Size(3)-1));
		end
	end
end
Y = Y / Pratio(1) / Pratio(2) / Pratio(3);
Vr.dim = V.dim;
Vr.mat = V.mat;
spm_write_vol(Vr, Y);
clear V Vr Pratio Y2 Size Y i j k
close(fim);
fim = vanuc_disp('rc1bcoT1.nii');
fim.Name = 'Coregistered GM image';
pause(0.1);

% Tissue probability map
% ----------------------------------------------------------------
VolC1 = spm_vol('rc1acoT1.nii');
M = zeros([VolC1.dim, 7]);
M(:,:,:,1) = spm_read_vols(VolC1);
VolC2 = spm_vol('rc1bcoT1.nii');
M(:,:,:,2) = spm_read_vols(VolC2);
VolC3 = spm_vol('rc2coT1.nii');
M(:,:,:,3) = spm_read_vols(VolC3);
VolC4 = spm_vol('rc3coT1.nii');
M(:,:,:,4) = spm_read_vols(VolC4);
VolC5 = spm_vol('rc4coT1.nii');
M(:,:,:,5) = spm_read_vols(VolC5);
VolC6 = spm_vol('rc5coT1.nii');
M(:,:,:,6) = spm_read_vols(VolC6);
VolC7 = spm_vol('rc6coT1.nii');
M(:,:,:,7) = spm_read_vols(VolC7);
Msum = sum(M, 4);
Mask = 1 - Msum;
Mask(find(Mask < 0)) = 0;
Msum = Msum + Mask;
M(:,:,:,7) = M(:,:,:,7) + Mask;
M(:,:,:,1) = M(:,:,:,1) ./ Msum;
M(:,:,:,2) = M(:,:,:,2) ./ Msum;
M(:,:,:,3) = M(:,:,:,3) ./ Msum;
M(:,:,:,4) = M(:,:,:,4) ./ Msum;
M(:,:,:,5) = M(:,:,:,5) ./ Msum;
M(:,:,:,6) = M(:,:,:,6) ./ Msum;
M(:,:,:,7) = M(:,:,:,7) ./ Msum;
clear Msum
spm_write_vol(VolC1, M(:,:,:,1));
spm_write_vol(VolC2, M(:,:,:,2));
spm_write_vol(VolC3, M(:,:,:,3));
spm_write_vol(VolC4, M(:,:,:,4));
spm_write_vol(VolC5, M(:,:,:,5));
spm_write_vol(VolC6, M(:,:,:,6));
spm_write_vol(VolC7, M(:,:,:,7));
mkdir('temp');
cd temp
save M.mat M
cd ..
clear M VolC1 VolC2 VolC3 VolC4 VolC5 VolC6 VolC7
delete rc1acoT1.nii
delete rc1bcoT1.nii
delete rc2coT1.nii
delete rc3coT1.nii
delete rc4coT1.nii
delete rc5coT1.nii
delete rc6coT1.nii
disp(' Moutput D O N E');

% Image-based parameters estimation
% ----------------------------------------------------------------
if strcmp(Apara, 'Image-based')
	if ~exist('xyLlimit', 'var');
		xyLlimit = 0;
		xyUlimit = 6;
		xyprecision = 0.01;
		zLlimit = 0;
		zUlimit = 6;
		zprecision = 0.01;
	end
	
	cd temp
	load G.mat
	load M.mat
	if ~(exist('Vh', 'var') && exist('Vw', 'var') && exist('Vt', 'var')) && exist('Ppsf.mat', 'file')
		load Ppsf.mat
	end
	if ~(exist('Vh', 'var') && exist('Vw', 'var') && exist('Vt', 'var'))
		savepara = 0;
	end
	[Sx, Sz] = vanuc_estimateresol(G, M, xyLlimit, xyUlimit, xyprecision, zLlimit, zUlimit, zprecision, 'narrow');
	Sy = Sx;
	
	if ~exist('savepara')
		savepara = 1;
	end
	if savepara
		Fh = Sx * Vh * 2 * sqrt(2 * log(2));
		Fw = Sy * Vw * 2 * sqrt(2 * log(2));
		Ft = Sz * Vt * 2 * sqrt(2 * log(2));
		Sh = Sx;
		Sw = Sy;
		St = Sz;
		if exist('manuinput') && manuinput
			Input = inputdlg({'FWHM_x of PSF (mm)   (0<fx<100)', 'FWHM_y of PSF(mm)   (0<fy<100)', 'FWHM_z of PSF(mm)   (0<fz<100)'}, 'Input Parameters', [1 30], {num2str(Fh), num2str(Fw), num2str(Ft)});
			if numel(Input)==0
				return
			end
			Fh = str2num(Input{1});
			Fw = str2num(Input{2});
			Ft = str2num(Input{3});
		end
		save('Ppsf.mat', 'Fh', 'Vh', 'Sh', 'Fw', 'Vw', 'Sw', 'Ft', 'Vt', 'St');
	end
	cd ..
end

close(fim);
Success = 1;

end