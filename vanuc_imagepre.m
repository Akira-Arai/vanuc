function [Vx, Vy, Vz, pathH] = vanuc_imagepre(outname, modality, pathH)
% Brain image processing for partial volume effect correction
% 
% Input:
% outname (char): File name of output image
% modality (char): Image modality name
% pathH (char): Directory displayed first when selecting files
% 
% Return:
% Vx (double): Width of pixel
% Vy (double): Height of pixel
% Vz (double): Thickness of slice
% pathH (char): Directory where the original image was saved
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Image file reading
% ----------------------------------------------------------------

if ~exist('modality')==1
	modality = 'AnyType';
end
if ~exist('pathH', 'var')
	pathH = pwd;
end
Dir = pwd;
Vx = 0;
Vy = 0;
Vz = 0;
skipcrop = 0;

while 1
	answer = questdlg(['File format of ' modality ' image'], 'Select', 'DICOM (*.dcm)', 'NIfTI (*.nii)', 'DICOM (*.dcm)');
	if strcmp(answer, 'DICOM (*.dcm)')
		[file, pathH] = vanuc_dcmconvert('single', pathH);
		if numel(file)==0
			return
		end
		file = file{1};
	elseif strcmp(answer, 'NIfTI (*.nii)')
		cd(pathH);
		[file, pathH] = uigetfile('*.nii', [modality ' Image']);
		cd(Dir);
		if ~ischar(file)
			continue
		else
			fnamelength = numel(file);
			if fnamelength<4 || ~strcmp(file(fnamelength - 3 : fnamelength), '.nii')
				continue
			end
		end
		if strcmp([pwd, '\'], pathH)
			V = spm_vol(file);
			if ~isequal(V.dt, [16 1])
				newfile = [file(1 : fnamelength - 4), '_n.nii'];
				copyfile([pathH, file], newfile);
				V = spm_vol(file);
				Y = spm_read_vols(V);
				V.dt = [16 1];
				V.fname = newfile;
				spm_write_vol(V, Y);
				file = newfile;
				clear Y
			end
		else
			copyfile([pathH, file], file);
			V = spm_vol(file);
			if ~isequal(V.dt, [16 1])
				newfile = [file(1 : fnamelength - 4), '_n.nii'];
				Y = spm_read_vols(V);
				V.dt = [16 1];
				V.fname = newfile;
				spm_write_vol(V, Y);
				file = newfile;
				clear Y
			end
		end
		clear V
		skipcrop = questdlg(['Do you perform Cropping?'], 'Select', 'OK (recommended)', 'Skip', 'OK (recommended)');
		skipcrop = strcmp(skipcrop, 'Skip');
	else
		return
	end
	break
end

[fileout, Vx, Vy, Vz] = vanuc_crop(file, skipcrop);
if exist(['ro' file])
	delete(['ro' file]);
end
if exist(['fro' file])
	delete(['fro' file]);
end
movefile(fileout, outname);
end