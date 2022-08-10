function [file, pathH] = vanuc_dcmconvert(mul, pathH)
% Convert DICOM files to NIfTI files (Convert to SUV units)
% 
% Input:
% mul - (default): Multiple series
%       'single': Single series
% pathH (char): Directory displayed first when selecting files
% 
% Return:
% file (char): Output file name
% pathH (char): Directory where the original image was saved
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Select files
% ----------------------------------------------------------------
file = {};
Dir = pwd;
if ~exist('pathH', 'var')
	pathH = Dir;
end
while 1
	cd(pathH);
	pathH = uigetdir;
	if ~ischar(pathH)
		return
	end
	cd(pathH);
	DCMfiles = dir('**/*.dcm');
	cd(Dir);
	if numel(DCMfiles)==0
		disp('No DICOM file');
	else
		break;
	end
end
for k = 1 : numel(DCMfiles)
	files{k} = [DCMfiles(k).folder '\' DCMfiles(k).name];
end
hdr = spm_dicom_headers(files);

for k=1:numel(hdr)
	SeriesNum(k) = hdr{k}.SeriesNumber;
end
[NumList, FFile] = unique(SeriesNum, 'stable');
for k = 1:numel(NumList)
	NameList{k} = [num2str(NumList(k)),': ',hdr{FFile(k)}.SeriesDescription];
end

if exist('mul') && ischar(mul) && strcmp(mul, 'single')
	ListSelect = listdlg('SelectionMode', 'single', 'ListString', NameList);
else
	ListSelect = listdlg('ListString', NameList);
end

if numel(ListSelect)==0
	return
end

% Convert DICOM to NIfTI
% ----------------------------------------------------------------
for k = 1:numel(ListSelect)
	List{k} = find(SeriesNum == NumList(ListSelect(k)));
end

for l = 1:numel(List)
	if isfield(hdr{List{l}(1)}, 'NumberOfFrames')
		Sorg = hdr{List{l}(1)}.NumberOfFrames;
	else
		Sorg = 1;
	end
	
	if Sorg > 1
		if isfield(hdr{List{l}(1)}, 'TransferSyntaxUID') && strcmp(hdr{List{l}(1)}.TransferSyntaxUID, '1.2.840.10008.1.2.2');
			fileO = fopen(hdr{List{l}(1)}.Filename, 'r', 'b');
		else
			fileO = fopen(hdr{List{l}(1)}.Filename, 'r', 'l');
		end
		Worg = hdr{List{l}(1)}.Columns;
		Horg = hdr{List{l}(1)}.Rows;
		fseek(fileO, hdr{List{l}(1)}.StartOfPixelData, 'bof');
		Y = reshape(fread(fileO, Worg * Horg * Sorg, strcat('ubit', num2str(hdr{List{l}(1)}.BitsAllocated), '=>unit32')), [Worg, Horg, Sorg]);
		fclose(fileO);
		if strcmp(hdr{List{l}(1)}.PatientOrientation(2),'\') && strcmp(hdr{List{l}(1)}.PatientOrientation(4),' ')
			if strcmp(hdr{List{l}(1)}.PatientOrientation(1),'L')
				Y = Y(Worg:-1:1, :, Sorg:-1:1);
			end
			if strcmp(hdr{List{l}(1)}.PatientOrientation(3),'P')
				Y = Y(:, Horg:-1:1, Sorg:-1:1);
			end
		end
		V.fname = strcat(hdr{List{l}(1)}.Modality, num2str(hdr{List{l}(1)}.SeriesNumber), '.nii');
		V.dim = [Worg, Horg, Sorg];
		V.dt = [4 1];
		V.pinfo = [1; 0; 352];
		Vox = [hdr{List{l}(1)}.PixelSpacing', hdr{List{l}(1)}.SpacingBetweenSlices];
		V.mat = spm_matrix([-0.5*V.dim.*Vox 0 0 0 Vox 0 0 0]);
		V.n = [1 1];
		spm_write_vol(V, Y);
		file{l} = V.fname;
		clear V
	else
		file{l} = [hdr{List{l}(1)}.Modality num2str(hdr{List{l}(1)}.SeriesNumber) '.nii'];
		for m = 1 : numel(List{l})
			if isfield(hdr{List{l}(m)}, 'MRScaleSlope')
				hdr{List{l}(m)}.MRScaleSlope = hdr{List{l}(m)}.MRScaleSlope(1);
			end
			if isfield(hdr{List{l}(m)}, 'MRScaleIntercept')
				hdr{List{l}(m)}.MRScaleIntercept = hdr{List{l}(m)}.MRScaleIntercept(1);
			end
		end
		NII = spm_dicom_convert(hdr(List{l}));
		Find = strfind(NII.files{1}, '\');
		Nstart = 1;
		if numel(Find)>0
			Nstart = max(Find)+1;
		end
		fileNII = NII.files{1}(Nstart:numel(NII.files{1}));
		movefile(fileNII, file{l});
	end
	
	if isfield(hdr{List{l}(1)}, 'RadiopharmaceuticalInformationSequence') && isfield(hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}, 'RadiopharmaceuticalStartTime') && isfield(hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}, 'RadionuclideHalfLife') && isfield(hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}, 'RadionuclideTotalDose') && isfield(hdr{List{l}(1)}, 'PatientWeight')
		V = spm_vol(file{l});
		Y = spm_read_vols(V);
		TIMEhl = (hdr{List{l}(1)}.SeriesTime - hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}.RadiopharmaceuticalStartTime) / hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}.RadionuclideHalfLife;
		Ksuv = 1000 * hdr{List{l}(1)}.PatientWeight / hdr{List{l}(1)}.RadiopharmaceuticalInformationSequence{1}.RadionuclideTotalDose * 2 ^ TIMEhl;
		Y = Y * Ksuv;
		V.fname = ['SUV_' file{l}];
		V.dt = [16 1];
		V.descrip = 'SUVbw';
		V.private.descrip = 'SUV map by A Arai';
		spm_write_vol(V, Y);
		file{l} = V.fname;
	end	
end
end