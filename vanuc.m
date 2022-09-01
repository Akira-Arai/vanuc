function vanuctest2
% VANUC: Voxel-wise Anatomical-region-based Non-Uniformity Correction
% 
% VANUC software is an image processing tool using VANUC
% method developed by Akira Arai (Kousei Sendai Clinic).
% VANUC is one of the partial volume effect correction method
% for PET or SPECT using MRI-derived anatomical information.
% It is based on probability theory assuming a non-uniformity
% and is robust for non-uniform distributions in tissues.
% 
% This software works with MATLAB and SPM12.
% Therefore, they are required to be installed prior to use.
% 
% This function launches the start menu of VANUC.
%     Start: This function starts partial volume effect correction.
%     dcm>nii: This function converts DICOM files to NIfTI files.
%     Orientation: This function re-orientates the directions.
%     Crop: This function crops excess area from the image and
%           leaves the head area.
%     Phantom: This function analyses phantom data and
%              determine parameters.
%     about VANUC
%     Mannual
%     Quit: This function closes the VANUC menu window.
% 
% Referrence:
%     Akira Arai, Yuriko Kagaya, Kentaro Takanami, Kei Takase.
%     A novel partial volume effect correction method which
%     allows for heterogeneous distribution: The potential
%     advantages in the white matter activity estimation on
%     FDG-PET. J Nucl Med. 2016 May 1;57(supplement 2):1924
% No papers on the algorithm of this VANUC method have
% been published to date (2022), because the developer has
% been prohibited from submitting all scientific papers to
% scientific journals since 2016.
% However, publication of results using this software is free for
% all users as long as it complies with the license agreement.
% 
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

PathV = which('vanuc.m');
PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
r = groot;
I=imread([PathV 'VANUClogoM.jpg']);
fmenu = figure;
imshow(I, 'Border', 'tight');
fmenu.Color = [0 0 0];
fmenu.Position = [20 r.MonitorPositions(4)-340 260 285];
fmenu.Name = 'Menu';
fmenu.NumberTitle = 'off';
fmenu.MenuBar='none';
fmenu.ToolBar='none';
text(11,260,'Partial Volume Effect Correction tool using SPM.','Color',[0.5 0.5 0.5],'FontSize',8);
text(167,142,'VANUC Ver. 2','Color',[0.9 0.9 0.9],'FontSize',9);
c1 = uicontrol;
c1.Position = [166 88 80 35];
c1.String = 'Start';
c1.FontSize = 19;
c1.Callback = @startVANUC;
c2a = uicontrol;
c2a.BackgroundColor = [0.7 0.7 0.7];
c2a.Position = [181 243 80 25];
c2a.String = 'Phantom';
c2a.FontSize = 11;
c2a.Callback = @startPhantom;
c2b = uicontrol;
c2b.BackgroundColor = [0.7 0.7 0.7];
c2b.Position = [1 243 60 25];
c2b.String = 'dcm>nii';
c2b.FontSize = 11;
c2b.Callback = @Convert;
c2c = uicontrol;
c2c.BackgroundColor = [0.7 0.7 0.7];
c2c.Position = [141 243 40 25];
c2c.String = 'Crop';
c2c.FontSize = 11;
c2c.Callback = @Cropping;
c2d = uicontrol;
c2d.BackgroundColor = [0.7 0.7 0.7];
c2d.Position = [61 243 80 25];
c2d.String = 'Orientation';
c2d.FontSize = 11;
c2d.Callback = @Reorientation;
c3 = uicontrol;
c3.BackgroundColor = [0.7 0.7 0.7];
c3.Position = [126 268 90 18];
c3.String = 'Mannual';
c3.FontSize = 10;
c3.Callback = @Help;
c4 = uicontrol;
c4.BackgroundColor = [0.7 0.7 0.7];
c4.Position = [216 268 45 18];
c4.String = 'Quit';
c4.FontSize = 10;
c4.Callback = @Quit;
c5 = uicontrol;
c5.BackgroundColor = [0.7 0.7 0.7];
c5.Position = [1 268 125 18];
c5.String = 'about VANUC';
c5.FontSize = 10;
c5.Callback = @Setting;
	function startVANUC(src, event)
		disp(datetime);
		Success = vanuc_pre;
		if ~Success
			return
		end
		disp(datetime);
		vanuc_pvec;
		disp(datetime);
		vanuc_post;
		disp(datetime);
	end
	function startPhantom(src, event)
		disp(datetime);
		[Sx, Sy, Sz] = vanuc_phantom;
		if Sx==0 || Sy==0 || Sz==0
			return
		end
		disp(datetime);
	end
	function Convert(src, event)
		disp(datetime);
		file = vanuc_dcmconvert;
		if numel(file)==0
			return
		end
		disp(datetime);
	end
	function Cropping(src, event)
		disp(datetime);
		[file, path] = uigetfile('*.nii', 'Image');
		if ~ischar(file)
			return;
		else
			fnamelength = numel(file);
			if fnamelength<4 || ~strcmp(file(fnamelength - 3 : fnamelength), '.nii')
				return;
			end
			if strcmp([pwd, '\'], path)
				V = spm_vol(file);
				if ~isequal(V.dt, [16 1])
					newfile = [file(1 : fnamelength - 4), '_n.nii'];
					copyfile([path, file], newfile);
					Y = spm_read_vols(V);
					V.dt = [16 1];
					V.fname = newfile;
					spm_write_vol(V, Y);
					file = newfile;
					clear Y
				end
			else
				copyfile([path, file], file);
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
		end
		vanuc_crop(file);
		disp(datetime);
	end
	function Reorientation(src, event)
		disp(datetime);
		[file, path] = uigetfile('*.nii', 'Image');
		if ~ischar(file)
			return;
		else
			fnamelength = numel(file);
			if fnamelength<4 || ~strcmp(file(fnamelength - 3 : fnamelength), '.nii')
				return;
			end
			if strcmp([pwd, '\'], path)
				V = spm_vol(file);
				if ~isequal(V.dt, [16 1])
					newfile = [file(1 : fnamelength - 4), '_n.nii'];
					copyfile([path, file], newfile);
					Y = spm_read_vols(V);
					V.dt = [16 1];
					V.fname = newfile;
					spm_write_vol(V, Y);
					file = newfile;
					clear Y
				end
			else
				copyfile([path, file], file);
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
		end
		vanuc_reorientation(file);
		disp(datetime);
	end
	function Quit(src, event)
		close(fmenu);
	end
	function Setting(src, event)
		PathV = which('vanuc.m');
		PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
		open([PathV 'vanuc_paper.pdf']);
	end
	function Help(src, event)
		PathV = which('vanuc.m');
		PathV = PathV(1 : strfind(PathV, 'vanuc.m') - 1);
		open([PathV 'usermanual.pdf']);
	end
end