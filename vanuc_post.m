function varargout = vanuc_post(varargin)
% Output of partial volume effect corrected images
% 
% This function outputs various kinds of images constructed based on 
% partial volume effect correction (PVEC) data calculated by 
% vanuc_pvec function. 
% 
% Outputs: 
% PT*.nii or NM*.nii: PET or SPECT raw image
% SUV_PT*.nii: SUV map (only for PET)
% trim_PET.nii: Cropped PET or SPECT
% MR*.nii: MRI raw image
% T1.nii: Cropped MRI
% coT1.nii: Isovoxel MRI
% coT1_seg8.mat: Segmentation data by SPM
% c1~6coT1.nii: Tissues' maps (GM/WM/CSF/skull/soft tissue/air)
% y_coT1.nii: warpfield from original MRI space to MNI space
% c1acoT1.nii: Tissue map of cerebellar GM
% c1bcoT1.nii: Tissue map of the other GM
% rcoT1.nii: MRI coresistered with PET or SPECT
% rtrim_PET.nii: PET or SPECT coresistered with MRI
% PVCvanuc.nii: PVEC image of 7 segments by VANUC
% PVCmg.nii: PVEC image of 7 segments by Muller-Gartner
% PVCvanuc.nii: PVEC image of 7 segments by RBV
% rPVC*.nii: PVEC image coregistered with MRI
% PVCC*.nii: Composite PVEC image of all 7 segments
% PVCC*.nii: Composite PVE-reduced image of all 7 segments
% MNIspace folder: Normalized images
% temp folder: Raw data of PVEC
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Output corrected images
% ----------------------------------------------------------------
% VANUC
VolPET = spm_vol('trim_PET.nii');
VolPET.descrip = 'PVE-corrected with VANUC';
VolPET.private.descrip = 'vanuc by A Arai';
cd temp
load Rvanuc.mat
load G.mat
Gmax = 10 * max(G(:));
cd ..
Rvanuc(find(Rvanuc<0)) = 0;
Rvanuc(find(Rvanuc>Gmax)) = Gmax;
VolPET.fname = 'PVCvanuc1.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,1));
VolPET.fname = 'PVCvanuc2.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,2));
VolPET.fname = 'PVCvanuc3.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,3));
VolPET.fname = 'PVCvanuc4.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,4));
VolPET.fname = 'PVCvanuc5.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,5));
VolPET.fname = 'PVCvanuc6.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,6));
VolPET.fname = 'PVCvanuc7.nii';
spm_write_vol(VolPET, Rvanuc(:,:,:,7));
clear Rvanuc
fim = vanuc_disp('PVCvanuc2.nii', 1);
% ===>
fim.Name = 'PVEC image by VANUC';
pause(0.1);

% mMG
VolPET.descrip = 'PVE-corrected with mMG';
cd temp
load Rmg.mat
cd ..
Rmg(find(Rmg<0)) = 0;
Rmg(find(Rmg>Gmax)) = Gmax;
VolPET.fname = 'PVCmg1.nii';
spm_write_vol(VolPET, Rmg(:,:,:,1));
VolPET.fname = 'PVCmg2.nii';
spm_write_vol(VolPET, Rmg(:,:,:,2));
VolPET.fname = 'PVCmg3.nii';
spm_write_vol(VolPET, Rmg(:,:,:,3));
VolPET.fname = 'PVCmg4.nii';
spm_write_vol(VolPET, Rmg(:,:,:,4));
VolPET.fname = 'PVCmg5.nii';
spm_write_vol(VolPET, Rmg(:,:,:,5));
VolPET.fname = 'PVCmg6.nii';
spm_write_vol(VolPET, Rmg(:,:,:,6));
VolPET.fname = 'PVCmg7.nii';
spm_write_vol(VolPET, Rmg(:,:,:,7));
clear Rmg

% RBV
VolPET.descrip = 'PVE-corrected with RBV';
cd temp
load Rrbv.mat
cd ..
Rrbv(find(Rrbv<0)) = 0;
Rrbv(find(Rrbv>Gmax)) = Gmax;
VolPET.fname = 'PVCrbv1.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,1));
VolPET.fname = 'PVCrbv2.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,2));
VolPET.fname = 'PVCrbv3.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,3));
VolPET.fname = 'PVCrbv4.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,4));
VolPET.fname = 'PVCrbv5.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,5));
VolPET.fname = 'PVCrbv6.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,6));
VolPET.fname = 'PVCrbv7.nii';
spm_write_vol(VolPET, Rrbv(:,:,:,7));
clear Rrbv VolPET
disp('PVCIoutput D O N E')

% Coregister with MRI (isovoxel)
% ----------------------------------------------------------------
copyfile trim_PET.nii trim_PETcopy.nii
copyfile PVCvanuc1.nii PVCvanuc1copy.nii
copyfile PVCvanuc2.nii PVCvanuc2copy.nii
copyfile PVCvanuc3.nii PVCvanuc3copy.nii
copyfile PVCvanuc4.nii PVCvanuc4copy.nii
copyfile PVCvanuc5.nii PVCvanuc5copy.nii
copyfile PVCvanuc6.nii PVCvanuc6copy.nii
copyfile PVCvanuc7.nii PVCvanuc7copy.nii
copyfile PVCmg1.nii PVCmg1copy.nii
copyfile PVCmg2.nii PVCmg2copy.nii
copyfile PVCmg3.nii PVCmg3copy.nii
copyfile PVCmg4.nii PVCmg4copy.nii
copyfile PVCmg5.nii PVCmg5copy.nii
copyfile PVCmg6.nii PVCmg6copy.nii
copyfile PVCmg7.nii PVCmg7copy.nii
copyfile PVCrbv1.nii PVCrbv1copy.nii
copyfile PVCrbv2.nii PVCrbv2copy.nii
copyfile PVCrbv3.nii PVCrbv3copy.nii
copyfile PVCrbv4.nii PVCrbv4copy.nii
copyfile PVCrbv5.nii PVCrbv5copy.nii
copyfile PVCrbv6.nii PVCrbv6copy.nii
copyfile PVCrbv7.nii PVCrbv7copy.nii

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'coT1.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'trim_PET.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {
                                                   'PVCvanuc1.nii,1'
                                                   'PVCvanuc2.nii,1'
                                                   'PVCvanuc3.nii,1'
                                                   'PVCvanuc4.nii,1'
                                                   'PVCvanuc5.nii,1'
                                                   'PVCvanuc6.nii,1'
                                                   'PVCvanuc7.nii,1'
                                                   'PVCmg1.nii,1'
                                                   'PVCmg2.nii,1'
                                                   'PVCmg3.nii,1'
                                                   'PVCmg4.nii,1'
                                                   'PVCmg5.nii,1'
                                                   'PVCmg6.nii,1'
                                                   'PVCmg7.nii,1'
                                                   'PVCrbv1.nii,1'
                                                   'PVCrbv2.nii,1'
                                                   'PVCrbv3.nii,1'
                                                   'PVCrbv4.nii,1'
                                                   'PVCrbv5.nii,1'
                                                   'PVCrbv6.nii,1'
                                                   'PVCrbv7.nii,1'
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
close(fim);
fim = vanuc_disp('rPVCvanuc2.nii', 1);
fim.Name = 'Coregistered PVEC image';
pause(0.1);
matlabbatch{1}.spm.util.defs.comp{1}.def = {'y_coT1.nii'};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {
                                                   'coT1.nii'
                                                   'rtrim_PET.nii'
                                                   'rPVCvanuc1.nii'
                                                   'rPVCvanuc2.nii'
                                                   'rPVCvanuc3.nii'
                                                   'rPVCmg1.nii'
                                                   'rPVCmg2.nii'
                                                   'rPVCmg3.nii'
                                                   'rPVCrbv1.nii'
                                                   'rPVCrbv2.nii'
                                                   'rPVCrbv3.nii'
                                                   };
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear matlabbatch
close(fim);
fim = vanuc_disp('wrPVCvanuc2.nii', 1);
fim.Name = 'Normalized to NMI space';
pause(0.1);

delete trim_PET.nii
delete PVCvanuc1.nii
delete PVCvanuc2.nii
delete PVCvanuc3.nii
delete PVCvanuc4.nii
delete PVCvanuc5.nii
delete PVCvanuc6.nii
delete PVCvanuc7.nii
delete PVCmg1.nii
delete PVCmg2.nii
delete PVCmg3.nii
delete PVCmg4.nii
delete PVCmg5.nii
delete PVCmg6.nii
delete PVCmg7.nii
delete PVCrbv1.nii
delete PVCrbv2.nii
delete PVCrbv3.nii
delete PVCrbv4.nii
delete PVCrbv5.nii
delete PVCrbv6.nii
delete PVCrbv7.nii
movefile trim_PETcopy.nii trim_PET.nii
movefile PVCvanuc1copy.nii PVCvanuc1.nii
movefile PVCvanuc2copy.nii PVCvanuc2.nii
movefile PVCvanuc3copy.nii PVCvanuc3.nii
movefile PVCvanuc4copy.nii PVCvanuc4.nii
movefile PVCvanuc5copy.nii PVCvanuc5.nii
movefile PVCvanuc6copy.nii PVCvanuc6.nii
movefile PVCvanuc7copy.nii PVCvanuc7.nii
movefile PVCmg1copy.nii PVCmg1.nii
movefile PVCmg2copy.nii PVCmg2.nii
movefile PVCmg3copy.nii PVCmg3.nii
movefile PVCmg4copy.nii PVCmg4.nii
movefile PVCmg5copy.nii PVCmg5.nii
movefile PVCmg6copy.nii PVCmg6.nii
movefile PVCmg7copy.nii PVCmg7.nii
movefile PVCrbv1copy.nii PVCrbv1.nii
movefile PVCrbv2copy.nii PVCrbv2.nii
movefile PVCrbv3copy.nii PVCrbv3.nii
movefile PVCrbv4copy.nii PVCrbv4.nii
movefile PVCrbv5copy.nii PVCrbv5.nii
movefile PVCrbv6copy.nii PVCrbv6.nii
movefile PVCrbv7copy.nii PVCrbv7.nii

% Count normalise
% ----------------------------------------------------------------
% Reference count
cd temp
load M.mat
Size = size(M);
VX = Size(1) * Size(2) * Size(3);
Voxel = sum(reshape(M, VX, Size(4)), 1);
clear M Size VX
save Voxel.mat Voxel
load Rgtm.mat
CBL = Rgtm(1);
GM = (Rgtm(1) * Voxel(1) + Rgtm(2) * Voxel(2)) / (Voxel(1) + Voxel(2));
WB = (Rgtm(1) * Voxel(1) + Rgtm(2) * Voxel(2) + Rgtm(3) * Voxel(3)) / (Voxel(1) + Voxel(2) + Voxel(3));
clear Voxel Rgtm
cd ..

% Non-corrected
Vol = spm_vol('wrtrim_PET.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrtrim_PET_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrtrim_PET_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrtrim_PET_wb.nii';
spm_write_vol(Vol, R/WB);

% VANUC
Vol = spm_vol('wrPVCvanuc1.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCvanuc1_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCvanuc1_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCvanuc1_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCvanuc2.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCvanuc2_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCvanuc2_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCvanuc2_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCvanuc3.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCvanuc3_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCvanuc3_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCvanuc3_wb.nii';
spm_write_vol(Vol, R/WB);

% mMG
Vol = spm_vol('wrPVCmg1.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCmg1_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCmg1_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCmg1_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCmg2.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCmg2_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCmg2_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCmg2_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCmg3.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCmg3_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCmg3_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCmg3_wb.nii';
spm_write_vol(Vol, R/WB);

% RBV
Vol = spm_vol('wrPVCrbv1.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCrbv1_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCrbv1_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCrbv1_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCrbv2.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCrbv2_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCrbv2_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCrbv2_wb.nii';
spm_write_vol(Vol, R/WB);
Vol = spm_vol('wrPVCrbv3.nii');
R = spm_read_vols(Vol);
Vol.fname = 'wrPVCrbv3_cbl.nii';
spm_write_vol(Vol, R/CBL);
Vol.fname = 'wrPVCrbv3_gm.nii';
spm_write_vol(Vol, R/GM);
Vol.fname = 'wrPVCrbv3_wb.nii';
spm_write_vol(Vol, R/WB);
clear Vol R CBL GM WB
disp('reference D O N E')

mkdir('MNIspace');
movefile wcoT1.nii MNIspace
movefile wrtrim_PET.nii MNIspace
movefile wrPVCvanuc1.nii MNIspace
movefile wrPVCvanuc2.nii MNIspace
movefile wrPVCvanuc3.nii MNIspace
movefile wrPVCmg1.nii MNIspace
movefile wrPVCmg2.nii MNIspace
movefile wrPVCmg3.nii MNIspace
movefile wrPVCrbv1.nii MNIspace
movefile wrPVCrbv2.nii MNIspace
movefile wrPVCrbv3.nii MNIspace
movefile wrtrim_PET_cbl.nii MNIspace
movefile wrtrim_PET_gm.nii MNIspace
movefile wrtrim_PET_wb.nii MNIspace
movefile wrPVCvanuc1_cbl.nii MNIspace
movefile wrPVCvanuc1_gm.nii MNIspace
movefile wrPVCvanuc1_wb.nii MNIspace
movefile wrPVCvanuc2_cbl.nii MNIspace
movefile wrPVCvanuc2_gm.nii MNIspace
movefile wrPVCvanuc2_wb.nii MNIspace
movefile wrPVCvanuc3_cbl.nii MNIspace
movefile wrPVCvanuc3_gm.nii MNIspace
movefile wrPVCvanuc3_wb.nii MNIspace
movefile wrPVCmg1_cbl.nii MNIspace
movefile wrPVCmg1_gm.nii MNIspace
movefile wrPVCmg1_wb.nii MNIspace
movefile wrPVCmg2_cbl.nii MNIspace
movefile wrPVCmg2_gm.nii MNIspace
movefile wrPVCmg2_wb.nii MNIspace
movefile wrPVCmg3_cbl.nii MNIspace
movefile wrPVCmg3_gm.nii MNIspace
movefile wrPVCmg3_wb.nii MNIspace
movefile wrPVCrbv1_cbl.nii MNIspace
movefile wrPVCrbv1_gm.nii MNIspace
movefile wrPVCrbv1_wb.nii MNIspace
movefile wrPVCrbv2_cbl.nii MNIspace
movefile wrPVCrbv2_gm.nii MNIspace
movefile wrPVCrbv2_wb.nii MNIspace
movefile wrPVCrbv3_cbl.nii MNIspace
movefile wrPVCrbv3_gm.nii MNIspace
movefile wrPVCrbv3_wb.nii MNIspace

% Output composite images
% ----------------------------------------------------------------
% VANUC
Vol = spm_vol('rPVCvanuc1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rtrim_PET.nii');
Vol.fname = 'PVCCvanuc.nii';
Vol.descrip = 'PVE-corrected with VANUC';
Vol.private.descrip = 'vanuc by A Arai';
spm_write_vol(Vol, Rcomp);
close(fim);
fim = vanuc_disp('PVCCvanuc.nii', 1);
fim.Name = 'PVE-corrected image by VANUC';
pause(0.1);

% mMG
Vol = spm_vol('rPVCmg1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('PVCCvanuc.nii');
Vol.fname = 'PVCCmg.nii';
Vol.descrip = 'PVE-corrected with mMG';
spm_write_vol(Vol, Rcomp);

% RBV
Vol = spm_vol('rPVCrbv1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('c6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('PVCCvanuc.nii');
Vol.fname = 'PVCCrbv.nii';
Vol.descrip = 'PVE-corrected with RBV';
spm_write_vol(Vol, Rcomp);
clear Vol Rpvc Rcomp
disp('PVCCoutput D O N E')

% Output PVE-reducted images
% ----------------------------------------------------------------
% Smoothing
FWHM = [2.5 2.5 2.5];
matlabbatch{1}.spm.spatial.smooth.data = {
                                          'c1acoT1.nii,1'
                                          'c1bcoT1.nii,1'
                                          'c2coT1.nii,1'
                                          'c3coT1.nii,1'
                                          'c4coT1.nii,1'
                                          'c5coT1.nii,1'
                                          'c6coT1.nii,1'
                                          };
matlabbatch{1}.spm.spatial.smooth.fwhm = FWHM;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'PET');
spm_jobman('run', matlabbatch);
clear FWHM matlabbatch

% VANUC
Vol = spm_vol('rPVCvanuc1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCvanuc7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('PVCCvanuc.nii');
Vol.fname = 'PVRCvanuc.nii';
Vol.descrip = 'PVE-reducted with VANUC';
spm_write_vol(Vol, Rcomp);
close(fim);
fim = vanuc_disp('PVRCvanuc.nii', 1);
fim.Name = 'PVE-reducted image by VANUC';
pause(0.1);

% mMG
Vol = spm_vol('rPVCmg1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCmg7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('PVCCvanuc.nii');
Vol.fname = 'PVRCmg.nii';
Vol.descrip = 'PVE-reducted with mMG';
spm_write_vol(Vol, Rcomp);

% RBVâÊëúèoóÕ
Vol = spm_vol('rPVCrbv1.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1acoT1.nii');
Rcomp = Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv2.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc1bcoT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv3.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc2coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv4.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc3coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv5.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc4coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv6.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc5coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('rPVCrbv7.nii');
Rpvc = spm_read_vols(Vol);
Vol = spm_vol('sc6coT1.nii');
Rcomp = Rcomp + Rpvc .* spm_read_vols(Vol);
Vol = spm_vol('PVCCvanuc.nii');
Vol.fname = 'PVRCrbv.nii';
Vol.descrip = 'PVE-reducted with RBV';
spm_write_vol(Vol, Rcomp);
clear Vol Rcomp Rpvc
delete 'sc1acoT1.nii'
delete 'sc1bcoT1.nii'
delete 'sc2coT1.nii'
delete 'sc3coT1.nii'
delete 'sc4coT1.nii'
delete 'sc5coT1.nii'
delete 'sc6coT1.nii'
disp('PVRCoutput D O N E')

spm_file_merge({'PVCvanuc1.nii';'PVCvanuc2.nii';'PVCvanuc3.nii';'PVCvanuc4.nii';'PVCvanuc5.nii';'PVCvanuc6.nii';'PVCvanuc7.nii'}, 'PVCvanuc.nii');
delete PVCvanuc1.nii
delete PVCvanuc2.nii
delete PVCvanuc3.nii
delete PVCvanuc4.nii
delete PVCvanuc5.nii
delete PVCvanuc6.nii
delete PVCvanuc7.nii
spm_file_merge({'PVCmg1.nii';'PVCmg2.nii';'PVCmg3.nii';'PVCmg4.nii';'PVCmg5.nii';'PVCmg6.nii';'PVCmg7.nii'}, 'PVCmg.nii');
delete PVCmg1.nii
delete PVCmg2.nii
delete PVCmg3.nii
delete PVCmg4.nii
delete PVCmg5.nii
delete PVCmg6.nii
delete PVCmg7.nii
spm_file_merge({'PVCrbv1.nii';'PVCrbv2.nii';'PVCrbv3.nii';'PVCrbv4.nii';'PVCrbv5.nii';'PVCrbv6.nii';'PVCrbv7.nii'}, 'PVCrbv.nii');
delete PVCrbv1.nii
delete PVCrbv2.nii
delete PVCrbv3.nii
delete PVCrbv4.nii
delete PVCrbv5.nii
delete PVCrbv6.nii
delete PVCrbv7.nii
spm_file_merge({'rPVCvanuc1.nii';'rPVCvanuc2.nii';'rPVCvanuc3.nii';'rPVCvanuc4.nii';'rPVCvanuc5.nii';'rPVCvanuc6.nii';'rPVCvanuc7.nii'}, 'rPVCvanuc.nii');
delete rPVCvanuc1.nii
delete rPVCvanuc2.nii
delete rPVCvanuc3.nii
delete rPVCvanuc4.nii
delete rPVCvanuc5.nii
delete rPVCvanuc6.nii
delete rPVCvanuc7.nii
spm_file_merge({'rPVCmg1.nii';'rPVCmg2.nii';'rPVCmg3.nii';'rPVCmg4.nii';'rPVCmg5.nii';'rPVCmg6.nii';'rPVCmg7.nii'}, 'rPVCmg.nii');
delete rPVCmg1.nii
delete rPVCmg2.nii
delete rPVCmg3.nii
delete rPVCmg4.nii
delete rPVCmg5.nii
delete rPVCmg6.nii
delete rPVCmg7.nii
spm_file_merge({'rPVCrbv1.nii';'rPVCrbv2.nii';'rPVCrbv3.nii';'rPVCrbv4.nii';'rPVCrbv5.nii';'rPVCrbv6.nii';'rPVCrbv7.nii'}, 'rPVCrbv.nii');
delete rPVCrbv1.nii
delete rPVCrbv2.nii
delete rPVCrbv3.nii
delete rPVCrbv4.nii
delete rPVCrbv5.nii
delete rPVCrbv6.nii
delete rPVCrbv7.nii

end