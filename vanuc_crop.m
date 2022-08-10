function [fileout, vx, vy, vz] = vanuc_crop(filein, skip)
% Reorientation and cropping the brain image
% 
% Input:
% filein (char): Input file name
% skip - 0: Reorientation and cropping
%        1: Skip cropping
% 
% Return:
% fileout (char): Output file name
% Vx (double): Width of pixel
% Vy (double): Height of pixel
% Vz (double): Thickness of slice
% 
% Output data:
% 'crop_*.nii': Cropped image
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Reorientation
% ----------------------------------------------------------------
rofilein = vanuc_reorientation(filein);
disp('start cropping');
disp(datetime);
V = spm_vol(rofilein);
Vec = spm_imatrix(V.mat);
vx = Vec(7);
vy = Vec(8);
vz = Vec(9);
if exist('skip') && skip
	fileout = ['crop_' rofilein];
	copyfile(rofilein, fileout);
	return
end

% Binarization
% ----------------------------------------------------------------
Y = spm_read_vols(V);
Scenter = ceil(size(Y, 3) / 2);
nx = ceil(10 / vx);
ny = ceil(10 / vy);
Per = median([Y(1 : nx, 1 : ny, Scenter) Y(size(Y,1) - nx + 1 : size(Y, 1), 1 : ny, Scenter) Y(1 : nx, size(Y, 2) - ny + 1 : size(Y, 2), Scenter) Y(size(Y,1) - nx + 1 : size(Y, 1), size(Y, 2) - ny + 1 : size(Y, 2), Scenter)], 1:2);
Cen = median(Y(ceil(size(Y, 1) / 2) - 3*nx + 1 : ceil(size(Y, 1) / 2) + 3*nx, ceil(size(Y, 2) / 2) - 3*ny + 1 : ceil(size(Y, 2) / 2) + 3*ny, Scenter), 1:2);
Threshold = (Per + Cen) / 2;
if Per < Cen
	Ybin = Y > Threshold;
elseif Per > Cen
	Ybin = Y < Threshold;
else
	Ybin = Y > 0;
end

% Determine center of the head
% ----------------------------------------------------------------
Zbin = sum(Ybin, 1 : 2);
z0 = Scenter;
S0 = Zbin(z0);
z1 = z0;
S1 = S0;
z2 = z1;
S2 = S1;
CP = 0;
for z = Scenter + 1 : numel(Zbin)
	if Zbin(z) >= S0
		z0 = z;
		S0 = Zbin(z);
		z1 = z0;
		S1 = S0;
		z2 = z1;
		S2 = S1;
		CP = 0;
		continue
	end
	if CP == 0 && Zbin(z) >= 0.75 * S0
		z1 = z;
		S1 = Zbin(z);
		z2 = z1;
		S2 = S1;
		continue
	end
	if CP ~= 2 && Zbin(z) >= 0.25 * S0
		z2 = z;
		S2 = Zbin(z);
		CP = 1;
		continue
	end
	CP = 2;
	if Zbin(z) == 0
		break
	end
end
z3 = z;
S3 = Zbin(z);
if CP ~= 2
	return
elseif z1 == z2
	return
end
d1 = (S1 - S0) / (z1 - z0);
d2 = (S2 - S1) / (z2 - z1);
dd = (d2 - d1) / (z2 - z0);
if dd >= 0
	return
end
zaxis = (z0 + z1  - d1 / dd) / 2;
Smax = S0 - dd * (z0 - zaxis) ^ 2;
zp = zaxis + sqrt(-1 * Smax / dd);
dT = 100 * (1 - sqrt(Smax * vx * vy / 18000)) / vz;
if dT > 0
	zp = zp + dT;
end
Zcenter = round(1.2 * zp - 110 / vz);
TOP = min(round(Zcenter + 105 / vz), numel(Zbin));
BASE = max(round(Zcenter - 115 / vz), 1);
clear Zbin
Xcenter = round(sum(sum(Ybin(:, :, Zcenter + round(15 / vz)), 2)' .* [1 : size(Ybin, 1)]) / sum(Ybin(:, :, Zcenter + round(15 / vz)), 1 : 2));
xstart = max(Xcenter - round(95 / vx), 1);
xend = min(Xcenter + round(95 / vx), size(Ybin, 1));
Ycenter = round(sum(sum(Ybin(:, :, Zcenter + round(15 / vz)), 1) .* [1 : size(Ybin, 2)]) / sum(Ybin(:, :, Zcenter + round(15 / vz)), 1 : 2));
ystart = max(Ycenter - round(110 / vy), 1);
yend = min(Ycenter + round(110 / vy), size(Ybin, 2));
Ycenter = min(Ycenter + round(17 / vy), size(Ybin, 2));

% Cropping and Saving
% ----------------------------------------------------------------
Y = Y(xstart : xend, ystart : yend, BASE : TOP);
V.dim = size(Y);
invM = inv(V.mat);
invM(1:3,4) = [Xcenter - xstart + 1, Ycenter - ystart + 1, Zcenter - BASE + 1];
V.mat = inv(invM);
V.fname = ['crop_' V.fname];
spm_write_vol(V, Y);
fileout = V.fname;

end