function [Sxy, Sz] = vanuc_estimateresol(G, M, xyLlimit, xyUlimit, xyprecision, zLlimit, zUlimit, zprecision, narrow, Nseg)
% Image-based estimation of sigma of PSF
% 
% (x-sigma equals y-sigma)
% 
% Input:
% G (3D double): Observed image (PET or SPECT)
% M (4D double): Tissue maps (multiple segments)
% xyLlimit (double): Lower limit of xy-sigma
% xyUlimit (double): Upper limit of xy-sigma
% xyprecision (double): Precision of xy-sigma
% zLlimit (double): Lower limit of z-sigma
% zUlimit (double): Upper limit of z-sigma
% zprecision (double): Precision of z-sigma
% narrow - 'narrow': Trimming before analysis
% Nseg (positive integer): Tissues with numbers less than this value
%                          will be left after trimming
%                          (default: segmentation number - 1)
% 
% Return:
% Sxy (double): xy-sigma
% Sz (double): z-sigma
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Trimming
% ----------------------------------------------------------------
disp(datetime)
sizeM = size(M);
Kseg = sizeM(4);
if ~exist('Nseg');
	Nseg = Kseg - 1;
end
R = sum(M(:, :, :, 1 : Nseg), 4);
if exist('narrow') && ischar(narrow) && strcmp(narrow, 'narrow')
	Mask = R >= 1;
	SUM = sum(sum(Mask, 2), 3);
	for Xmin = 1 : size(Mask, 1)
		if SUM(Xmin) > 0
			break
		end
	end
	for Xmax = size(Mask, 1) : -1 : Xmin
		if SUM(Xmax) > 0
			break
		end
	end
	SUM = sum(sum(Mask, 1), 3);
	for Ymin = 1 : size(Mask, 2)
		if SUM(Ymin) > 0
			break
		end
	end
	for Ymax = size(Mask, 2) : -1 : Ymin
		if SUM(Ymax) > 0
			break
		end
	end
	SUM = sum(sum(Mask, 1), 2);
	for Zmin = 1 : size(Mask, 3)
		if SUM(Zmin) > 0
			break
		end
	end
	for Zmax = size(Mask, 3) : -1 : Zmin
		if SUM(Zmax) > 0
			break
		end
	end
	G = G(Xmin : Xmax, Ymin : Ymax, Zmin : Zmax);
	M = M(Xmin : Xmax, Ymin : Ymax, Zmin : Zmax, :);
end
clear R

% Initial value
% ----------------------------------------------------------------
dxy = xyUlimit - xyLlimit;
Sxy = xyLlimit + dxy * (sqrt(5) - 1) * 0.5;
SSxy = Sxy;
dz = zUlimit - zLlimit;
Sz = zLlimit + dz * (sqrt(5) - 1) * 0.5;
SSz = Sz;
MSEmin = vanuc_GTMMSE(G, M, Sxy, Sxy, Sz);
fplot = figure;
fplot.Name = 'Plot';
fplot.NumberTitle = 'off';
plot(SSxy, SSz, '-o');
xlim([xyLlimit xyUlimit]);
ylim([zLlimit zUlimit]);
pause(0.001);

% Initial optimization
% ----------------------------------------------------------------
sign = -1;
c = dxy * (3 - sqrt(5)) * 0.5;
while dxy >= xyprecision
	dxy = dxy * (sqrt(5) - 1) * 0.5;
	c = c * (sqrt(5) - 1) * 0.5;
	XY = Sxy + sign * c;
	MSE = vanuc_GTMMSE(G, M, XY, XY, Sz);
	SSxy = [SSxy ; XY];
	SSz = [SSz ; Sz];
	if MSE < MSEmin
		Sxy = XY;
		MSEmin = MSE;
	else
		sign = - sign;
	end
	figure(fplot);
	plot(SSxy, SSz, '-o');
	xlim([xyLlimit xyUlimit]);
	ylim([zLlimit zUlimit]);
	pause(0.001);
end

sign = -1;
c = dz * (3 - sqrt(5)) * 0.5;
zMove = Sz;
while dz >= zprecision
	dz = dz * (sqrt(5) - 1) * 0.5;
	c = c * (sqrt(5) - 1) * 0.5;
	Z = Sz + sign * c;
	MSE = vanuc_GTMMSE(G, M, Sxy, Sxy, Z);
	SSxy = [SSxy ; Sxy];
	SSz = [SSz ; Z];
	if MSE < MSEmin
		Sz = Z;
		MSEmin = MSE;
	else
		sign = - sign;
	end
	figure(fplot);
	plot(SSxy, SSz, '-o');
	xlim([xyLlimit xyUlimit]);
	ylim([zLlimit zUlimit]);
	pause(0.001);
end
zMove = Sz - zMove;

% Iterative optimization
% ----------------------------------------------------------------
xyMlim = (xyUlimit - xyLlimit) * 0.5;
if Sxy < xyMlim
	sign = 1;
	c = (xyUlimit - Sxy) * (sqrt(5) - 1) * 0.5;
else
	sign = -1;
	c = (Sxy - xyLlimit) * (sqrt(5) - 1) * 0.5;
end
dxy = c * (sqrt(5) + 3) * 0.5;
xyMove = Sxy;
while dxy >= xyprecision
	dxy = dxy * (sqrt(5) - 1) * 0.5;
	c = c * (sqrt(5) - 1) * 0.5;
	XY = Sxy + sign * c;
	MSE = vanuc_GTMMSE(G, M, XY, XY, Sz);
	SSxy = [SSxy ; XY];
	SSz = [SSz ; Sz];
	if MSE < MSEmin && XY > xyLlimit && XY < xyUlimit
		Sxy = XY;
		MSEmin = MSE;
	else
		sign = - sign;
	end
	figure(fplot);
	plot(SSxy, SSz, '-o');
	xlim([xyLlimit xyUlimit]);
	ylim([zLlimit zUlimit]);
	pause(0.001);
end
xyMove = Sxy - xyMove;

while xyMove ~= 0 && zMove ~= 0
	c = zMove;
	zMove = Sz;
	Z = Sz + c;
	if Z > zUlimit || Z < zLlimit
		c = - c * (3 - sqrt(5)) * 0.5;
		Z = Sz + c;
	end
	MSE =  vanuc_GTMMSE(G, M, Sxy, Sxy, Z);
	SSxy = [SSxy ; Sxy];
	SSz = [SSz ; Z];
	figure(fplot);
	plot(SSxy, SSz, '-o');
	xlim([xyLlimit xyUlimit]);
	ylim([zLlimit zUlimit]);
	pause(0.001);
	if MSE < MSEmin
		Sz = Z;
		MSEmin = MSE;
	else
		c = - c * (3 - sqrt(5)) * 0.5;
	end
	while 1
		c = c * (sqrt(5) + 1) * 0.5;
		Z = Sz + c;
		if Z > zUlimit || Z < zLlimit
			dz = abs(c) * (sqrt(5) + 1) * 0.5;
			c = c * (sqrt(5) - 1) * 0.5;
			break;
		end
		MSE =  vanuc_GTMMSE(G, M, Sxy, Sxy, Z);
		SSxy = [SSxy ; Sxy];
		SSz = [SSz ; Z];
		figure(fplot);
		plot(SSxy, SSz, '-o');
		xlim([xyLlimit xyUlimit]);
		ylim([zLlimit zUlimit]);
		pause(0.001);
		if MSE < MSEmin
			Sz = Z;
			MSEmin = MSE;
		else
			dz = abs(c) * (sqrt(5) + 1) * 0.5;
			c = c * (sqrt(5) - 1) * 0.5;
			break;
		end
	end
	while dz >= zprecision
		dz = dz * (sqrt(5) - 1) * 0.5;
		c = c * (sqrt(5) - 1) * 0.5;
		Z = Sz + c;
		MSE = vanuc_GTMMSE(G, M, Sxy, Sxy, Z);
		SSxy = [SSxy ; Sxy];
		SSz = [SSz ; Z];
		if MSE < MSEmin && Z > zLlimit && Z < zUlimit
			Sz = Z;
			MSEmin = MSE;
		else
			c = - c;
		end
		figure(fplot);
		plot(SSxy, SSz, '-o');
		xlim([xyLlimit xyUlimit]);
		ylim([zLlimit zUlimit]);
		pause(0.001);
	end
	zMove = Sz - zMove;

	c = xyMove;
	xyMove = Sxy;
	XY = Sxy + c;
	if XY > xyUlimit || XY < xyLlimit
		c = - c * (3 - sqrt(5)) * 0.5;
		XY = Sxy + c;
	end
	MSE =  vanuc_GTMMSE(G, M, XY, XY, Sz);
	SSxy = [SSxy ; XY];
	SSz = [SSz ; Sz];
	figure(fplot);
	plot(SSxy, SSz, '-o');
	xlim([xyLlimit xyUlimit]);
	ylim([zLlimit zUlimit]);
	pause(0.001);
	if MSE < MSEmin
		Sxy = XY;
		MSEmin = MSE;
	else
		c = - c * (3 - sqrt(5)) * 0.5;
	end
	while 1
		c = c * (sqrt(5) + 1) * 0.5;
		XY = Sxy + c;
		if XY > xyUlimit || XY < xyLlimit
			dxy = abs(c) * (sqrt(5) + 1) * 0.5;
			c = c * (sqrt(5) - 1) * 0.5;
			break;
		end
		MSE = vanuc_GTMMSE(G, M, XY, XY, Sxy);
		SSxy = [SSxy ; XY];
		SSz = [SSz ; Sz];
		figure(fplot);
		plot(SSxy, SSz, '-o');
		xlim([xyLlimit xyUlimit]);
		ylim([zLlimit zUlimit]);
		pause(0.001);
		if MSE < MSEmin
			Sxy = XY;
			MSEmin = MSE;
		else
			dxy = abs(c) * (sqrt(5) + 1) * 0.5;
			c = c * (sqrt(5) - 1) * 0.5;
			break;
		end
	end
	while dxy >= xyprecision
		dxy = dxy * (sqrt(5) - 1) * 0.5;
		c = c * (sqrt(5) - 1) * 0.5;
		XY = Sxy + c;
		MSE = vanuc_GTMMSE(G, M, XY, XY, Sz);
		SSxy = [SSxy ; XY];
		SSz = [SSz ; Sz];
		if MSE < MSEmin && XY > xyLlimit && XY < xyUlimit
			Sxy = XY;
			MSEmin = MSE;
		else
			c = - c;
		end
		figure(fplot);
		plot(SSxy, SSz, '-o');
		xlim([xyLlimit xyUlimit]);
		ylim([zLlimit zUlimit]);
		pause(0.001);
	end
	xyMove = Sxy - xyMove;
end

disp('Sxy      Sz       MSE');
disp([num2str(Sxy, '%.6f') ' ' num2str(Sz, '%.6f') ' ' num2str(MSEmin, '%.6f')]);
disp(datetime);
close(fplot);

end