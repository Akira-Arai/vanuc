function [Sxy, Sz] = vanuc_measureresol(G, R, xyLlimit, xyUlimit, xyprecision, zLlimit, zUlimit, zprecision)
% Image-based estimation of sigma of PSF
% 
% (x-sigma equals y-sigma)
% 
% Input:
% G (3D double): Observed image (PET or SPECT)
% R (3D double): True distribution
% xyLlimit (double): Lower limit of xy-sigma
% xyUlimit (double): Upper limit of xy-sigma
% xyprecision (double): Precision of xy-sigma
% zLlimit (double): Lower limit of z-sigma
% zUlimit (double): Upper limit of z-sigma
% zprecision (double): Precision of z-sigma
% 
% Return:
% Sxy (double): xy-sigma
% Sz (double): z-sigma
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Initial value
% ----------------------------------------------------------------
dxy = xyUlimit - xyLlimit;
Sxy = xyLlimit + dxy * (sqrt(5) - 1) * 0.5;
SSxy = Sxy;
dz = zUlimit - zLlimit;
Sz = zLlimit + dz * (sqrt(5) - 1) * 0.5;
SSz = Sz;
MSEmin = vanuc_trueMSE(G, R, Sxy, Sxy, Sz, 'narrow');
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
	MSE = vanuc_trueMSE(G, R, XY, XY, Sz, 'narrow');
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
	MSE = vanuc_trueMSE(G, R, Sxy, Sxy, Z, 'narrow');
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
	MSE = vanuc_trueMSE(G, R, XY, XY, Sz, 'narrow');
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
	MSE =  vanuc_trueMSE(G, R, Sxy, Sxy, Z, 'narrow');
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
		MSE =  vanuc_trueMSE(G, R, Sxy, Sxy, Z, 'narrow');
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
		MSE = vanuc_trueMSE(G, R, Sxy, Sxy, Z, 'narrow');
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
	MSE =  vanuc_trueMSE(G, R, XY, XY, Sz, 'narrow');
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
		MSE = vanuc_trueMSE(G, R, XY, XY, Sxy, 'narrow');
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
		MSE = vanuc_trueMSE(G, R, XY, XY, Sz, 'narrow');
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
disp('datetime');
close(fplot);

end