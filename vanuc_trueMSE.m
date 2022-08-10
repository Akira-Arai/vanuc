function MSE = vanuc_trueMSE(G, R, Sh, Sw, St, narrow)
% Estimate difference of convolution model from the observed image
% 
% Input:
% G (3D double): Observed image (PET or SPECT)
% R (3D double): True distribution
% Sh (double): x-sigma
% Sw (double): y-sigma
% St (double): z-sigma
% narrow - 'narrow': Trimming before analysis
% 
% Return:
% MSE (double): MSE of convolution model from the observed image
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

sizeM = size(G);
Worg = sizeM(2);
Horg = sizeM(1);
Sorg = sizeM(3);

% trimming
% ----------------------------------------------------------------
if exist('narrow') && ischar(narrow) && strcmp(narrow, 'narrow')
	M = R > 0;
	SUM = sum(sum(M, 2), 3);
	for Xmin = 1 : size(M, 1)
		if SUM(Xmin) > 0
			break
		end
	end
	for Xmax = Xmin : size(M, 1) - 1
		if SUM(Xmax + 1) == 0
			break
		end
	end
	SUM = sum(sum(M, 1), 3);
	for Ymin = 1 : size(M, 2)
		if SUM(Ymin) > 0
			break
		end
	end
	for Ymax = Ymin : size(M, 2) - 1
		if SUM(Ymax + 1) == 0
			break
		end
	end
	SUM = sum(sum(M, 1), 2);
	for Zmin = 1 : size(M, 3)
		if SUM(Zmin) > 0
			break
		end
	end
	for Zmax = Zmin : size(M, 3) - 1
		if SUM(Zmax + 1) == 0
			break
		end
	end
end

% Calculate point spread function (PSF)
% ----------------------------------------------------------------
PSFw = [1:Worg];
PSFh = [1:Horg];
PSFt = [1:Sorg];
c = exp(-1 / 2 / Sw^2);
PSFw = c.^((PSFw-1).^2);
c = exp(-1 / 2 / Sh^2);
PSFh = c.^((PSFh-1).^2);
PSFh = reshape(PSFh, [Horg,1]);
c = exp(-1 / 2 / St^2);
PSFt = c.^((PSFt-1).^2);
PSFt = reshape(PSFt, [1,1,Sorg]);
clear c
PSFw = PSFw / (sum(PSFw) * 2 - 1);
PSFh = PSFh / (sum(PSFh) * 2 - 1);
PSFt = PSFt / (sum(PSFt) * 2 - 1);
Wpsf = sum(sign(PSFw));
PSFw = PSFw(1:Wpsf);
Hpsf = sum(sign(PSFh));
PSFh = PSFh(1:Hpsf);
Tpsf = sum(sign(PSFt));
PSFt = PSFt(:,:,1:Tpsf);

% Convolution
% ----------------------------------------------------------------
Vk = zeros(Horg, Worg+Wpsf*2-2, Sorg);
Vk(:, Wpsf:Worg+Wpsf-1, :) = R * PSFw(1);
for n=2:Wpsf
	Vk(:, Wpsf-n+1:Worg+Wpsf-n, :) = Vk(:, Wpsf-n+1:Worg+Wpsf-n, :) + R * PSFw(n);
	Vk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) = Vk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) + R * PSFw(n);
end
R = Vk(:, Wpsf:Worg+Wpsf-1, :);
R(:, 1:Wpsf-1, :) = R(:, 1:Wpsf-1, :) + Vk(:, Wpsf-1:-1:1, :);
R(:, Worg-Wpsf+2:Worg, :) = R(:, Worg-Wpsf+2:Worg, :) + Vk(:, Worg+Wpsf*2-2:-1:Worg+Wpsf, :);
Vk = zeros(Horg+Hpsf*2-2, Worg, Sorg);
Vk(Hpsf:Horg+Hpsf-1, :, :) = R * PSFh(1);
for n=2:Hpsf
	Vk(Hpsf-n+1:Horg+Hpsf-n, :, :) = Vk(Hpsf-n+1:Horg+Hpsf-n, :, :) + R * PSFh(n);
	Vk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) = Vk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) + R * PSFh(n);
end
R = Vk(Hpsf:Horg+Hpsf-1, :, :);
R(1:Hpsf-1, :, :) = R(1:Hpsf-1, :, :) + Vk(Hpsf-1:-1:1, :, :);
R(Horg-Hpsf+2:Horg, :, :) = R(Horg-Hpsf+2:Horg, :, :) + Vk(Horg+Hpsf*2-2:-1:Horg+Hpsf, :, :);
Vk = zeros(Horg, Worg, Sorg+Tpsf*2-2);
Vk(:, :, Tpsf:Sorg+Tpsf-1) = R * PSFt(1);
for n=2:Tpsf
	Vk(:, :, Tpsf-n+1:Sorg+Tpsf-n) = Vk(:, :, Tpsf-n+1:Sorg+Tpsf-n) + R * PSFt(n);
	Vk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) = Vk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) + R * PSFt(n);
end
R = Vk(:, :, Tpsf:Sorg+Tpsf-1);
R(:, :, 1:Tpsf-1) = R(:, :, 1:Tpsf-1) + Vk(:, :, Tpsf-1:-1:1);
R(:, :, Sorg-Tpsf+2:Sorg) = R(:, :, Sorg-Tpsf+2:Sorg) + Vk(:, :, Sorg+Tpsf*2-2:-1:Sorg+Tpsf);

% MSE
% ----------------------------------------------------------------
if exist('narrow') && ischar(narrow) && strcmp(narrow, 'narrow')
	G = G(Xmin : Xmax, Ymin : Ymax, Zmin : Zmax);
	R = R(Xmin : Xmax, Ymin : Ymax, Zmin : Zmax);
end
SE = (G - R) .^ 2;
MSE = mean(SE(:));
end