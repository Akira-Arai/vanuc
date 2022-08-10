function [MSE, Rgtm, GM] = vanuc_GTMMSE(G, M, Sh, Sw, St)
% Estimate difference of GTM-based model from the observed image
% 
% Input:
% G (3D double): Observed image (PET or SPECT)
% M (4D double): Tissue maps (multiple segments)
% Sh (double): x-sigma
% Sw (double): y-sigma
% St (double): z-sigma
% 
% Return:
% MSE (double): MSE of GTM-based model from the observed image
% Rgtm (1xn double): GTM-corrected values of n tissues
% GM (1xn double): uncorrected values of n tissues
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

sizeM = size(M);
Kseg = sizeM(4);
Worg = sizeM(2);
Horg = sizeM(1);
Sorg = sizeM(3);

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
V = zeros(Horg, Worg, Sorg, Kseg);
for k=1:Kseg
	Mk = M(:, :, :, k);
	Vk = zeros(Horg, Worg+Wpsf*2-2, Sorg);
	Vk(:, Wpsf:Worg+Wpsf-1, :) = Mk * PSFw(1);
	for n=2:Wpsf
		Vk(:, Wpsf-n+1:Worg+Wpsf-n, :) = Vk(:, Wpsf-n+1:Worg+Wpsf-n, :) + Mk * PSFw(n);
		Vk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) = Vk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) + Mk * PSFw(n);
	end
	Mk = Vk(:, Wpsf:Worg+Wpsf-1, :);
	Mk(:, 1:Wpsf-1, :) = Mk(:, 1:Wpsf-1, :) + Vk(:, Wpsf-1:-1:1, :);
	Mk(:, Worg-Wpsf+2:Worg, :) = Mk(:, Worg-Wpsf+2:Worg, :) + Vk(:, Worg+Wpsf*2-2:-1:Worg+Wpsf, :);
	Vk = zeros(Horg+Hpsf*2-2, Worg, Sorg);
	Vk(Hpsf:Horg+Hpsf-1, :, :) = Mk * PSFh(1);
	for n=2:Hpsf
		Vk(Hpsf-n+1:Horg+Hpsf-n, :, :) = Vk(Hpsf-n+1:Horg+Hpsf-n, :, :) + Mk * PSFh(n);
		Vk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) = Vk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) + Mk * PSFh(n);
	end
	Mk = Vk(Hpsf:Horg+Hpsf-1, :, :);
	Mk(1:Hpsf-1, :, :) = Mk(1:Hpsf-1, :, :) + Vk(Hpsf-1:-1:1, :, :);
	Mk(Horg-Hpsf+2:Horg, :, :) = Mk(Horg-Hpsf+2:Horg, :, :) + Vk(Horg+Hpsf*2-2:-1:Horg+Hpsf, :, :);
	Vk = zeros(Horg, Worg, Sorg+Tpsf*2-2);
	Vk(:, :, Tpsf:Sorg+Tpsf-1) = Mk * PSFt(1);
	for n=2:Tpsf
		Vk(:, :, Tpsf-n+1:Sorg+Tpsf-n) = Vk(:, :, Tpsf-n+1:Sorg+Tpsf-n) + Mk * PSFt(n);
		Vk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) = Vk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) + Mk * PSFt(n);
	end
	Mk = Vk(:, :, Tpsf:Sorg+Tpsf-1);
	Mk(:, :, 1:Tpsf-1) = Mk(:, :, 1:Tpsf-1) + Vk(:, :, Tpsf-1:-1:1);
	Mk(:, :, Sorg-Tpsf+2:Sorg) = Mk(:, :, Sorg-Tpsf+2:Sorg) + Vk(:, :, Sorg+Tpsf*2-2:-1:Sorg+Tpsf);
	V(:, :, :, k ) = Mk;
end
clear Mk

% Geometry transfer matrix (GTM)
% ----------------------------------------------------------------
VM = zeros(Kseg);
GM = zeros(Kseg, 1);
for k=1:Kseg
	VM(k, :) = sum(reshape(M(:,:,:,k).*V, [Worg*Horg*Sorg, Kseg]));
	MkG = M(:,:,:,k).*G;
	GM(k) = sum(MkG(:));
end
clear M MkG
Rgtm = VM \ GM;

% MSE
% ----------------------------------------------------------------
for k=1:Kseg
	G = G - Rgtm(k) * V(:, :, :, k);
end
clear V
G = G.^2;
MSE = mean(G(:));
end