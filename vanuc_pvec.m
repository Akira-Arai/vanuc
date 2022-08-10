function varargout = vanuc_pvec(varargin)
% Main part of partial volume effect correction using VANUC method
% 
% This function performs partial volume effect correction of image 
% (G.mat file) using morphological information (M.mat file) and 
% parameters (Ppsf.mat, Ploc.mat, Plam.mat files). 
% Three methods of PVEC are used: VANUC, Muller-Gartner 
% (Muller-Gartner, 1992) and RBV (Thomas et al, 2011).
% VANUC was developed by Akira Arai (Kousei Sendai Clinic), and 
% is based on probability theory assuming a non-uniformity. 
% All these methods also use the GTM method (Rousset et al, 1998) 
% in the process. 
% 
% Referrence:
%  1. H W Muller-Gartner, J M Links, J L Prince, R N Bryan, 
%     E McVeigh, J P Leal, C Davatzikos, J J Frost. Measurement of 
%     radiotracer concentration in brain gray matter using positron 
%     emission tomography: MRI-based correction for partial volume 
%     effects. J Cereb Blood Flow Metab. 1992 Jul;12(4):571-83
%  2. Benjamin A Thomas, Kjell Erlandsson, Marc Modat, 
%     Lennart Thurfjell, Rik Vandenberghe, Sebastien Ourselin, 
%     Brian F Hutton. The importance of appropriate partial volume 
%     correction for PET quantification in Alzheimer's disease. 
%     Eur J Nucl Med Mol Imaging. 2011 Jun;38(6):1104-19
%  3. Akira Arai, Yuriko Kagaya, Kentaro Takanami, Kei Takase.
%     A novel partial volume effect correction method which
%     allows for heterogeneous distribution: The potential
%     advantages in the white matter activity estimation on
%     FDG-PET. J Nucl Med. 2016 May 1;57(supplement 2):1924
%  4. O G Rousset, Y Ma, A C Evans. Correction for partial volume 
%     effects in PET: principle and validation. J Nucl Med. 
%     1998 May;39(5):904-11
% ----------------------------------------------------------------
% Akira Arai (Kousei Sendai Clinic)

% Import tissue probability maps (M)
% ----------------------------------------------------------------
cd temp
load 'M.mat'
sizeM = size(M);
Horg = sizeM(1);
Worg = sizeM(2);
Sorg = sizeM(3);
Kseg = sizeM(4);
clear sizeM
for a=1:Kseg
	NAMEmk = strcat('M', num2str(a), '.mat');
	Mk = M(:,:,:,a);
	save Mk
	movefile('Mk.mat', NAMEmk);
end
clear M a Mk NAMEmk
disp(' File Import D O N E')

% Calculate point spread function (PSF)
% ----------------------------------------------------------------
load 'Ppsf.mat'
clear Fh Ft Fw
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
save 'PSF.mat' PSFw Wpsf PSFh Hpsf PSFt Tpsf
clear Sh Sw St PSFw Wpsf PSFh Hpsf PSFt Tpsf
disp(' Point Spread Function Calculation D O N E')

% Convolution
% ----------------------------------------------------------------
load 'PSF.mat'
V = zeros(Horg, Worg, Sorg, Kseg);
for k=1:Kseg
	NAMEmk = strcat('M', num2str(k), '.mat');
	load (NAMEmk)
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
clear k NAMEmk Mk Vk PSFw PSFh PSFt Wpsf Hpsf Tpsf
save 'V.mat' V
clear V
disp(' Region Spread Function Calculation D O N E')

% Geometry transfer matrix (GTM)
% ----------------------------------------------------------------
VM = zeros(Kseg);
GM = zeros(Kseg, 1);
load 'V.mat'
load 'G.mat'
for k=1:Kseg
	NAMEmk = strcat('M', num2str(k), '.mat');
	load (NAMEmk)
	VM(k, :) = sum(reshape(Mk.*V, [Worg*Horg*Sorg, Kseg]));
	MkG = Mk.*G;
	GM(k) = sum(MkG(:));
end
clear NAMEmk Mk MkG
Rgtm = VM \ GM;
for k=1:Kseg
	G = G - Rgtm(k) * V(:, :, :, k);
end
clear V k GM VM
G = G.^2;
Var = mean(G(:));
clear G
save 'Rgtm.mat' Rgtm
disp(' GTM D O N E');

load M.mat
RM = sum(reshape(Rgtm, [1 1 1 Kseg]) .* M, 4);
fim = vanuc_disp(RM, 1, Vh, Vw, Vt);
fim.Name = 'PVEC by GTM';
pause(0.1);

clear Rgtm M

% Ordinary least squares (OLS)
% ----------------------------------------------------------------
% Extension
load 'Ploc.mat'
load 'V.mat'
Vh = zeros(Horg+Lh*2, Worg, Sorg, Kseg);
Vh(Lh+1:Horg+Lh, :, :, :) = V;
Vh(1:Lh, :, :, :) = V(Lh:-1:1, :, :, :);
Vh(Horg+Lh+1:Horg+Lh*2, :, :, :) = V(Horg:-1:Horg-Lh+1, :, :, :);
clear V
Vhw = zeros(Horg+Lh*2, Worg+Lw*2, Sorg, Kseg);
Vhw(:, Lw+1:Worg+Lw, :, :) = Vh;
Vhw(:, 1:Lw, :, :) = Vh(:, Lw:-1:1, :, :);
Vhw(:, Worg+Lw+1:Worg+Lw*2, :, :) = Vh(:, Worg:-1:Worg-Lw+1, :, :);
clear Vh Vw Vt
Vhws = zeros(Horg+Lh*2, Worg+Lw*2, Sorg+Lt*2, Kseg);
Vhws(:, :, Lt+1:Sorg+Lt, :) = Vhw;
Vhws(:, :, 1:Lt, :) = Vhw(:, :, Lt:-1:1, :);
Vhws(:, :, Sorg+Lt+1:Sorg+Lt*2, :) = Vhw(:, :, Sorg:-1:Sorg-Lt+1, :);
clear Vhw
load 'G.mat'
Gh = zeros(Horg+Lh*2, Worg, Sorg);
Gh(Lh+1:Horg+Lh, :, :) = G;
Gh(1:Lh, :, :) = G(Lh:-1:1, :, :);
Gh(Horg+Lh+1:Horg+Lh*2, :, :) = G(Horg:-1:Horg-Lh+1, :, :);
clear G
Ghw = zeros(Horg+Lh*2, Worg+Lw*2, Sorg);
Ghw(:, Lw+1:Worg+Lw, :) = Gh;
Ghw(:, 1:Lw, :) = Gh(:, Lw:-1:1, :);
Ghw(:, Worg+Lw+1:Worg+Lw*2, :) = Gh(:, Worg:-1:Worg-Lw+1, :);
clear Gh
Ghws = zeros(Horg+Lh*2, Worg+Lw*2, Sorg+Lt*2);
Ghws(:, :, Lt+1:Sorg+Lt) = Ghw;
Ghws(:, :, 1:Lt) = Ghw(:, :, Lt:-1:1);
Ghws(:, :, Sorg+Lt+1:Sorg+Lt*2) = Ghw(:, :, Sorg:-1:Sorg-Lt+1);
clear Ghw

% OLS
load 'Plam.mat'
lamVE = (Lam ^ 2) * (Var ^ 2) * eye(Kseg);
clear Lam Var
Rols = zeros(Horg, Worg, Sorg, Kseg);
llw = 2 * Lw;
llh = 2 * Lh;
llt = 2 * Lt;
for t=1:Sorg
	for h=1:Horg
		for w=1:Worg
			Gloc = Ghws(h:h+llh, w:w+llw, t:t+llt);
			Gloc = Gloc(:);
			Vloc = reshape(Vhws(h:h+llh, w:w+llw, t:t+llt, :), [Lsize, Kseg]);
			Rols(h,w,t,:) = reshape((Vloc' * Vloc + lamVE) \ (Vloc' * Gloc), [1, 1, 1, Kseg]);
		end
	end
end
clear Lh Lw Lt Lsize lamVE Vhws Ghws llw llh llt t h w Gloc Vloc
save 'Rols.mat' Rols
clear Rols
disp(' OLS D O N E')

% Calcurate region squared function
% ----------------------------------------------------------------
load 'PSF.mat'
delete 'PSF.mat'
PSFw = PSFw .^ 2;
PSFw = PSFw / (sum(PSFw) * 2 - PSFw(1));
PSFh = PSFh .^ 2;
PSFh = PSFh / (sum(PSFh) * 2 - PSFh(1));
PSFt = PSFt .^ 2;
PSFt = PSFt / (sum(PSFt) * 2 - PSFt(1));
W = zeros(Horg, Worg, Sorg, Kseg);
for k=1:Kseg
	NAMEmk = strcat('M', num2str(k), '.mat');
	load (NAMEmk)
	Wk = zeros(Horg, Worg+Wpsf*2-2, Sorg);
	Wk(:, Wpsf:Worg+Wpsf-1, :) = Mk * PSFw(1);
	for n=2:Wpsf
		Wk(:, Wpsf-n+1:Worg+Wpsf-n, :) = Wk(:, Wpsf-n+1:Worg+Wpsf-n, :) + Mk * PSFw(n);
		Wk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) = Wk(:, Wpsf+n-1:Worg+Wpsf+n-2, :) + Mk * PSFw(n);
	end
	Mk = Wk(:, Wpsf:Worg+Wpsf-1, :);
	Mk(:, 1:Wpsf-1, :) = Mk(:, 1:Wpsf-1, :) + Wk(:, Wpsf-1:-1:1, :);
	Mk(:, Worg-Wpsf+2:Worg, :) = Mk(:, Worg-Wpsf+2:Worg, :) + Wk(:, Worg+Wpsf*2-2:-1:Worg+Wpsf, :);
	Wk = zeros(Horg+Hpsf*2-2, Worg, Sorg);
	Wk(Hpsf:Horg+Hpsf-1, :, :) = Mk * PSFh(1);
	for n=2:Hpsf
		Wk(Hpsf-n+1:Horg+Hpsf-n, :, :) = Wk(Hpsf-n+1:Horg+Hpsf-n, :, :) + Mk * PSFh(n);
		Wk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) = Wk(Hpsf+n-1:Horg+Hpsf+n-2, :, :) + Mk * PSFh(n);
	end
	Mk = Wk(Hpsf:Horg+Hpsf-1, :, :);
	Mk(1:Hpsf-1, :, :) = Mk(1:Hpsf-1, :, :) + Wk(Hpsf-1:-1:1, :, :);
	Mk(Horg-Hpsf+2:Horg, :, :) = Mk(Horg-Hpsf+2:Horg, :, :) + Wk(Horg+Hpsf*2-2:-1:Horg+Hpsf, :, :);
	Wk = zeros(Horg, Worg, Sorg+Tpsf*2-2);
	Wk(:, :, Tpsf:Sorg+Tpsf-1) = Mk * PSFt(1);
	for n=2:Tpsf
		Wk(:, :, Tpsf-n+1:Sorg+Tpsf-n) = Wk(:, :, Tpsf-n+1:Sorg+Tpsf-n) + Mk * PSFt(n);
		Wk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) = Wk(:, :, Tpsf+n-1:Sorg+Tpsf+n-2) + Mk * PSFt(n);
	end
	Mk = Wk(:, :, Tpsf:Sorg+Tpsf-1);
	Mk(:, :, 1:Tpsf-1) = Mk(:, :, 1:Tpsf-1) + Wk(:, :, Tpsf-1:-1:1);
	Mk(:, :, Sorg-Tpsf+2:Sorg) = Mk(:, :, Sorg-Tpsf+2:Sorg) + Wk(:, :, Sorg+Tpsf*2-2:-1:Sorg+Tpsf);
	W(:, :, :, k ) = Mk;
end
clear k n NAMEmk Mk Wk PSFw PSFh PSFt Wpsf Hpsf Tpsf
save 'W.mat' W
clear W
disp(' Region Squared Function Calculation D O N E')

% Calcurate regional variances
% ----------------------------------------------------------------
load 'Rols.mat'
load 'V.mat'
E = zeros(Horg, Worg, Sorg);
for k=1:Kseg
	E = E + Rols(:, :, :, k) .* V(:, :, :, k);
end
clear Rols V
load 'G.mat'
E = (G - E) .^ 2;
clear G
load 'W.mat'
WM = zeros(Kseg);
EM = zeros(Kseg, 1);
for k=1:Kseg
	NAMEmk = strcat('M', num2str(k), '.mat');
	load(NAMEmk);
	delete(NAMEmk);
	WM(k, :) = sum(reshape(Mk .* W, [Worg*Horg*Sorg, Kseg]));
	MkE = Mk .* E;
	EM(k) = sum(MkE(:));
end
clear E W k NAMEmk Mk MkE
Vols = WM \ EM;
clear WM EM
disp(' Variance (OLS) Calculation D O N E')

% Maximum likelihood estimation (MLE)
% ----------------------------------------------------------------
load 'Rols.mat'
load 'V.mat'
RV = Rols .* V;
clear Rols
E = sum(RV, 4);
clear RV
load 'G.mat'
E = G - E;
clear G
load 'W.mat'
delete 'W.mat'
Sd = reshape(Vols, [1, 1, 1, Kseg]) .* W;
Rd = reshape(Vols, [1, 1, 1, Kseg]) ./ sum(Sd, 4) .* E;
clear Vols E Sd V W
Rd(find(isfinite(Rd)==0))=0;
load 'Rols.mat'
Rvanuc = Rols + Rd;
save 'Rvanuc.mat' Rvanuc
clear Rols Rd Rvanuc
disp(' MLE (OLS) D O N E')

% modified Muller-Gartner method (mMG)
% ----------------------------------------------------------------
load 'G.mat'
load 'Rgtm.mat'
load 'V.mat'
E = G;
clear G
for k=1:Kseg
	E = E - Rgtm(k) * V(:, :, :, k);
end
Rmg = V;
clear V
for k=1:Kseg
	Rmg(:, :, :, k ) = E ./ Rmg(:, :, :, k) + Rgtm(k);
end
Rmg(find(isfinite(Rmg)==0))=0;
save 'Rmg.mat' Rmg
clear Rmg
disp(' mMG D O N E')

% region-based voxel-wise method (RBV)
% ----------------------------------------------------------------
load 'G.mat'
R = G ./ (G - E);
clear G E
R(find(isfinite(R)==0))=0;
Rrbv = zeros(Horg, Worg, Sorg, Kseg);
for k=1:Kseg
	Rrbv(:,:,:,k ) = R .* Rgtm(k);
end
save 'Rrbv.mat' Rrbv
clear k R Rrbv Rgtm
disp(' RBV D O N E')
cd ..
close(fim);
disp('PVEC completed')

end