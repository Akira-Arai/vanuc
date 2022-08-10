# VANUC: Voxel-wise Anatomical-region-based Non-Uniformity Correction for partial volume effects

<img src="/fig1b.jpg" width="100%">  

[English Page](/README.md)  

# 概要
<img src="/VANUClogoM.jpg" width="15%">  

**VANUC** はVANUC法[^1]を用いてPET画像やSPECT画像の部分容積効果補正を行うプログラムです。  
[^1]:[Akira Arai, Yuriko Kagaya, Kentaro Takanami, Kei Takase. A novel partial volume effect correction method which allows for heterogeneous distribution: The potential advantages in the white matter activity estimation on FDG-PET. J Nucl Med. 2016;57(supplement 2):1924](https://jnm.snmjournals.org/content/57/supplement_2/1924)  

* 簡単なボタン操作のみで、画像データの読み込みから、位置合わせ、部分容積効果補正、統計解析に必要な標準脳座標画像の出力も可能  
* Hoffmann脳ファントム画像を入力すれば、処理に必要なパラメーターの自動算出も可能  

本プログラムはMATLAB上で動作し、SPMのインストールが必要    
（開発者：荒井 晃）  

# 目次
1. [使用前の準備](#1-使用前の準備)
    1. [使用環境](#1-1-使用環境)
    1. [VANUCのダウンロード](#1-2-vanucのダウンロード)
    1. [VANUCのパス設定](#1-3-vanucのパス設定)
    1. [画像の準備](#1-4-画像の準備)
1. [使用方法](#2-使用方法)
    1. [VANUCプログラムの実行](#2-1-vanucプログラムの実行)
    1. [出力データ一覧](#2-2-出力データ一覧)
1. [原理](#3-原理)
    1. [プログラムの処理フロー](#3-1-プログラムの処理フロー)
    1. [部分容積効果補正法の理論](#3-2-部分容積効果補正法の理論)
    1. [過去の検討結果](#3-3-過去の検討結果)
1. [ライセンス](#4-ライセンス)
1. [参考文献](#5-参考文献)

# 1. 使用前の準備
## 1-1. 使用環境
本ソフトウェアは数値計算ソフトウェアMATLAB上で動作します。  
* MATLABの購入は[こちらから](https://jp.mathworks.com/products/matlab.html)（外部サイト）  

SPM12の関数を使用するため、SPM12のダウンロードとパス設定が必要です。  
* SPM12のダウンロードは[こちらから](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)（外部サイト）  
* ダウンロード後に、MATLABでSPM12のパス設定を行って下さい。  
   1. MATLABを立ち上げる。  
   1. ツールバーの「環境」＞「パス設定」を開く。  
   1. 「フォルダーを追加...」をクリックし、「spm」フォルダーを選択し、「保存」をクリックする。  
   1. パス設定画面を閉じ、MATLABコマンドウィンドウに「spm」と入力し、SPM12が立ち上がることを確認する。  

## 1-2. VANUCのダウンロード
1. [こちらのページ](../../)のCodeをクリックし、Download ZIPをクリックしてリポジトリ内のファイルを一括ダウンロードする。  
1. ダウンロードしたzipファイルを解凍する。  
1. vanuc-mainフォルダーごと保存する。エラーを回避するため、C:ドライブの直下などへの保存が望ましい。  

## 1-3. VANUCのパス設定
1. MATLABを立ち上げる。  
1. ツールバーの「環境」＞「パス設定」を開く。  
1. 「フォルダーを追加...」をクリックし、「vanuc-main」フォルダー（うち最も最下位のフォルダー）を選択し、「保存」をクリックする。  
1. パス設定画面を閉じ、MATLABコマンドウィンドウに「vanuc」と入力し、vanucが立ち上がることを確認する。  

## 1-4. 画像の準備
* 部分容積効果補正処理など  
   * 脳PETまたはSPECTと、脳3D-T1強調画像が必要  
   * DICOM形式またはNIfTI形式に対応  
      * PET画像でのSUVへの変換は、DICOM形式の画像でのみ可能  
* ファントム画像
   * 現時点ではDICOM形式のみ対応  
* その他の形式、本ソフトウェアで読み込めない場合  
   * 他のソフトウェア等でNIfTI形式への変換を行って下さい。   

# 2. 使用方法
## 2-1. VANUCプログラムの実行
* MATLABコマンドウィンドウに「vanuc」と入力し、vanucのGUIを開く。  
* 処理画像などは、原則としてMATLABの「現在のフォルダー」に出力されるため、適宜変更。

### 部分容積効果補正ツール
* [VANUC法](#vanuc)、[Muller-Gartner法](#muller-gartner)、[RBV法](#rbv)による部分容積効果補正を行い、補正データや画像を出力。   
1. **「Start」** ボタンをクリック。  
1. **PETまたはSPECT画像を選択。** PETでSUV値に変換したい場合は「DICOM」データを読み込む。  
   * 「DICOM」を選択した場合：  
      1. DICOMファイルが保存されているフォルダーを選択。  
      1. フォルダー及びサブフォルダー内のDICOMデータのリストから、対象画像を選択。  

   * 「NIfTI」を選択した場合：
      1. NIfTIファイルを選択。
      1. Cropping（頭部以外の部分をカット）を行うかどうか選択。通常は「OK」を、既にcroppingが行われている場合には「Skip」を選択。  

1. **画像の方向を合わせる。** 画像の方向と、H（上）、F（下）、A（前）、P（後）、R（右）、L（左）の表示が、合っていれば「Confirm directions」をクリック。誤っていれば「Rotate/Flip」をクリックして、上方向、前方向、左右の順に表示に従って合わせる。  
1. **Cropping** （頭部以外の部分をカット）後の画像がカラー表示される。  
   * Croppingがうまくいかない場合の対処:   
      1. プログラムを中断。  
      1. 3.で出力されたNIfTIファイル（PETの場合はSUV_で始まるファイル）を、必要に応じて他のソフトウェアなどでトリミングを行う。  
      1. 再度処理を実行し、PETまたはSPECT画像として、作成した画像ファイルを選択して、CroppingをSkipする。

1. **PETパラメータの設定。** 予め設定するか（Preset）、画像から推定するか（Image-based）を選択する。通常は「Preset」が推奨される。  
   * 「Preset」を選択した場合：  
      * 最初に前回使用したパラメータ設定値（初回はデフォルト値）が表示される。変更が必要な場合は、予め保存したプロトコルから選択するか（Other Protocol）、直接入力するか（Direct Imput）を選ぶ。「Direct Imput」で入力した値を新たなプロトコルとして保存しておくと、次回から選択できる。また、ファントムデータからパラメータを決定したい場合は[下記](#ファントムデータ解析ツール)。  
      * 設定値の説明：  
         <details>
         <summary>FWHM of PSF: </summary>

         PET画像の点広がり関数の半値幅（mm）。（xは左右方向、yは前後方向、zは上下方向）  
         </details>
         <details>
         <summary>Distance from center of local region: </summary>
         VANUC法で用いられるカーネルの、中心ボクセルから辺縁までの距離（ボクセル単位）。通常はdefault値でよいが、分解能が低くノイズが多い場合には大きめに設定する。値が大きいと処理時間が長くなるので注意。  
         </details>

         <details>
         <summary>Regularization parameter: </summary>
         VANUC法で用いられる正則化パラメータ。臨床のPET・SPECT画像では1～100程度。ノイズが大きい場合には大きめの値とする。厳密にはL-curve解析などで決定する。  
         </details>

1. **3D T1強調画像** を選択し、画像の方向を合わせる。Croppingが行われる。（上記のPETまたはSPECTの場合と同様）  
1. 部分容積効果補正処理が始まり、10～30分で終わる。（途中、TPMファイル（組織確率マップ）がSPMのフォルダー内の所定の場所に見当たらない場合に、ファイルの場所を尋ねられる）  

### ファントムデータ解析ツール
* 部分容積効果補正に必要なPETパラメータを決定するためのツール。
* 原理については[下記](#vanucにおけるパラメーター決定)を参照。
1. 「Phantom」をクリック。  
1. ファントム画像（DICOM）が保存されているフォルダーを選択し、リストから対象のファントム画像を選択。  
1. ファントム画像の方向を合わせる。（上記と同様）  
1. デジタルファントムを選択。Hoffman 3D脳ファントムのデータは予め準備されている。  
1. 解析に数分かかる。  
1. FWHM of PSFの値が得られるので、他のパラメータやプロトコル名も入力して、新しいプロトコルとして保存。保存しない場合は「Cancel」をクリックして終了。保存すると、次回の部分容積効果補正から保存したプロトコルのパラメータ設定を使用できる。  

### その他の解析ツール
* dcm>nii: DICOMからNIfTIへの変換（PETの場合はSUV値への変換を含む）  
* Orientation: 画像の方向の変換（回転や左右反転など）  
* Crop: 頭部以外の部分のカット（先にOrientationの処理が行われる）  

## 2-2. 出力データ一覧
### 部分容積効果補正出力データ
<details>
<summary>現在のフォルダー内リスト（処理途中の画像や部分容積効果補正画像などの主な出力画像）：</summary>

|データ名|内容|
|:---|:---|
|PT~.niiまたはNM~.nii|PETまたはSPECT元画像|
|SUV_PT~.nii|SUV値への変換画像（PETの場合のみ）|
|trim_PET.nii|Cropping後のPET（またはSPECT）画像|
|MR~.nii|MRI元画像|
|T1.nii|Cropping後のMRI画像|
|coT1.nii|等方ボクセル化後のMRI画像|
|coT1_seg8.mat|SPMによるsegmentationデータ|
|c1~6coT1.nii|灰白質、白質、脳脊髄液、頭蓋骨、軟部組織、空気のマップ|
|y_coT1.nii|MRI座標からMNI標準脳座標への変換場|
|c1acoT1.nii|小脳灰白質マップ|
|c1bcoT1.nii|小脳以外の灰白質のマップ|
|rcoT1.nii|PET（またはSPECT）座標に位置合わせをしたMRI画像|
|rtrim_PET.nii|MRI座標に位置合わせをしたPET（またはSPECT）画像|
|PVCvanuc.nii|VANUC法による部分容積効果補正画像（7領域）|
|PVCmg.nii|Muller-Gartner法による部分容積効果補正画像（7領域）|
|PVCrbv.nii|RBV法による部分容積効果補正画像（7領域）|
|rPVC~.nii（3種）|MRI座標に位置合わせをした部分容積効果補正画像（7領域）|
|PVCC~.nii（3種）|全7領域を合成した部分容積効果補正画像|
|PVRC~.nii（3種）|全7領域を合成した部分容積効果"軽減"画像（FWHM 2.5mm）|
</details>

<details>
<summary>MNIspaceフォルダー内リスト（MNI標準脳座標に変換した画像。SPMによる統計画像解析などに利用。）：</summary>

|データ名|内容|
|:---|:---|
|wcoT1.nii|解剖学的標準化MRI画像|
|wrtrim_PET.nii|解剖学的標準化PET（またはSPECT）画像（PETはSUV値、SPECTは元の値）|
|wrPVC~.nii（3手法×3領域）|解剖学的標準化部分容積効果補正画像（小脳灰白質、その他の灰白質、白質）|
|wrPVC~\_cbl.nii（3手法×3領域）|小脳灰白質正規化解剖学的標準化部分容積効果補正画像|
|wrPVC~\_gm.nii（3手法×3領域）|灰白質正規化解剖学的標準化部分容積効果補正画像|
|wrPVC~\_wb.nii（3手法×3領域）|全脳正規化解剖学的標準化部分容積効果補正画像|
</details>

<details>
<summary>tempフォルダー内リスト（部分容積効果補正の生データ）：</summary>

|データ名|内容|
|:---|:---|
|Ppsf.mat|点広がり関数に関するパラメータ|
|Ploc.mat|VANUC法におけるカーネルサイズに関するパラメータ|
|Plam.mat|VANUC法における正則化に関するパラメータ|
|G.mat|PET（またはSPECT）画素値データ|
|M.mat|MRIの領域分割により得られた7領域のマップデータ|
|V.mat|7領域の領域広がり関数（RSF）|
|Rgtm.mat|GTM法による7領域の部分容積効果補正値|
|Rols.mat|VANUC法における7領域の期待値成分画像データ|
|Rvanuc.mat|VANUC法による7領域の部分容積効果補正画像データ|
|Rmg.mat|Muller-Gartner法による7領域の部分容積効果補正画像データ|
|Rrbv.mat|RBV法による7領域の部分容積効果補正画像データ|
|Voxel.mat|ボクセルサイズに関するパラメータ|
</details>

### ファントムデータ解析出力データ
（準備中）

# 3. 原理
## 3-1. プログラムの処理フロー
部分容積効果補正ツールの処理の流れを示す。
1. 元のPETまたはSPECT画像、3D T1強調画像、いくつかのパラメーターを入力。  
1. 3D T1強調画像のsegmentation  
1. 3D T1強調画像のMNI標準脳座標への解剖学的標準化  
1. 小脳マスクの個人脳座標への逆変換  
1. 得られた各組織領域画像（小脳灰白質、その他の灰白質、白質、脳脊髄液、頭蓋骨、軟部組織、空気）のPETまたはSPECT画像への位置合わせ  
1. 部分容積効果補正  
1. 各種画像の出力（出力画像については[上記](#部分容積効果補正出力データ)）  

## 3-2. 部分容積効果補正法の理論
多くの部分容積効果補正法は、次のような畳み込み積分を用いたモデルに基づいている。すなわち、PETまたはSPECTの画像 $G$ は、真の放射能濃度分布を $\rho$ 、PETまたはSPECTの点拡がり関数を $psf$ として、次式で表されるものとする。  

$$G = \rho \otimes psf$$  

このモデルに基づき、画像 $G$ から真の放射能濃度分布 $\rho$ を求めたいが、現実には $G$ はモデルからのずれや統計ノイズにより劣化しているため、単純に逆畳み込み積分により $\rho$ を推定することは難しい。そこで様々な手法が提案されてきたが、以下に紹介する手法では、MRI等の形態画像から得られる組織領域のマップ $M_k, (i=1,2,...,n)$ を用いることにより、画像劣化によって失われた組織領域の境界に関する情報を補っている。各組織領域のマップ $M_k$ は、ボクセル内に占める当該組織成分が占める体積割合を表したものであり、例えばSPMのsegmentationで得られるような画像である。また、多くの場合はそうであるように、VANUCプログラムでは、 $psf$ を3次元Gaussian関数とする。  

### GTM
GTM (geometry transfer matrix) 法は、Roussetらによって報告された手法である[^2]。GTM法では、各組織領域内の放射能濃度が一定であると仮定している。  
この仮定に基づいて、各領域内の放射能濃度から、各領域内の平均画素値への変換行列 $\mathbf{GTM}$ を求める。この変換行列を用いて、実際の画像から得られた平均画素値 $\mathbf{g}$ から逆変換によって真の放射能濃度 $\mathbf{\rho}$ を推定する。  
[^2]:[O G Rousset, Y Ma, A C Evans. Correction for partial volume effects in PET: principle and validation. J Nucl Med. 1998 May;39(5):904-11](https://pubmed.ncbi.nlm.nih.gov/9591599/)  

$$\mathbf{\hat{\rho_{GTM}}} = \mathbf{GTM^{-1}} \mathbf{g}, GTM_{i,j} = \dfrac{\displaystyle \sum M_i RSF_j}{\displaystyle \sum M_i}$$  

ここで $RSF$ は領域拡がり関数と呼ばれる。  

$$RSF_k = M_k \otimes psf$$  

### Muller-Gartner
Muller-Gartner法は、Muller-Gartnerによって報告された手法である[^3]。Muller-Gartner法も、組織内の放射能濃度が一定であると仮定する手法である。  
この仮定に基づいて、まず周囲の他の組織からのspill-inを取り除き、残りの画素値を領域拡がり関数で割ることにより、組織本来の放射能濃度を得る。  

$$\hat{\rho_{MG}} (k) = \dfrac{G - \displaystyle \sum_{i \neq k} C_i RSF_i}{RSF_k}$$  

ここで、周囲の他の組織の放射能濃度 $C_i$ を与えなければならないが、VANUCプログラムではGTM法で求めた値 $\hat{\rho_{GTM}}(i) $を使っている。  
[^3]:[H W Muller-Gartner, J M Links, J L Prince, R N Bryan, E McVeigh, J P Leal, C Davatzikos, J J Frost. Measurement of radiotracer concentration in brain gray matter using positron emission tomography: MRI-based correction for partial volume effects. J Cereb Blood Flow Metab. 1992 Jul;12(4):571-83](https://pubmed.ncbi.nlm.nih.gov/1618936/)  

### RBV
RBV (region-based voxel-wise correction) 法は、Thomasによって報告された手法である[^4]。RBV法も、組織内の放射能濃度が一定であると仮定する手法である。まず、GTM法で推定された放射能濃度分布を、これに $psf$ を畳み込み積分して得られた分布画像で割ることにより、ボクセルごとの補正係数を求める。このボクセルごとの補正係数を、元の画像 $G$ に掛けることによって、組織本来の放射能濃度を得る。  

$$\hat{\rho_{RBV}} = \dfrac{\displaystyle \sum_i^n \hat{\rho_{GTM}}(i) M_i}{\displaystyle \sum_i^n \hat{\rho_{GTM}}(i) RSF_i} G$$  

ここで、周囲の他の組織の放射能濃度 $C_i$ を与えなければならないが、VANUCプログラムではGTM法で求めた値を使っている。  
[^4]:[Benjamin A Thomas, Kjell Erlandsson, Marc Modat, Lennart Thurfjell, Rik Vandenberghe, Sebastien Ourselin, Brian F Hutton. The importance of appropriate partial volume correction for PET quantification in Alzheimer's disease. Eur J Nucl Med Mol Imaging. 2011 Jun;38(6):1104-19](https://pubmed.ncbi.nlm.nih.gov/21336694/)  

### VANUC
VANUC法[^1]は上述の手法とは異なり、組織内の放射能濃度の均一性を前提としない。その代わりに、注目するボクセルの近傍領域において、組織内の各ボクセルにおける放射能濃度はある決まった確率分布に従うものとする。部分容積効果によって複数の組織成分が混合した画素値 $G$ を、その確率分布に基づき、最も確率の高い割合で各組織成分に分配することによって、組織本来の放射能濃度を得る。確率分布が一定の期待値と分散をもつ正規分布 $N(\mu, \sigma^2)$ に従うものとみなすと、期待値成分 $\mu$ とそこからの偏差 $r$ との和として以下のように推定できる。  

$$\hat{\rho_{VANUC}}(k) = \hat{\mu}(k) + \hat{r}(k), $$  

$$\hat{r}(k) = \dfrac{\hat{{\sigma_k}^2} B_k}{\displaystyle \sum_t^n {\sigma_t}^2 B_t} \left(\dfrac{G}{RSF_k} - \hat{\mu}(k) \right), $$  

$$B_k = {M_k}^2 \otimes psf^2$$  

ここで期待値 $\mu$ は、近傍領域内の複数のボクセルデータを元に、一般化最小二乗法で推定されるが、VANUCプログラムでは簡単のために分散共分散行列 $\mathbf{\Omega}$ として $\sigma^2 \mathbf{I}$ を用いる。  

$$\begin{align}\mathbf{\hat{\mu}} &= (\mathbf{RSF}^T \mathbf{\Omega}^{-1} \mathbf{RSF} + \lambda^2 \mathbf{I})^{-1} \mathbf{RSF}^T \mathbf{\Omega}^{-1} \mathbf{G}\\
&= (\mathbf{RSF}^T \mathbf{RSF} + \lambda^2 \sigma^2 \mathbf{I})^{-1} \mathbf{RSF}^T \mathbf{G} \end{align}$$  

また、確率分布の分散 $\sigma^2$ は、組織領域内において一定であるとみなし、上述のGTM法と似た原理で、画像 $G$ での各領域内の画素値の分散 $s^2$ から推定する。 

$$\mathbf{\hat{V}} = \mathbf{A^{-1}} \mathbf{S}, A_{i,j} = \dfrac{\displaystyle \sum M_i B_j}{\displaystyle \sum M_i}, $$  

$$\mathbf{V} = \begin{pmatrix} {\sigma_1}^2\\
\vdots\\
{\sigma_n}^2 \end{pmatrix}, \mathbf{S} = \begin{pmatrix} {s_1}^2\\
\vdots\\
{s_n}^2 \end{pmatrix}$$  

### 部分容積効果軽減画像
VANUCプログラムでは、通常の部分容積効果補正画像に加えて、部分容積効果"軽減"画像を出力する。これは、部分容積効果補正の原理を応用して得られる、特定の分解能（本ソフトウェアではFWHM = 2.5mm）を持つ画像である。これによって、より分解能の低い元のPETまたはSPECT画像の部分容積効果を低減するだけでなく、分解能の異なる様々な画像から仮想的に等しい分解能を持つ画像を得ることができる。部分容積効果"低減"は、単に部分容積効果補正画像の平滑化ではない。  
具体的には、Muller-Gartner法やRBV法に基づく部分容積効果"低減"画像 $H$ は、設定された分解能の領域拡がり関数 $P$ を用いて以下のようにして得られる。  

$$H(k) = \hat{\rho}(k) P_k$$  

一方VANUC法の場合には、期待値成分に対して同様の演算を行う。  

$$H(k) = \hat{\mu}(k) P_k + \hat{r}(k)$$  

### VANUCにおけるパラメーター決定
VANUCソフトウェアは、過去の検討結果に基づき、Hoffman脳ファントムを用いた分解能測定を推奨するが、他のファントムに対しても応用可能である。  
解析フローは以下の通り。  
1. ファントムを撮像して得られた実際のPETまたはSPECT画像に、ファントムの真の放射能濃度分布画像の位置を合わせる。  
1. 真の分布画像に、一定の半値幅をもった3次元Gaussian関数を畳み込み積分してモデル画像を作成する。  
1. モデル画像と実際のPETまたはSPECT画像との差を、平均二乗誤差MSEとして評価する。  
1. 2.～3.を反復して、MSEを最小化するような半値幅を探索する。  

## 3-3. 過去の検討結果
（準備中）  

# 4. 使用上の注意
本ソフトウェアを使用して得られた成果を、科学雑誌や学会等で発表する際には、ソフトウェア名などを明記してください。  
開発者は、本ソフトウェアの使用や不具合によって生じた損害に対する責任を負いません。  
ライセンス条項については[こちら](/LICENSE)をご確認下さい。  

# 5. 参考文献
