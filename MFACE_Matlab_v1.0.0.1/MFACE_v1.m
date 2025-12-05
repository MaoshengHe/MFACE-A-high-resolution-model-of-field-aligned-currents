function [J, ACC_MLat, EOF1, EOF2, EOF0] = MFACE_v1(MLT, MLAT, varargin)
% -------------------------------------------------------------------------------------------------------------------------------------------------------- 
% A high-resolution model of field-aligned currents through empirical orthogonal functions analysis (MFACE) 
% -------------------------------------------------------------------------------------------------------------------------------------------------------- 
% Version 1.0   
% MFACE is a high-resolution field-aligned currents model derived from ten years of CHAMP data with a novel 
% technique based on EOF decomposition. Predictor variables include MLat, MLT, DoY, IMF conditions, Vsw 
% and AE index. Compared to existing models, MFACE yields significantly better spatial resolution, reproduces 
% typically observed FAC thickness and intensity, improves on the magnetic local time distribution, and gives the 
% seasonal dependence of FAC latitudes and the NBZ current signature. Our EOF decomposition reveals that 
% EOF1 represents mainly the Bz-controlled large scale R1/R2 current pattern whereas EOF2 represents the By-
% controlled cusp current signature. 
% --------------------- 
% Reference
%  He, M., et al. (2012). "A high-resolution model of field-aligned currents through empirical orthogonal functions analysis (MFACE)." Geophys. Res. Lett. 39(18): L18105.
%  http://dx.doi.org/10.1029/2012GL053168
%------------------------
%Contents
%------------------------
% MFACE_v1.m                                            :  the code of MFACE  v1.0
% MFACE_30min_Lag_N_v1.0                        :  MFACE v1.0 coefficients for the Northern Hemisphere
% MFACE_30min_Lag_S_v1.0                        :  MFACE v1.0 coefficients for the Southern Hemisphere
% AEmodel.mat                                             : model coefficients for AE index
%--%
% MFACE_Example.m                                   : the code of one exmaple
% exampleData.mat                                       : the data for the exmaple
% MFACE_Example.png                                : the output of the exmaple
%--%
% FAC_in_the_*_Hemisphere_from_*_to_*.gif  : Movies for few days. 
% require_MFACE_output.txt                          : example message to require MFACE Movies or Figures
% £¡TO REQUIRE SIMILAR MOVIES OR FIGURES FOR ANY SPECIFIC PERIODS PLEASE CONTACT MAOSHENG HE.
%------------------------
%Communication to:  
%------------------------
% Maosheng He
% Jacobs University Bremen, 
% Campus Ring 1
% 28725 Bremen, Germany
% +49-421-200-3209
% m.he@jacobs-university.de  /  hmq512@gmail.com
% Sept. 20, 2013 
%%  
%-----------------------
% Syntax
%------------------------
% MFACE_v1(MLT,MLAT) returns the downward FAC density
% MFACE_v1(...,'OptionalVariableName',OptionalVariableValue,...)
%OptionalVariableName is  case-insensitive & oder-insensitive.
%OptionalVariableName list: doy, imfBx,imfBy,imfBz,Vsw,AE,BsInt1h
%------------------------
% Inputs
%------------------------
% MLT      :  Geomagnetic local time in hours
% MLAT    :  Geomagnetic latitude in degrees
% doy       :  Time in the format of  Day of Year in range 0 - 366, could be fraction. Default is 81.
% imfBy    :  IMF By component in GSM in nT, 30 min beforehand. Default is 3nT.
% imfBz    :  IMF Bz component in GSM in nT, 30 min beforehand. Default is 3nT.
% imfBx    :  IMF Bx component in GSM in nT, 30 min beforehand. Default is 3nT.
% Vsw      :  Solar wind velocity in km/s, 30 min beforehand. Default is 450km/s.
% AE        :  AE index, default value is modelled according to other solar winds conditions input
% BsInt1h   :  Integral of Bs (is Bz when Bz < 0, is 0 when Bz > 0) in the preceding hour.
%             Default value is estimated by assuming Bz in the preceding hour equals to the input instantaneous Bz.
%------------------------
% Outputs
%------------------------
% ACC_MLat :  Geomagnetic latitude at reference point of AAC in degrees
% J            :  Downward FAC density in uA/m2 at 110 km altitude
% EOF1     :  Current density in uA/m2 represented by EOF1 ( EOF1 gives the basic Region-1/Region-2 pattern varying mainly with the interplanetary magnetic fieldBz component)
% EOF2     :  Current density in uA/m2 represented byEOF2 ( EOF2 captures separately the cusp current signature and By-related variability)
% EOF0     :  Mean current density in uA/m2
%------------------------
%Examples
%------------------------
%---% 1st example 
% MLT = [12,13]; MLAT = [60 76]; doy = [321,1]; imfBx = 0*MLT; imfBy = 0*MLT; imfBz = 0*MLT; Vsw = 300*MLT;
% [J, ACC_MLat, EOF1,EOF2,EOF0] = MFACE_v1(MLT,MLAT,'doy',doy,'imfBy',imfBy,'imfBz',imfBz,'Vsw',Vsw) %,'imfBx',imfBx
%---% 2nd example
% clear;
% MLT = rand(1,1e5)*24; MLAT =(rand(1,1e5)-.5)*2*90 ; doy = 0*MLT+81; imfBx = 0*MLT; imfBy = 0*MLT; imfBz = 0*MLT-10; Vsw = 300+0*MLT;
% [J, ACC_MLat, EOF1,EOF2,EOF0] = MFACE_v1(MLT,MLAT,'doy',doy,'imfBy',imfBy,'imfBz',imfBz,'Vsw',Vsw); %,'imfBx',imfBx
% figure; scatter(MLT,MLAT,5,J,'Marker','x');colormap(BWRg);
% xlabel('MLT(hour)');ylabel('MLat(\circ)');set(gca,'xlim',[0 24],'ylim',[-1 1]*90)
% figure;  scatterOnPolar(MLT/24*360,MLAT,J,3,1);colorbar;%colormap(BWRg)   
% title('Downward FAC (\muAm^-^2)');%axis equal 
% try 
%     load exampleData.mat COLORMAP
%      colormap(COLORMAP);
% end
%---% 3rd example
% run MFACE_Example
%% 
% Copyright (c) 2012, Jacobs University Bremen 
% All rights reserved. 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
% following conditions are met: 
%  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the 
% following disclaimer. 
%  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the 
% following disclaimer in the documentation and/or other materials provided with the distribution. 
%  
% THIS SOFTWARE IS PROVIDED BY JACOBS UNIVERSITY BREMEN "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL JACOBS UNIVERSITY BREMEN BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

%%
%% (1) >>load variables
% step 1/3: load optional input variables
nin = length(varargin);
if nin > 1,
    for i = max(1,nin-1):-1:1,
        if ~ischar(varargin{i}), continue, end;
        varargin{i} = lower(varargin{i});
        if ~isnumeric(varargin{i+1}),error(['bad varar for ',varargin{i}]); end;
        expression = [varargin{i},'=varargin{i+1}(:);'];
        eval(expression);
        varargin(i+1) = [];
        nin = nin-1;
    end
end

% step 2/3:  set the missing variables to the default values.
if ~exist( 'doy','var'),
    doy = 81 + MLT*0;
end
if ~exist( 'imfbx','var'),
    imfbx = MLT*0 + 3;
end
if ~exist( 'imfby','var'),
    imfby = MLT*0 + 3;
end
if ~exist( 'imfbz','var'),
    imfbz = MLT*0 + 3;
end
if ~exist( 'vsw','var'),
    vsw = 450 + MLT*0;
end

% step 3/3: check the finite input values
doy(~isfinite(doy)) = 81;
imfbx(~isfinite(imfbx)) = 3;
imfby(~isfinite(imfby)) = 0;
imfbz(~isfinite(imfbz)) = 0;
vsw(~isfinite(vsw)) = 450;
%%  << end (1)
%% (2) >> Prepare regressors
imfB_scale = sqrt(imfby(:).^2 + imfbz(:).^2 + imfbx(:).^2); % IMF magnitude
SIZE = size(MLT);
MLT = MLT(:); MLAT = MLAT(:); doy = doy(:);
imfby = imfby(:); imfbz = imfbz(:);
imfB_scale = imfB_scale(:);
vsw = vsw(:);
MLT_in_radian = MLT/24*2*pi;
IMFB = imfby + imfbz*1i;
B_t = abs(IMFB); % Tangential IMF component in GSM y-z plane
clock_IMF = angle(IMFB); % IMF clock angle
elev_IMF = acos(imfbz./imfB_scale);% IMF elev. angle
PotentialI = vsw.^2*1e-4 + 11.7*imfB_scale.*sin(elev_IMF/2).^3;
doy_in_radian = (doy+10)/365.25*2*pi; doy_in_radian = mod(doy_in_radian,2*pi);
% % a seprate fucntion for bsint1h
%bsint1h=datebsint1h(dateN-40/24/60);
%% >>> the following  section is eastimating the AE index, which is active only AE is missing from the input or contains NaN 
if ~exist( 'ae','var')||any(isnan(ae(:))),
    run MFACE_AE
end
%%
ae= ae(:);
%{'SLT1';'CLT1';'SLT2';'CLT2';'SLT3';'CLT3';'Sdoy1';'Cdoy1';'Sdoy2';'Cdoy2';'Simf';'Cimf';'Simf2';'Cimf2';'imfBt';'AE1h';'Potential';}
X0 = [sin(MLT_in_radian), cos(MLT_in_radian), sin(MLT_in_radian*2), cos(MLT_in_radian*2), sin(MLT_in_radian*3), cos(MLT_in_radian*3), sin(MLT_in_radian*4), cos(MLT_in_radian*4), ...
    sin(doy_in_radian), cos(doy_in_radian), sin(doy_in_radian*2), cos(doy_in_radian*2), ...
    sin(clock_IMF), cos(clock_IMF), sin(clock_IMF*2), cos(clock_IMF*2), sin(clock_IMF*3),cos(clock_IMF*3), ...
    B_t, ae, PotentialI];
%%  << end (2)
%% (3) >> run MFACE
% northen Hemisphere
NS_ = MLAT >= 0;
if sum(NS_)>0,[J(NS_), ACC_MLat(NS_)] = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:));
    if nargout > 2,
        EOF1(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),1);
        EOF2(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),2);
        EOF0(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),0);
    end
end
% southern Hemisphere
NS_ = ~NS_;
if sum(NS_) > 0, [J(NS_), ACC_MLat(NS_)] = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:)); J(NS_) = -J(NS_);
    if nargout >2,
        EOF1(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),1); EOF1(NS_) = -EOF1(NS_);
        EOF2(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),2); EOF2(NS_) = -EOF2(NS_);
        EOF0(NS_) = MFACE_NS_JUB(MLAT(NS_),X0(NS_,:),0); EOF0(NS_) = -EOF0(NS_);
    end
end;
%%  << end (3)
J = reshape(J,SIZE);
ACC_MLat = reshape(ACC_MLat,SIZE);
if nargout >2,
    EOF1 = reshape(EOF1,SIZE); EOF2 = reshape(EOF2,SIZE); EOF0 = reshape(EOF0,SIZE);
end
%%
function [J, ACC_MLat] = MFACE_NS_JUB(MLAT, X0, EOFid)
%% inputs
% X0: indepedent regressors
% MLAT: geomagnetic latitude
% EOFid: id of EOF(s) taken into account.  Default is [0, 1, 2, ...12].
%% outputs
%J: Downward FAC density in uA/m2 at 110 km altitude
%ACC_MLat:  Geomagnetic latitude at reference point of AAC in degrees
%% load model coefficient
if all(MLAT>=0)
    FAC_= load('MFACE_30min_Lag_N_v1.0.mat');%     
elseif all(MLAT<=0)
    FAC_= load('MFACE_30min_Lag_S_v1.0.mat');
end
%--% list of coefficients
%     Regressors_Name_List: names of regressors %
%     coef_ACC_MLat: coef. of Auroral Oval Magnetic latitude % 
%     index_selected_regressors: index selecting the regressors % 
%     logical_index_organizing_regressors: index organizing the regressors %   
%     linear_index_organizing_regressors: %  
%     avarage_current: j, avarage current density profile %
%     EOFs:EOFs %
%     variance2all: contribution of PCs to the total variance %
%     coef_J_EOFs: coef. of Auroral Oval current density for all EOFs %
%     delta_MLat_grid: latitude to ACC %
%     average_X: the average of regressors %
%     std_X: stand deviation of regressors % 
%     Lag2SolarWind: time lag to solar winds on bow shock nose,   in minu
%     stats_Reg_EOF:  stats for the EOF regressions % 
%     stats_ACC_MLat:  stats for the  regression of Auroral Oval latitude % 
%--% note stats_ACC_MLat and each column of stats_Reg_EOF contains four elements, corresponding to, in order, 
%--% the R2 statistic, the F statistic and its p value, and an estimate of the error variance
%% 
J0_taken_into_account=1;
if nargin>2,
    J0_taken_into_account=0;
    if any(EOFid==0), J0_taken_into_account=1;EOFid(EOFid==0)=[];end
    EOF_single=FAC_.coef_J_EOFs(:,EOFid);
    FAC_.coef_J_EOFs=FAC_.coef_J_EOFs*0;
    FAC_.coef_J_EOFs(:,EOFid)=EOF_single;
    J0=0;
end
%% gernerate interaction and squared terms
X2 = x2fx(X0,FAC_.logical_index_organizing_regressors);
%% Standardized Regressors
X2z=(X2-repmat(FAC_.average_X,[size(X0,1),1]))./repmat(FAC_.std_X,[size(X0,1),1]);
X2z(:,1)=1;
%% model  ACC_MLat
ACC_MLat=X2z*FAC_.coef_ACC_MLat;
rMLAT=MLAT-ACC_MLat;
OKlat_=rMLAT>min(FAC_.delta_MLat_grid)&rMLAT<max(FAC_.delta_MLat_grid);
J(~OKlat_)=0;
if sum(OKlat_)>0,
    EOFs_interp=interp1(FAC_.delta_MLat_grid,FAC_.EOFs,rMLAT(OKlat_));    
    %% model  J
    Yr(OKlat_,:)=X2z(OKlat_,:)*FAC_.coef_J_EOFs;
    J(OKlat_)=sum(EOFs_interp.*Yr(OKlat_,:),2);
    if J0_taken_into_account,
        J0=interp1(FAC_.delta_MLat_grid,FAC_.avarage_current,rMLAT(OKlat_));
        J0=J0';
        J(OKlat_)=J(OKlat_)+J0;
    end
end