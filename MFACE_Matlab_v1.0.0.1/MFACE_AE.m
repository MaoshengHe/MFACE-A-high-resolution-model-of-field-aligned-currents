%% This script estimates the AE index using the other inputs. An alternative approach is to build the model with the other inputs directly.
try
    ae_= isnan(ae);
    aegood=ae(~ae_);
end
%AE=120+MLT*0;
%thetaakasofu = acos(imfbz./imfBt);
%akasofu = vsw.*imfBscale.^2.*sin(thetakasofu/2).^4
if ~exist( 'bsint1h','var')||any(isnan(bsint1h)),
    try,
        bs_= isnan(bsint1h);
        bsgood=bsint1h(~bs_);
    end
    bsint1h=imfbz.*(imfbz<0);%bsint1h=-1.5+MLT*0;
    if exist( 'bs_','var')
        bsint1h(~bs_)=bsgood;
    end
end

bsint1h = bsint1h(:);
X0AE = [sin(doy_in_radian), cos(doy_in_radian), sin(doy_in_radian*2), cos(doy_in_radian*2), ...
    sin(clock_IMF), cos(clock_IMF), sin(clock_IMF*2), cos(clock_IMF*2), ...
    sin(clock_IMF*3), cos(clock_IMF*3), B_t, PotentialI, bsint1h]; % Akasofu
AEmodel = load('AEmodel.mat');
X2AE = x2fx(X0AE, AEmodel.modelab);
X2AEz =(X2AE-repmat(AEmodel.mX,[size(X0AE,1),1]))./repmat(AEmodel.sigmaX,[size(X0AE,1),1]);
X2AEz(:,1) = 1; % assign the 1st column to 1's
ae = X2AEz*AEmodel.b21;
if exist( 'ae_','var')
    ae(~ae_)=aegood;
end
%% <<< the end of the eastimation of AE index