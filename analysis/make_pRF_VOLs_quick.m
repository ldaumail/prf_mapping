addpath('/Applications/freesurfer/7.4.1/matlab')

subj='sub-F019_230718'; % e.g. '026'

% measured 1/21/2020

ScreenVertRes = 768; % pixels
ScreenWidth = 17; %42.7; % cm
SceeenHorzRes = 1024; % pixles
ViewDistance = 48; %99; % cm

%ppd = (SceeenHorzRes/2)*(2*pi/360)*(1/atan(ScreenWidth/(2*ViewDistance)));
%%not needed but good to have the equation somewhere
ppc = SceeenHorzRes/ScreenWidth;
ScreenHeight = ScreenVertRes/ppc;%32; % cm

VisAngle=(2*atan2(ScreenHeight/2,ViewDistance))*(180/pi); % degrees
ppd=round(ScreenVertRes/VisAngle); %pixels per degree
%%%%%%%%%%%%%%%
volDir = ['/Users/tonglab/Documents/Loic/prfMapping/analysis/',subj,'/analysis/topUp/multibar_only/pRF_Vols'];
mkdir(volDir);

load(['/Users/tonglab/Documents/Loic/prfMapping/analysis/',subj,'/analysis/topUp/multibar_only/results_quick_multibar.mat']);


template=MRIread(['/Users/tonglab/Documents/Loic/prfMapping/analysis/' subj '/data/sess01/bold/template.nii.gz']);

ang_vol=template; ang_vol.vol=results_quick.ang;
MRIwrite(ang_vol,[volDir '/pRFquick_ang.nii.gz'])

ecc_vol=template; ecc_vol.vol=results_quick.ecc;
% rescale to degrees, and add cut-off 110% of radius visual angle
ecc_vol.vol=ecc_vol.vol./ppd;
ecc_vol.vol(ecc_vol.vol > (1.1*(ScreenVertRes/2)/ppd))=0;
MRIwrite(ecc_vol,[volDir '/pRFquick_ecc.nii.gz'])

expt_vol=template; expt_vol.vol=results_quick.expt;
MRIwrite(expt_vol,[volDir '/pRFquick_expt.nii.gz'])

rfsize_vol=template; rfsize_vol.vol=results_quick.rfsize;
% rescale to degrees, and add cut-off
rfsize_vol.vol=rfsize_vol.vol./ppd;
rfsize_vol.vol(rfsize_vol.vol > (1.1*(ScreenVertRes/2)/ppd))=0;
MRIwrite(rfsize_vol,[volDir '/pRFquick_rfsize.nii.gz'])

R2_vol=template; R2_vol.vol=results_quick.R2;
MRIwrite(R2_vol,[volDir '/pRFquick_R2.nii.gz'])

gain_vol=template; gain_vol.vol=results_quick.gain;
MRIwrite(gain_vol,[volDir '/pRFquick_gain.nii.gz'])

meanvol_vol=template; meanvol_vol.vol=results_quick.meanvol;
MRIwrite(meanvol_vol,[volDir '/pRFquick_meanvol.nii.gz'])
