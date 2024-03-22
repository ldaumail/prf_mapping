
addpath('/Applications/freesurfer/7.4.1/matlab')
addpath(genpath('/Users/tonglab/Documents/Loic/prfMapping/analysis/sub-F019_230718/code/utils/KK_tools_dev')) % source dev version
addpath('/Users/tonglab/Documents/Loic/prfMapping/analysis/sub-F019_230718/code/analysis')

subj='sub-F019_230718';
outpath=['/Users/tonglab/Documents/Loic/prfMapping/analysis/',subj,'/data/'];
TR=2;
conds=[93 93 93 93]; %, 94 94 94 94];  % 93 - multibar, 94 - wedgerings
runs=['001'; '002'; '003'; '004']; % '005'; '006'; '007'; '008']; % FSFAST directories, order matched condition list

%% Make stimuli
% showMasks=0; % set to zero when doing actual analysis on cluster, 1 to view masks locally
% tps=(7*TR):(15*TR):4500;



% load('/Users/tonglab/Documents/Loic/prfMapping/analysis/sub-F019_230718/code/utils/KK_tools_dev/stimuli.mat');
% 
% mystimuli=zeros([size(conds,2),768,768,size(tps,2)]);
% 
% % play stimuli
% for c=1:size(conds,2)
%     if showMasks
%         figure
%     end
%     for f=1:size(tps,2)
%         switch conds(c)
%             case 93
%                 if multibarindices(tps(f))
%                     if showMasks
%                         imshow(masks(:,:,multibarindices(tps(f))));
%                     end
%                     mystimuli(c,:,:,f)=masks(:,:,multibarindices(tps(f)));
%                 else
%                     if showMasks
%                         imshow(zeros(768,768));
%                     end
%                 end
%             case 94
%                 if wedgeringindices(tps(f))
%                     if showMasks
%                         imshow(masks(:,:,wedgeringindices(tps(f))));
%                     end
%                     mystimuli(c,:,:,f)=masks(:,:,wedgeringindices(tps(f)));
%                 else
%                     if showMasks
%                         imshow(zeros(768,768));
%                     end
%                 end
%         end
%         if showMasks
%             title(['Run' num2str(c) ' Cond' num2str(conds(c)) ': TR ' num2str(f) '/' num2str(size(tps,2))])
%             pause(0.5);
%         end
%     end
%     if showMasks
%         close all
%         disp('View next time-series?')
%         pause
%     end
% end
% % 
% % % Only want first two stimulus masks, assuming you are averaging BOLD runs
% % % within each condition (93-first and 94-second)
% finalstimuli={squeeze(mystimuli(1,:,:,:)), squeeze(mystimuli(2,:,:,:))};
% %Loic Daumail
% finalstimuli={squeeze(mystimuli(2,:,:,:))};

% stimulus and acquisition parameters are hardcoded
downSampleFactor = 1; % stimulus masks can be spatially downsampled
prfRadius = 7; % in degrees (saved out during stim pres)
ppd = 51; % pixels per degree (saved out during stim pres)
[multibarscreens, wedgescreens] = createPrfScreens(downSampleFactor, prfRadius, ppd);
finalstimuli={multibarscreens}; %wedgescreens};
%% Make concatenated and averaged BOLD data

% load ribbon label
ribbon_vol=MRIread([outpath 'sess01/bold/Yeo7-net1.VOL.nii.gz']);
ribbon_mask=find(ribbon_vol.vol == 1);

bar_runs=find(conds == 93);
wedge_runs=find(conds == 94);

% Cond 93
for r=1:size(bar_runs,2)
    currVol=MRIread([outpath 'sess01/bold/' runs(bar_runs(r),:) '/fmc.up.nii.gz']); %Loic 06/12/2023
    disp([outpath 'sess01/bold/' runs(bar_runs(r),:) '/fmc.up.nii.gz']) %Loic 06/12/2023
    if r == 1
        cond93Vol = zeros(currVol.height,currVol.width,currVol.depth,currVol.nframes);
    end
    cond93Vol=cond93Vol + currVol.vol;
end
final_cond93Vol=cond93Vol ./ size(bar_runs,2); % average across runs

% Cond 94
% for r=1:size(wedge_runs,2)
%     currVol=MRIread([outpath 'sess01/bold/' runs(wedge_runs(r),:) '/fmc.up.nii.gz']); %Loic 06/12/2023
%     disp([outpath 'sess01/bold/' runs(wedge_runs(r),:) '/fmc.up.nii.gz']) %Loic 06/12/2023
%     if r == 1
%         cond94Vol = zeros(currVol.height,currVol.width,currVol.depth,currVol.nframes); % init to data dimensions if first run
%     end
%     cond94Vol=cond94Vol + currVol.vol;
% end
% final_cond94Vol=cond94Vol ./ size(wedge_runs,2); % average across runs
% %interpDat= tseriesinterp(final_cond94Vol,2,1,4);
% % Concatenate [93 94] averaged runs, order matches finalstimuli
% finaldata=cat(4,final_cond93Vol,final_cond94Vol); 
finaldata={final_cond93Vol};%,final_cond94Vol}; %Loic 06/12/2023
%finaldata={interpDat}; %Loic 06/12/2023
%finaldata={final_cond94Vol};
%options.vxs=ribbon_mask;
options.display='off';  %Loic Daumail
options.seedmode=[-2]; %Loic Daumail

clearvars -except finalstimuli finaldata options outpath TR %avoid too much memory usage
results_quick=analyzePRF(finalstimuli,finaldata,TR,options);
%results = analyzePRF(stimulus,data,1,struct('seedmode',[0 1],'display','off'));
save([outpath 'results_quick_multibar.mat'], 'results_quick');
