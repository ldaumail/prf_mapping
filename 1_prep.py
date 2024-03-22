#!/usr/bin/python
'''
prepares directory structure, runs brain extractions, prepares b0 fieldmaps,
and all other steps necessary before FSL design files are created
'''

import os
import os.path as op
import glob
import datetime
import shutil

# get scan info from experiment file
#os.chdir('/Users/tong_processor/Desktop/Loic/retinotopy/dave_retinotopy/sub-F019/code/pipeline')
utils = '/Users/tonglab/Documents/Loic/prfMapping/analysis/sub-F019_230718/code/utils'
sys.path.append(op.expanduser(f'{utils}'))
from experiment import experiment
#from resize_topup_runs import resize_topup_runs

for fieldStrength in experiment['scanInfo'].keys():
    for subject in experiment['scanInfo'][fieldStrength].keys():
        for s, session in enumerate(experiment['scanInfo'][fieldStrength][subject].keys()):

            sessID = experiment['scanInfo'][fieldStrength][subject][session]['sessID']
            sessDir = os.path.join(f'/Users/tonglab/Documents/Loic/prfMapping/analysis/{subject}_{session}/data')
            rawDir = os.path.join(sessDir, 'sourceData/NIFTI')

            # copy over anatomical scan
            freesurferSubjectDir = os.path.join('/Users/tonglab/Documents/Loic/all_mri_data/freesurfer/', f'{subject}_{session}')
            # if not os.path.isdir(freesurferSubjectDir):
            anatScan = experiment['scanInfo'][fieldStrength][subject][session]['anatScan']
            inFile = glob.glob(os.path.join(rawDir, f'*{sessID}.{anatScan:02}*.nii'))[0]

            # print('Running surface segmentation...')
            # freesurferCommand = f'recon-all -subject {subject} -all -i {inFile}'
            #os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')
            # os.system(f'bash "{freesurferCommand}"')
            # make copies of anatomical and run brain extraction
            shutil.copy(inFile, os.path.join(freesurferSubjectDir, 'mri/orig/anatomical.nii'))
            os.system(f'bet {inFile} {os.path.join(freesurferSubjectDir, "mri/orig/anatomical_brain.nii.gz")}')

            # keep local copy of anat used for segmentation in subject dir
            anatDir = os.path.join(sessDir, 'anatomical')
            os.makedirs(anatDir, exist_ok=True)
            shutil.copy(os.path.join(freesurferSubjectDir, 'mri/orig/anatomical.nii'), os.path.join(anatDir, 'anatomical.nii') )
            shutil.copy(os.path.join(freesurferSubjectDir, 'mri/orig/anatomical_brain.nii.gz'), os.path.join(anatDir, 'anatomical_brain.nii.gz') )

            # copy over funcNoEPI scan
            outFile = None
            if 'funcNoEPIscan' in experiment['scanInfo'][fieldStrength][subject][session].keys():
                FNEdir = os.path.join(sessDir, 'funcNoEPI')
                os.makedirs(FNEdir, exist_ok=True)
                FNEscan = experiment['scanInfo'][fieldStrength][subject][session]['funcNoEPIscan']
                inFile = glob.glob(os.path.join(rawDir, f'*{sessID}.{FNEscan:02}*.nii'))[0]
                outFile = os.path.join(FNEdir, 'funcNoEPI.nii')
                if not os.path.exists(outFile):
                    shutil.copyfile(inFile, outFile)

            for scan in experiment['scanInfo'][fieldStrength][subject][session]['funcScans'].keys():

                for r, run in enumerate(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan]):

                    print(f'{datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")} | Field Strength: {fieldStrength} | Subject: {subject} | Session: {session} | Scan: {scan} | Run: {r+1}')
                    funcDir = os.path.join(sessDir, 'functional', scan, f'run{r+1:02}')
                    os.makedirs(os.path.join(funcDir, 'inputs'), exist_ok=True)

                    # move raw data
                    inFile = glob.glob(os.path.join(rawDir, f'*{sessID}.{run:02}*.nii'))[0]
                    outFile = os.path.join(funcDir, 'inputs/rawData.nii')
                    if not os.path.exists(outFile):
                        shutil.copyfile(inFile, outFile)

                    # move top up file
                    if fieldStrength == '7T':
                        inFile = glob.glob(os.path.join(rawDir, f'*{sessID}.{run+1:02}*.nii'))[0]
                        outFile = os.path.join(funcDir, 'inputs/oppPE.nii')
                        if not os.path.exists(outFile):
                            shutil.copyfile(inFile, outFile)
                        #change oppPE size to the corresponding functional run size
                        #resize_topup_runs(os.path.join(funcDir, 'inputs'))
                        # apply Top Up
                        inFile = os.path.join(funcDir, 'inputs/rawData.nii')
                        topupFile = os.path.join(funcDir, 'inputs/oppPE.nii')
                        outFile = os.path.join(funcDir, 'inputs/rawData_topUp.nii.gz')  # set final filename for top up output
                        if not os.path.isfile(outFile):
                            print('Running Top Up...')
                            os.system(f'python2 {utils}/fsl_TOPUP_call.py {inFile} {topupFile} 90 270')
                            toppedUpFile = glob.glob(f'{inFile[:-4]}*TUcorrected.nii.gz')[0]
                            os.rename(toppedUpFile, outFile)
                        extraFiles = glob.glob(os.path.join(funcDir, 'inputs/*topup*'))
                        for extrafile in extraFiles:
                            os.remove(extrafile)
                        extraFiles = glob.glob('*topup*')
                        for extrafile in extraFiles:
                            os.remove(extrafile)

                    # if 'b0Scan' in experiment['scanInfo'][fieldStrength][subject][session].keys():  # only run if b0 map exists
                    #
                    #     for topup, withWithout in zip(['_topUp', ''], ['with', 'without']):
                    #
                    #         inFile = glob.glob(os.path.join(funcDir, f'inputs/rawData{topup}.nii*'))[0]
                    #         outFile = os.path.join(funcDir, f'inputs/rawData{topup}_b0.nii.gz')  # set final filename for b0 output
                    #
                    #         if os.path.isfile(inFile) and not os.path.isfile(outFile):
                    #             print(f'b0 correcting scan {withWithout} topup')
                    #             os.system(f'fugue -i {inFile} --dwell={ees} --loadfmap={realFileCopy[:-7]}_rads_reg_unwrapped.nii.gz --unwarpdir=y- --asym={te} --despike -u {outFile}')

                # get ref image for preprocessing, register to surf, create cortical labels
                if scan == 'multibar':
                    nRuns = len(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan])
                    midRun = int(nRuns/2)
                    funcDir = os.path.join(sessDir, 'functional', scan, f'run{midRun:02}')

                    for topup, topupString in zip(['topUp', 'noTopUp'], ['_topUp', '']):
                        for b0, b0String in zip(['b0', 'noB0'], ['_b0', '']):
                            inFile = glob.glob(f'{funcDir}/inputs/rawData{topupString}{b0String}.nii*')
                            if len(inFile) == 1:
                                outDir = os.path.join(sessDir, f'functional/{scan}/allRuns', topup, b0, 'reg')
                                os.makedirs(outDir, exist_ok=True)

                                # get ref image
                                refFile = os.path.join(outDir, 'refFuncImage.nii.gz')
                                if not os.path.isfile(refFile):
                                    os.system(f'fslmaths {inFile[0]} -Tmean {refFile}')

                                fsAnatFileMGZ = f'{freesurferSubjectDir}/mri/orig.mgz'  # do not confuse with original nifti
                                fsAnatFileNII = f'{outDir}/anatomical_fs.nii'
                                if not os.path.isfile(fsAnatFileNII):
                                    freesurferCommand = f'mri_convert --in_type mgz --out_type nii --out_orientation RAS {fsAnatFileMGZ} {fsAnatFileNII}'
                                    os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')


                                # register freesurfer anat to orig anat (copy of anatomical saved by recon-all is not same as original nifti)
                                origAnatFile = f'{freesurferSubjectDir}/mri/orig/anatomical.nii'
                                fsAnat2origAnatMat = f'{outDir}/fsAnat2origAnat.mat'
                                if not os.path.isfile(fsAnat2origAnatMat):
                                    os.system(f'flirt -in {fsAnatFileNII} -ref {origAnatFile} -omat {fsAnat2origAnatMat}')

                                # make copy of orginal anats in reg dir
                                shutil.copy(origAnatFile, f'{outDir}/anatomical_orig.nii')
                                shutil.copy(f'{origAnatFile[:-4]}_brain.nii.gz', f'{outDir}/anatomical_orig_brain.nii.gz')

                                # register orig anat to func
                                if fieldStrength == '3T':
                                    origAnat4reg = f'{outDir}/anatomical_orig.nii'
                                elif fieldStrength == '7T':
                                    origAnat4reg = f'{outDir}/anatomical_orig_brain.nii.gz'

                                origAnat2funcMat = os.path.join(outDir, 'origAnat2func.mat')
                                if not os.path.isfile(origAnat2funcMat):
                                    os.system(f'flirt -in {origAnat4reg} -ref {refFile} -omat {origAnat2funcMat}')

                                # convert cortical mask to func space
                                for hemi in ['lh', 'rh']:

                                    # cortex mgz to nifti
                                    inFile = f'{freesurferSubjectDir}/mri/{hemi}.ribbon.mgz'
                                    outFile = f'{outDir}/cortex_{hemi}_FSanat.nii.gz'
                                    if not os.path.exists(outFile):
                                        freesurferCommand = f'mri_convert --in_type mgz --out_type nii --out_orientation RAS {inFile} {outFile}'
                                        os.system(f'bash {utils}/callFreesurferFunction.sh -s "{freesurferCommand}"')

                                # combine across hemis
                                lhFile = f'{outDir}/cortex_lh_FSanat.nii.gz'
                                rhFile = f'{outDir}/cortex_rh_FSanat.nii.gz'
                                biFile = f'{outDir}/cortex_bi_FSanat.nii.gz'
                                os.system(f'fslmaths {lhFile} -add {rhFile} -bin {biFile}')


                                # freesurfer anat space to orig anat space
                                inFile = biFile
                                outFile = f'{outDir}/cortex_bi_origAnat.nii.gz'
                                if not os.path.exists(outFile):
                                    os.system(f'flirt -in {inFile} -ref {origAnatFile} -applyxfm -init {fsAnat2origAnatMat} -interp nearestneighbour -out {outFile}')

                                # orig anat space to func space
                                inFile = outFile
                                outFile = f'{outDir}/cortex_bi_func.nii.gz'
                                if not os.path.exists(outFile):
                                    os.system(f'flirt -in {inFile} -ref {refFile} -applyxfm -init {origAnat2funcMat} -interp nearestneighbour -out {outFile}')

                                # dilate ROI for complete cortical coverage
                                inFile = outFile
                                outFile = f'{outDir}/finalAnalysisMask.nii.gz'
                                if not os.path.exists(outFile):
                                    os.system(f'fslmaths {inFile} -kernel box 6 -dilM -roi 0 -1 0 33 0 -1 0 1 {outFile}') # dilates and removes anterior voxels

                                # set up dirs for preprocessing
                                for r, run in enumerate(experiment['scanInfo'][fieldStrength][subject][session]['funcScans'][scan]):
                                    funcDir = os.path.join(sessDir, 'functional', scan, f'run{r + 1:02}')
                                    preprocOutDir = os.path.join(funcDir, 'outputs', topup, b0)
                                    os.makedirs(preprocOutDir, exist_ok=True)

print('FINISHED')