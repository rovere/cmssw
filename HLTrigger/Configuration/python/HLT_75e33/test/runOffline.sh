#! /bin/bash

# Exit the script immediately if any command fails
set -e

# Enable pipefail to propagate the exit status of the entire pipeline
set -o pipefail

###MR FOLDER_FILES="/data/user/${USER}/"
###MR #LOCALPATH="/data/user/rovere/store/relval/CMSSW_14_1_0_pre6/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_141X_mcRun4_realistic_v1_STD_2026D110_PU-v3/2810000/skimmedOffline"
###MR LOCALPATH=" /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000"
###MR echo "Local repository: |${LOCALPATH}|"
###MR LOCALFILES=$(ls -1 ${LOCALPATH})
###MR ALL_FILES=""
###MR for f in ${LOCALFILES[@]}; do
###MR   ALL_FILES+="file:${LOCALPATH}/${f},"
###MR done
###MR # Remove the last character
###MR ALL_FILES="${ALL_FILES%?}"
###MR echo "Discovered files: $ALL_FILES"

###MR cmsDriver.py Phase2 -s RAW2DIGI,RECO --processName=HLTX \
###MR   --conditions auto:phase2_realistic_T33 --geometry Extended2026D110 \
###MR   --era Phase2C17I13M9 \
###MR   --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 \
###MR   --eventcontent RECO \
###MR   --filein=${ALL_FILES} \
###MR   --mc --nThreads 28 --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
###MR   -n 10 --no_exec

if [ -e 'Phase2_RAW2DIGI_RECO_FullDump_L1ASkim200PU.py' ]; then
  if [ ! -d 'patatrack-scripts' ]; then
    git clone https://github.com/cms-patatrack/patatrack-scripts --depth 1
  fi
  patatrack-scripts/benchmark -j 2 -t 64 -s 64 -e 400 --run-io-benchmark --event-skip 100 --event-resolution 25 --slot c=0-63,128-191 -k Phase2TimingOffline_resources.json -- Phase2_RAW2DIGI_RECO_FullDump_L1ASkim200PU.py
  mergeResourcesJson.py logs/step*/pid*/Phase2TimingOffline_resources.json > Phase2TimingOffline_resources.json
#  if [ -e "$(dirname $0)/augmentResources.py" ]; then
#    python3 $(dirname $0)/augmentResources.py
#  fi
fi
