#! /bin/bash

# Exit the script immediately if any command fails
set -e

# Enable pipefail to propagate the exit status of the entire pipeline
set -o pipefail

###MR FOLDER_FILES="/data/user/${USER}/"
###MR DATASET="/RelValTTbar_14TeV/CMSSW_14_1_0_pre6-PU_141X_mcRun4_realistic_v1_STD_2026D110_PU-v3/GEN-SIM-DIGI-RAW"
###MR FILES=( $(dasgoclient -query="file dataset=${DATASET}" --limit=-1 | sort | head -4) )
###MR 
###MR for f in ${FILES[@]}; do
###MR   # Create full MYPATH if it does not exist
###MR   MYPATH=$(dirname ${f})
###MR   if [ ! -d "${FOLDER_FILES}${MYPATH}" ]; then
###MR     echo "mkdir -p ${FOLDER_FILES}${MYPATH}"
###MR     mkdir -p ${FOLDER_FILES}${MYPATH}
###MR   fi
###MR   if [ -e "/eos/cms/${f}" ]; then
###MR     if [ ! -e "${FOLDER_FILES}${f}" ]; then
###MR       echo "cp /eos/cms/$f ${FOLDER_FILES}${MYPATH}"
###MR       cp /eos/cms/$f ${FOLDER_FILES}${MYPATH}
###MR     fi
###MR   fi
###MR done
 
###MR LOCALPATH=${FOLDER_FILES}$(dirname ${FILES[0]})
###MR LOCALPATH="/data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000"
###MR ###MR LOCALPATH="/data/user/rovere/store/relval/CMSSW_14_1_0_pre6/RelValMinBias_14TeV/GEN-SIM-DIGI-RAW/PU_141X_mcRun4_realistic_v1_STD_2026D110_PU-v3"
###MR echo "Local repository: |${LOCALPATH}|"
###MR LOCALFILES=$(ls -1 ${LOCALPATH})
###MR ALL_FILES=""
###MR for f in ${LOCALFILES[@]}; do
###MR   ALL_FILES+="file:${LOCALPATH}/${f},"
###MR done
###MR # Remove the last character
###MR ALL_FILES="${ALL_FILES%?}"
###MR echo "Discovered files: $ALL_FILES"
 
###MR cmsDriver.py Phase2 -s L1P2GT,HLT:75e33_timing --processName=HLTX \
###MR   --conditions auto:phase2_realistic_T33 --geometry Extended2026D110 \
###MR   --era Phase2C17I13M9 \
###MR   --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 \
###MR   --eventcontent FEVTDEBUGHLT \
###MR   --filein=${ALL_FILES} \
###MR   --mc --nThreads 8 --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
###MR   -n 1000 --no_exec --output={} --procModifier alpaka

if [ -e 'Phase2_L1P2GT_HLT_L1ASkimReducedMenu_200PU.py' ]; then
  if [ ! -d 'patatrack-scripts' ]; then
    git clone https://github.com/cms-patatrack/patatrack-scripts --depth 1
  fi
  patatrack-scripts/benchmark -j 32 -t 8 -s 8 -e 1000 --run-io-benchmark --event-skip 100 --event-resolution 25 -k Phase2Timing_resources.json -- Phase2_L1P2GT_HLT_L1ASkimReducedMenu_200PU.py
  mergeResourcesJson.py logs/step*/pid*/Phase2Timing_resources.json > Phase2Timing_resources.json
  if [ -e "$(dirname $0)/augmentResources.py" ]; then
    python3 $(dirname $0)/augmentResources.py
  fi
fi
