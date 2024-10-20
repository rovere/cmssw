#!/bin/bash

function die { echo $1: status $2; exit $2; }


# Check if a configuration name is provided as a command-line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <config_name>"
    echo "config_name: l1a140, l1a200, ttb200, ofp2l1a140, ofp2l1a200, or ofp2ttb200."
    exit 1
fi

# Get the configuration name from the command-line argument
config_name=$1

generate_hlt_cfg() {
  local input_file=$1
  local suffix=$2
  echo "Generating ${suffix} configuration"
  echo "Using ${input_file} as input"

  COMMON_CUSTOM='process.options.numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(64); process.options.numberOfConcurrentRuns = cms.untracked.uint32(64); process.options.wantSummary=True'

  echo "Adding the following customisations: ${COMMON_CUSTOM}"

  cmsDriver.py Phase2 -s L1P2GT,HLT:75e33_timing --processName=HLTX \
    --conditions auto:phase2_realistic_T33 --geometry Extended2026D110 \
    --era Phase2C17I13M9 \
    --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 \
    --eventcontent FEVTDEBUGHLT \
    --filein=file:$input_file \
    --mc --nThreads 4 --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
    -n 1 --no_exec --output={} \
    --python_filename Phase2_L1P2GT_HLT_${suffix}.py \
    --customise_commands="${COMMON_CUSTOM}" \
    --dump_python || die 'failed to produce a valid configuration...' $?
}

generate_ofp2_cfg() {
  local input_file=$1
  local suffix=$2
  echo "Generating ${suffix} configuration"

  COMMON_CUSTOM='process.options.numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(64); process.options.numberOfConcurrentRuns = cms.untracked.uint32(64); process.options.wantSummary=True;'

  echo "Adding the following customisations: ${COMMON_CUSTOM}"

  cmsDriver.py Phase2 -s RAW2DIGI,RECO --processName=OFFP2 \
    --conditions auto:phase2_realistic_T33 --geometry Extended2026D110 \
    --era Phase2C17I13M9 \
    --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 \
    --eventcontent RECO \
    --filein=file:$input_file \
    --mc --nThreads 4 --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
    -n 1 --no_exec \
    --python_filename Phase2_RAW2DIGI_RECO_${suffix}.py \
    --customise_commands="${COMMON_CUSTOM}" \
    --customise HLTrigger/Configuration/HLT_75e33/test/customize_genericConsumerOffline.customize_genericConsumerOffline,HLTrigger/Configuration/HLT_75e33/test/customize_removeMCFromOFP2.customize_removeMCFromOFP2 \
    --dump_python || die 'failed to produce a valid configuration...' $?
}

# Generate the specified configuration based on the command-line argument
case $config_name in
    l1a140)
      generate_hlt_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU140_Trk1GeV_140X_mcRun4_realistic_v4-v1/120000/output_Phase2_L1ASkimReducedMenu_Spring24_FixDiTau52_140PU_merged.root $config_name
        ;;
    l1a200)
      generate_hlt_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200ALCA_140X_mcRun4_realistic_v4-v2/120000/output_Phase2_L1ASkimReducedMenu_Spring24_FixDiTau52_200PU_merged.root $config_name
        ;;
    ttb200)
      generate_hlt_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/output_Phase2_L1T_Reduced.root $config_name
        ;;
    ofp2l1a140)
      generate_ofp2_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU140_Trk1GeV_140X_mcRun4_realistic_v4-v1/120000/output_Phase2_L1ASkimReducedMenu_Spring24_FixDiTau52_140PU_merged.root $config_name
        ;;
    ofp2l1a200)
      generate_ofp2_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/MinBias_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200ALCA_140X_mcRun4_realistic_v4-v2/120000/output_Phase2_L1ASkimReducedMenu_Spring24_FixDiTau52_200PU_merged.root $config_name
        ;;
    ofp2ttb200)
      generate_ofp2_cfg /data/user/rovere/store/mc/Phase2Spring24DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_AllTP_140X_mcRun4_realistic_v4-v1/2560000/output_Phase2_L1T_Reduced.root $config_name
        ;;
    *)
        echo "Invalid configuration name. Please choose l1a140, l1a200, ttb200, ofp2l1a140, ofp2l1a200 or ofp2ttb200."
        ;;
esac

echo "Script execution completed."

