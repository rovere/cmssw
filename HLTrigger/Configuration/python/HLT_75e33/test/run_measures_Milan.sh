#! /bin/bash

# Exit the script immediately if any command fails
set -e

# Enable pipefail to propagate the exit status of the entire pipeline
set -o pipefail

function die { echo $1: status $2; exit $2; }

# hlt_config_names=("l1a140" "l1a200" "ttb200")
hlt_config_names=("l1a200" "ttb200")
ofp2_config_names=("ofp2l1a140" "ofp2l1a200" "ofp2ttb200")

# Generate each configuration by calling the original script in a loop
echo "Generating all HLT configurations..."

for config_name in "${hlt_config_names[@]}"; do
  jobs=32
  threads=8
  streams=8
  events=1000
  logdir="logs.Milan.HLT.$config_name.${jobs}j.${threads}t.${streams}s"
  cfg="Phase2_L1P2GT_HLT_${config_name}.py"
  # Generate the current configuration name
  source generate_cfgs.sh $config_name || die "Cannot generate $config_name configuration."

  # Download patatrack-scripts, if they are not there already.
  if [ ! -d 'patatrack-scripts' ]; then
    git clone https://github.com/cms-patatrack/patatrack-scripts --depth 1
  fi

  if [ ! -d "${logdir}" ]; then
    mkdir -p ${logdir}
  fi

  patatrack-scripts/benchmark -j ${jobs} -t ${threads} -s ${streams} -e ${events} --run-io-benchmark \
    -k Phase2Timing_resources.json --event-skip 100 --event-resolution 25 --wait 30 \
    --logdir ${logdir} -- ${cfg} | tee ${logdir}/output.log

  mergeResourcesJson.py ${logdir}/step*/pid*/Phase2Timing_resources.json > Phase2Timing_resources.json

  if [ -e "$(dirname $0)/augmentResources.py" ]; then
    python3 $(dirname $0)/augmentResources.py
  fi

  cp -p Phase2Timing_resources.json ${logdir}
  cp -p ${cfg} ${logdir}

  #Compress log files
  find ${logdir} -iname stderr -exec gzip {} \;
done

echo "All HLT configurations have been processed successfully."

for config_name in "${ofp2_config_names[@]}"; do
  jobs=4
  threads=8
  streams=8
  events=1000
  logdir="logs.Milan.OFP2.$config_name.${jobs}j.${threads}t.${streams}s"
  cfg="Phase2_RAW2DIGI_RECO_${config_name}.py"
  # Generate the current configuration name
  source generate_cfgs.sh $config_name || die "Cannot generate $config_name configuration."

  # Download patatrack-scripts, if they are not there already.
  if [ ! -d 'patatrack-scripts' ]; then
    git clone https://github.com/cms-patatrack/patatrack-scripts --depth 1
  fi

  if [ ! -d "${logdir}" ]; then
    mkdir -p ${logdir}
  fi

  patatrack-scripts/benchmark -j ${jobs} -t ${threads} -s ${streams} -e ${events} --run-io-benchmark \
    -k Phase2Timing_resources.json --event-skip 100 --event-resolution 25 --wait 30 \
    --logdir ${logdir} \
    --slot n=2:m=2:nv= \
    --slot n=2:m=2:nv= -- ${cfg} | tee ${logdir}/output.log

  mergeResourcesJson.py ${logdir}/step*/pid*/Phase2TimingOffline_resources.json > Phase2TimingOffline_resources.json
  cp -p Phase2TimingOffline_resources.json ${logdir}
  cp -p ${cfg} ${logdir}

  #Compress log files
  find ${logdir} -iname stderr -exec gzip {} \;
done

echo "All OFP2 configurations have been processed successfully."

