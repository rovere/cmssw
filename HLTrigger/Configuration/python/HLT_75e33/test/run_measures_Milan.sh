#! /bin/bash

# Exit the script immediately if any command fails
set -e

# Enable pipefail to propagate the exit status of the entire pipeline
set -o pipefail

function die { echo $1: status $2; exit $2; }

# Check if what to execute is provided as a command-line argument
if [ -z "$1" ]; then
    echo "Usage: $0 <config_name>"
    echo "what: hlt, ofp2, ofr3 or all."
    exit 1
fi

what=$1

hlt_config_names=("l1a140" "l1a200" "ttb200")
ofp2_config_names=("ofp2l1a140" "ofp2l1a200" "ofp2ttb200")
ofr3_config_names=("ofr3ephhlt")

run_hlt() {
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
    cp -p Phase2Timing_resources_abs.json ${logdir}
    cp -p ${cfg} ${logdir}

    #Compress log files
    find ${logdir} -iname stderr -exec gzip {} \;
  done

  echo "All HLT configurations have been processed successfully."
}

run_ofp2() {
  for config_name in "${ofp2_config_names[@]}"; do
    jobs=2
    threads=64
    streams=64
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
      -k Phase2TimingOffline_resources.json --event-skip 100 --event-resolution 25 --wait 30 \
      --logdir ${logdir} \
      --slot n=0,1:nv= \
      --slot n=2,3:nv=  -- ${cfg} | tee ${logdir}/output.log

    mergeResourcesJson.py ${logdir}/step*/pid*/Phase2TimingOffline_resources.json > Phase2TimingOffline_resources.json
    cp -p Phase2TimingOffline_resources.json ${logdir}
    cp -p ${cfg} ${logdir}

    #Compress log files
    find ${logdir} -iname stderr -exec gzip {} \;
  done

  echo "All OFP2 configurations have been processed successfully."
}

run_ofr3() {
    for config_name in "${ofp2_config_names[@]}"; do
	jobs=2
	threads=64
	streams=64
	events=1000
	logdir="logs.Milan.OFR3.$config_name.${jobs}j.${threads}t.${streams}s"
	cfg="Run3_RAW2DIGI_RECO_${config_name}.py"
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
				    -k Run3TimingOffline_resources.json --event-skip 100 --event-resolution 25 --wait 30 \
				    --logdir ${logdir} \
				    --slot n=0,1:nv= \
				    --slot n=2,3:nv=  -- ${cfg} | tee ${logdir}/output.log

	mergeResourcesJson.py ${logdir}/step*/pid*/Run3TimingOffline_resources.json > Run3TimingOffline_resources.json
	cp -p Run3TimingOffline_resources.json ${logdir}
	cp -p ${cfg} ${logdir}

	#Compress log files
	find ${logdir} -iname stderr -exec gzip {} \;
    done

    echo "All OFR3 configurations have been processed successfully."
}

# Run the specified measurements
case ${what} in
  hlt)
    run_hlt
    ;;
  ofp2)
    run_ofp2
    ;;
  ofr3)
    run_ofr3
    ;;  
  all)
    run_hlt
    run_ofp2
    run_ofr3
    ;;
  *)
    echo "Invalide configuration. Please choose a measurement between hlt, ofp2, or all."
    ;;
esac

