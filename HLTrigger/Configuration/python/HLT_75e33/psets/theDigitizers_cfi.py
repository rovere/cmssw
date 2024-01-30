import FWCore.ParameterSet.Config as cms

theDigitizers = cms.PSet(
    MC_fastTimingLayer = cms.PSet(
        HepMCProductLabel = cms.InputTag("generatorSmeared"),
        MaxPseudoRapidity = cms.double(5.0),
        MinEnergy = cms.double(0.5),
        accumulatorType = cms.string('MtdTruthAccumulator'),
        allowDifferentSimHitProcesses = cms.bool(False),
        bunchspace = cms.uint32(25),
        genParticleCollection = cms.InputTag("genParticles"),
        maximumPreviousBunchCrossing = cms.uint32(0),
        maximumSubsequentBunchCrossing = cms.uint32(0),
        premixStage1 = cms.bool(False),
        simHitCollections = cms.PSet(
            mtdCollections = cms.VInputTag(cms.InputTag("g4SimHits","FastTimerHitsBarrel"), cms.InputTag("g4SimHits","FastTimerHitsEndcap"))
        ),
        simTrackCollection = cms.InputTag("g4SimHits"),
        simVertexCollection = cms.InputTag("g4SimHits")
    ),
    calotruth = cms.PSet(
        HepMCProductLabel = cms.InputTag("generatorSmeared"),
        MaxPseudoRapidity = cms.double(5.0),
        MinEnergy = cms.double(0.5),
        accumulatorType = cms.string('CaloTruthAccumulator'),
        allowDifferentSimHitProcesses = cms.bool(False),
        doHGCAL = cms.bool(True),
        genParticleCollection = cms.InputTag("genParticles"),
        maximumPreviousBunchCrossing = cms.uint32(0),
        maximumSubsequentBunchCrossing = cms.uint32(0),
        premixStage1 = cms.bool(False),
        simHitCollections = cms.PSet(
            hgc = cms.VInputTag(cms.InputTag("g4SimHits","HGCHitsEE"), cms.InputTag("g4SimHits","HGCHitsHEfront"), cms.InputTag("g4SimHits","HGCHitsHEback"))
        ),
        simTrackCollection = cms.InputTag("g4SimHits"),
        simVertexCollection = cms.InputTag("g4SimHits")
    ),
    ecal = cms.PSet(
        ConstantTerm = cms.double(0.003),
        EBCorrNoiseMatrixG01 = cms.vdouble(
            1.0, 0.73354, 0.64442, 0.58851, 0.55425,
            0.53082, 0.51916, 0.51097, 0.50732, 0.50409
        ),
        EBCorrNoiseMatrixG06 = cms.vdouble(
            1.0, 0.70946, 0.58021, 0.49846, 0.45006,
            0.41366, 0.39699, 0.38478, 0.37847, 0.37055
        ),
        EBCorrNoiseMatrixG12 = cms.vdouble(
            1.0, 0.71073, 0.55721, 0.46089, 0.40449,
            0.35931, 0.33924, 0.32439, 0.31581, 0.30481
        ),
        EBdigiCollection = cms.string(''),
        EBs25notContainment = cms.double(0.9675),
        EECorrNoiseMatrixG01 = cms.vdouble(
            1.0, 0.72698, 0.62048, 0.55691, 0.51848,
            0.49147, 0.47813, 0.47007, 0.46621, 0.46265
        ),
        EECorrNoiseMatrixG06 = cms.vdouble(
            1.0, 0.71217, 0.47464, 0.34056, 0.26282,
            0.20287, 0.17734, 0.16256, 0.15618, 0.14443
        ),
        EECorrNoiseMatrixG12 = cms.vdouble(
            1.0, 0.71373, 0.44825, 0.30152, 0.21609,
            0.14786, 0.11772, 0.10165, 0.09465, 0.08098
        ),
        EEdigiCollection = cms.string(''),
        EEs25notContainment = cms.double(0.968),
        ESdigiCollection = cms.string(''),
        EcalPreMixStage1 = cms.bool(False),
        EcalPreMixStage2 = cms.bool(False),
        UseLCcorrection = cms.untracked.bool(True),
        accumulatorType = cms.string('EcalDigiProducer'),
        apdAddToBarrel = cms.bool(False),
        apdDigiTag = cms.string('APD'),
        apdDoPEStats = cms.bool(True),
        apdNonlParms = cms.vdouble(
            1.48, -3.75, 1.81, 1.26, 2.0,
            45, 1.0
        ),
        apdSeparateDigi = cms.bool(True),
        apdShapeTau = cms.double(40.5),
        apdShapeTstart = cms.double(74.5),
        apdSimToPEHigh = cms.double(88200000.0),
        apdSimToPELow = cms.double(2450000.0),
        apdTimeOffWidth = cms.double(0.8),
        apdTimeOffset = cms.double(-13.5),
        applyConstantTerm = cms.bool(True),
        binOfMaximum = cms.int32(6),
        componentAddToBarrel = cms.bool(False),
        componentDigiTag = cms.string('Component'),
        componentSeparateDigi = cms.bool(False),
        componentTimePhase = cms.double(0.0),
        componentTimeTag = cms.string('Component'),
        cosmicsPhase = cms.bool(False),
        cosmicsShift = cms.double(0.0),
        doEB = cms.bool(True),
        doEE = cms.bool(False),
        doENoise = cms.bool(True),
        doES = cms.bool(False),
        doESNoise = cms.bool(True),
        doFast = cms.bool(True),
        doPhotostatistics = cms.bool(True),
        hitsProducer = cms.string('g4SimHits'),
        makeDigiSimLinks = cms.untracked.bool(False),
        photoelectronsToAnalogBarrel = cms.double(0.000444444),
        photoelectronsToAnalogEndcap = cms.double(0.000555555),
        samplingFactor = cms.double(1.0),
        simHitToPhotoelectronsBarrel = cms.double(2250.0),
        simHitToPhotoelectronsEndcap = cms.double(1800.0),
        syncPhase = cms.bool(True),
        timeDependent = cms.bool(False),
        timePhase = cms.double(0.0)
    ),
    ecalTime = cms.PSet(
        EBtimeDigiCollection = cms.string('EBTimeDigi'),
        EEtimeDigiCollection = cms.string('EETimeDigi'),
        accumulatorType = cms.string('EcalTimeDigiProducer'),
        componentWaveform = cms.bool(False),
        hitsProducerEB = cms.InputTag("g4SimHits","EcalHitsEB"),
        hitsProducerEE = cms.InputTag("g4SimHits","EcalHitsEE"),
        timeLayerBarrel = cms.int32(7),
        timeLayerEndcap = cms.int32(3)
    ),
    fastTimingLayer = cms.PSet(
        accumulatorType = cms.string('MTDDigiProducer'),
        barrelDigitizer = cms.PSet(
            DeviceSimulation = cms.PSet(
                LCEpositionSlope = cms.double(0.071),
                LightCollectionEff = cms.double(0.25),
                LightCollectionSlope = cms.double(0.075),
                LightYield = cms.double(40000.0),
                PhotonDetectionEff = cms.double(0.2),
                bxTime = cms.double(25)
            ),
            ElectronicsSimulation = cms.PSet(
                ChannelTimeOffset = cms.double(0.0),
                CorrelationCoefficient = cms.double(1.0),
                DarkCountRate = cms.double(10.0),
                EnergyThreshold = cms.double(4.0),
                LightCollectionEff = cms.double(0.25),
                LightYield = cms.double(40000.0),
                Npe_to_V = cms.double(0.0064),
                Npe_to_pC = cms.double(0.016),
                PhotonDetectionEff = cms.double(0.2),
                ReferencePulseNpe = cms.double(100.0),
                ScintillatorDecayTime = cms.double(40.0),
                ScintillatorRiseTime = cms.double(1.1),
                SigmaClock = cms.double(0.015),
                SigmaElectronicNoise = cms.double(1.0),
                SigmaRelTOFHIRenergy = cms.vdouble(0.139, -4.35e-05, 3.315e-09, -1.2e-13, 1.67e-18),
                SinglePhotonTimeResolution = cms.double(0.06),
                SmearTimeForOOTtails = cms.bool(True),
                TestBeamMIPTimeRes = cms.double(4.293),
                TimeThreshold1 = cms.double(20.0),
                TimeThreshold2 = cms.double(50.0),
                adcNbits = cms.uint32(10),
                adcSaturation_MIP = cms.double(600.0),
                adcThreshold_MIP = cms.double(0.064),
                bxTime = cms.double(25),
                smearChannelTimeOffset = cms.double(0.0),
                tdcNbits = cms.uint32(10),
                toaLSB_ns = cms.double(0.02)
            ),
            digiCollectionTag = cms.string('FTLBarrel'),
            digitizerName = cms.string('BTLDigitizer'),
            inputSimHits = cms.InputTag("g4SimHits","FastTimerHitsBarrel"),
            maxSimHitsAccTime = cms.uint32(100),
            premixStage1 = cms.bool(False),
            premixStage1MaxCharge = cms.double(1000000.0),
            premixStage1MinCharge = cms.double(0.0001)
        ),
        endcapDigitizer = cms.PSet(
            DeviceSimulation = cms.PSet(
                bxTime = cms.double(25),
                meVPerMIP = cms.double(0.015),
                tofDelay = cms.double(1)
            ),
            ElectronicsSimulation = cms.PSet(
                FluenceVsRadius = cms.string('1.937*TMath::Power(x,-1.706)'),
                IntegratedLuminosity = cms.double(1000.0),
                LGADGainVsFluence = cms.string('TMath::Min(15.,30.-x)'),
                TimeResolution2 = cms.string('0.0225/x'),
                adcNbits = cms.uint32(8),
                adcSaturation_MIP = cms.double(25),
                adcThreshold_MIP = cms.double(0.025),
                bxTime = cms.double(25),
                tdcNbits = cms.uint32(11),
                toaLSB_ns = cms.double(0.013)
            ),
            digiCollectionTag = cms.string('FTLEndcap'),
            digitizerName = cms.string('ETLDigitizer'),
            inputSimHits = cms.InputTag("g4SimHits","FastTimerHitsEndcap"),
            maxSimHitsAccTime = cms.uint32(100),
            premixStage1 = cms.bool(False),
            premixStage1MaxCharge = cms.double(1000000.0),
            premixStage1MinCharge = cms.double(0.0001)
        ),
        makeDigiSimLinks = cms.bool(False),
        verbosity = cms.untracked.uint32(0)
    ),
    hcal = cms.PSet(
        DelivLuminosity = cms.double(0),
        HBDarkening = cms.bool(False),
        HEDarkening = cms.bool(False),
        HFDarkening = cms.bool(False),
        HFRecalParameterBlock = cms.PSet(
            HFdepthOneParameterA = cms.vdouble(
                0.004123, 0.00602, 0.008201, 0.010489, 0.013379,
                0.016997, 0.021464, 0.027371, 0.034195, 0.044807,
                0.058939, 0.125497
            ),
            HFdepthOneParameterB = cms.vdouble(
                -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05,
                2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107,
                0.000425, 0.000209
            ),
            HFdepthTwoParameterA = cms.vdouble(
                0.002861, 0.004168, 0.0064, 0.008388, 0.011601,
                0.014425, 0.018633, 0.023232, 0.028274, 0.035447,
                0.051579, 0.086593
            ),
            HFdepthTwoParameterB = cms.vdouble(
                -2e-06, -0.0, -7e-06, -6e-06, -2e-06,
                1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05,
                0.000157, -3e-06
            )
        ),
        HcalPreMixStage1 = cms.bool(False),
        HcalPreMixStage2 = cms.bool(False),
        TestNumbering = cms.bool(True),
        accumulatorType = cms.string('HcalDigiProducer'),
        debugCaloSamples = cms.bool(False),
        doEmpty = cms.bool(True),
        doHFWindow = cms.bool(False),
        doIonFeedback = cms.bool(True),
        doNeutralDensityFilter = cms.bool(True),
        doNoise = cms.bool(True),
        doThermalNoise = cms.bool(True),
        doTimeSlew = cms.bool(True),
        hb = cms.PSet(
            binOfMaximum = cms.int32(6),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(True),
            firstRing = cms.int32(1),
            readoutFrameSize = cms.int32(10),
            samplingFactors = cms.vdouble(
                125.44, 125.54, 125.32, 125.13, 124.46,
                125.01, 125.22, 125.48, 124.45, 125.9,
                125.83, 127.01, 126.82, 129.73, 131.83,
                143.52
            ),
            simHitToPhotoelectrons = cms.double(2000.0),
            sipmTau = cms.double(10.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(18.7),
            timePhase = cms.double(6.0),
            timeSmearing = cms.bool(True)
        ),
        he = cms.PSet(
            binOfMaximum = cms.int32(6),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(True),
            firstRing = cms.int32(16),
            readoutFrameSize = cms.int32(10),
            samplingFactors = cms.vdouble(
                210.55, 197.93, 186.12, 189.64, 189.63,
                189.96, 190.03, 190.11, 190.18, 190.25,
                190.32, 190.4, 190.47, 190.54, 190.61,
                190.69, 190.83, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94,
                190.94, 190.94, 190.94, 190.94, 190.94
            ),
            simHitToPhotoelectrons = cms.double(2000.0),
            sipmTau = cms.double(10.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(18.7),
            timePhase = cms.double(6.0),
            timeSmearing = cms.bool(True)
        ),
        hf1 = cms.PSet(
            binOfMaximum = cms.int32(2),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            photoelectronsToAnalog = cms.double(2.79),
            readoutFrameSize = cms.int32(3),
            samplingFactor = cms.double(0.37),
            simHitToPhotoelectrons = cms.double(6.0),
            sipmTau = cms.double(0.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(3.0),
            timePhase = cms.double(9.0)
        ),
        hf2 = cms.PSet(
            binOfMaximum = cms.int32(2),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            photoelectronsToAnalog = cms.double(1.843),
            readoutFrameSize = cms.int32(3),
            samplingFactor = cms.double(0.37),
            simHitToPhotoelectrons = cms.double(6.0),
            sipmTau = cms.double(0.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(3.0),
            timePhase = cms.double(8.0)
        ),
        hitsProducer = cms.string('g4SimHits'),
        ho = cms.PSet(
            binOfMaximum = cms.int32(5),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            firstRing = cms.int32(1),
            readoutFrameSize = cms.int32(10),
            samplingFactors = cms.vdouble(
                231.0, 231.0, 231.0, 231.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0
            ),
            siPMCode = cms.int32(1),
            simHitToPhotoelectrons = cms.double(4000.0),
            sipmTau = cms.double(5.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(-999.0),
            timePhase = cms.double(5.0),
            timeSmearing = cms.bool(False)
        ),
        hoHamamatsu = cms.PSet(
            binOfMaximum = cms.int32(5),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            firstRing = cms.int32(1),
            readoutFrameSize = cms.int32(10),
            samplingFactors = cms.vdouble(
                231.0, 231.0, 231.0, 231.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0
            ),
            siPMCode = cms.int32(2),
            simHitToPhotoelectrons = cms.double(4000.0),
            sipmTau = cms.double(5.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(-999.0),
            timePhase = cms.double(5.0),
            timeSmearing = cms.bool(False)
        ),
        hoZecotek = cms.PSet(
            binOfMaximum = cms.int32(5),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            firstRing = cms.int32(1),
            readoutFrameSize = cms.int32(10),
            samplingFactors = cms.vdouble(
                231.0, 231.0, 231.0, 231.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0,
                360.0, 360.0, 360.0, 360.0, 360.0
            ),
            siPMCode = cms.int32(2),
            simHitToPhotoelectrons = cms.double(4000.0),
            sipmTau = cms.double(5.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(-999.0),
            timePhase = cms.double(5.0),
            timeSmearing = cms.bool(False)
        ),
        ignoreGeantTime = cms.bool(False),
        injectTestHits = cms.bool(False),
        injectTestHitsCells = cms.vint32(),
        injectTestHitsEnergy = cms.vdouble(),
        injectTestHitsTime = cms.vdouble(),
        killHE = cms.bool(True),
        makeDigiSimLinks = cms.untracked.bool(False),
        minFCToDelay = cms.double(5.0),
        zdc = cms.PSet(
            binOfMaximum = cms.int32(5),
            delayQIE = cms.int32(-999),
            doPhotoStatistics = cms.bool(True),
            doSiPMSmearing = cms.bool(False),
            photoelectronsToAnalog = cms.double(1.843),
            readoutFrameSize = cms.int32(10),
            samplingFactor = cms.double(1.0),
            simHitToPhotoelectrons = cms.double(6.0),
            sipmTau = cms.double(0.0),
            syncPhase = cms.bool(True),
            threshold_currentTDC = cms.double(-999.0),
            timePhase = cms.double(-4.0)
        )
    ),
    hgceeDigitizer = cms.PSet(
        NoiseGeneration_Method = cms.bool(True),
        accumulatorType = cms.string('HGCDigiProducer'),
        bxTime = cms.double(25),
        digiCfg = cms.PSet(
            cceParams = cms.PSet(
                refToPSet_ = cms.string('HGCAL_cceParams_toUse')
            ),
            chargeCollectionEfficiencies = cms.PSet(
                refToPSet_ = cms.string('HGCAL_chargeCollectionEfficiencies')
            ),
            doTimeSamples = cms.bool(False),
            feCfg = cms.PSet(
                adcNbits = cms.uint32(10),
                adcPulse = cms.vdouble(
                    0.0, 0.017, 0.817, 0.163, 0.003,
                    0.0
                ),
                adcSaturation_fC = cms.double(100),
                adcThreshold_fC = cms.double(0.672),
                fwVersion = cms.uint32(2),
                jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02),
                jitterNoise_ns = cms.vdouble(5.0, 5.0, 5.0),
                pulseAvgT = cms.vdouble(
                    0.0, 23.42298, 13.16733, 6.41062, 5.03946,
                    4.532
                ),
                targetMIPvalue_ADC = cms.uint32(10),
                tdcChargeDrainParameterisation = cms.vdouble(
                    -919.13, 365.36, -14.1, 0.2, -21.85,
                    49.39, 22.21, 0.8, -0.28, 27.14,
                    43.95, 3.89048
                ),
                tdcForToAOnset_fC = cms.vdouble(12.0, 12.0, 12.0),
                tdcNbits = cms.uint32(12),
                tdcOnset_fC = cms.double(60),
                tdcResolutionInPs = cms.double(0.001),
                tdcSaturation_fC = cms.double(10000),
                toaLSB_ns = cms.double(0.0244),
                toaMode = cms.uint32(1)
            ),
            ileakParam = cms.PSet(
                refToPSet_ = cms.string('HGCAL_ileakParam_toUse')
            ),
            keV2fC = cms.double(0.044259),
            noise_fC = cms.PSet(
                refToPSet_ = cms.string('HGCAL_noise_fC')
            ),
            thresholdFollowsMIP = cms.bool(True)
        ),
        digiCollection = cms.string('HGCDigisEE'),
        digitizationType = cms.uint32(0),
        digitizer = cms.string('HGCEEDigitizer'),
        eVPerEleHolePair = cms.double(3.62),
        hitCollection = cms.string('HGCHitsEE'),
        hitsProducer = cms.string('g4SimHits'),
        makeDigiSimLinks = cms.bool(False),
        maxSimHitsAccTime = cms.uint32(100),
        premixStage1 = cms.bool(False),
        premixStage1MaxCharge = cms.double(1000000.0),
        premixStage1MinCharge = cms.double(0),
        tofDelay = cms.double(-9),
        useAllChannels = cms.bool(True),
        verbosity = cms.untracked.uint32(0)
    ),
    hgchebackDigitizer = cms.PSet(
        NoiseGeneration_Method = cms.bool(True),
        accumulatorType = cms.string('HGCDigiProducer'),
        bxTime = cms.double(25),
        digiCfg = cms.PSet(
            algo = cms.uint32(2),
            doTimeSamples = cms.bool(False),
            feCfg = cms.PSet(
                adcNbits = cms.uint32(10),
                adcPulse = cms.vdouble(
                    0.0, 0.017, 0.817, 0.163, 0.003,
                    0.0
                ),
                adcSaturation_fC = cms.double(68.75),
                adcThreshold_fC = cms.double(0.5),
                fwVersion = cms.uint32(2),
                jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02),
                jitterNoise_ns = cms.vdouble(5.0, 5.0, 5.0),
                pulseAvgT = cms.vdouble(
                    0.0, 23.42298, 13.16733, 6.41062, 5.03946,
                    4.532
                ),
                targetMIPvalue_ADC = cms.uint32(15),
                tdcChargeDrainParameterisation = cms.vdouble(
                    -919.13, 365.36, -14.1, 0.2, -21.85,
                    49.39, 22.21, 0.8, -0.28, 27.14,
                    43.95, 3.89048
                ),
                tdcForToAOnset_fC = cms.vdouble(12.0, 12.0, 12.0),
                tdcNbits = cms.uint32(12),
                tdcOnset_fC = cms.double(55),
                tdcResolutionInPs = cms.double(0.001),
                tdcSaturation_fC = cms.double(1000),
                toaLSB_ns = cms.double(0.0244),
                toaMode = cms.uint32(1)
            ),
            keV2MIP = cms.double(0.0014814814814814814),
            nPEperMIP = cms.double(21.0),
            nTotalPE = cms.double(7500),
            noise = cms.PSet(
                refToPSet_ = cms.string('HGCAL_noise_heback')
            ),
            sdPixels = cms.double(1e-06),
            thresholdFollowsMIP = cms.bool(True)
        ),
        digiCollection = cms.string('HGCDigisHEback'),
        digitizationType = cms.uint32(1),
        digitizer = cms.string('HGCHEbackDigitizer'),
        hitCollection = cms.string('HGCHitsHEback'),
        hitsProducer = cms.string('g4SimHits'),
        makeDigiSimLinks = cms.bool(False),
        maxSimHitsAccTime = cms.uint32(100),
        premixStage1 = cms.bool(False),
        premixStage1MaxCharge = cms.double(1000000.0),
        premixStage1MinCharge = cms.double(0),
        tofDelay = cms.double(-14),
        useAllChannels = cms.bool(True),
        verbosity = cms.untracked.uint32(0)
    ),
    hgchefrontDigitizer = cms.PSet(
        NoiseGeneration_Method = cms.bool(True),
        accumulatorType = cms.string('HGCDigiProducer'),
        bxTime = cms.double(25),
        digiCfg = cms.PSet(
            cceParams = cms.PSet(
                refToPSet_ = cms.string('HGCAL_cceParams_toUse')
            ),
            chargeCollectionEfficiencies = cms.PSet(
                refToPSet_ = cms.string('HGCAL_chargeCollectionEfficiencies')
            ),
            doTimeSamples = cms.bool(False),
            feCfg = cms.PSet(
                adcNbits = cms.uint32(10),
                adcPulse = cms.vdouble(
                    0.0, 0.017, 0.817, 0.163, 0.003,
                    0.0
                ),
                adcSaturation_fC = cms.double(100),
                adcThreshold_fC = cms.double(0.672),
                fwVersion = cms.uint32(2),
                jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02),
                jitterNoise_ns = cms.vdouble(5.0, 5.0, 5.0),
                pulseAvgT = cms.vdouble(
                    0.0, 23.42298, 13.16733, 6.41062, 5.03946,
                    4.532
                ),
                targetMIPvalue_ADC = cms.uint32(10),
                tdcChargeDrainParameterisation = cms.vdouble(
                    -919.13, 365.36, -14.1, 0.2, -21.85,
                    49.39, 22.21, 0.8, -0.28, 27.14,
                    43.95, 3.89048
                ),
                tdcForToAOnset_fC = cms.vdouble(12.0, 12.0, 12.0),
                tdcNbits = cms.uint32(12),
                tdcOnset_fC = cms.double(60),
                tdcResolutionInPs = cms.double(0.001),
                tdcSaturation_fC = cms.double(10000),
                toaLSB_ns = cms.double(0.0244),
                toaMode = cms.uint32(1)
            ),
            ileakParam = cms.PSet(
                refToPSet_ = cms.string('HGCAL_ileakParam_toUse')
            ),
            keV2fC = cms.double(0.044259),
            noise_fC = cms.PSet(
                refToPSet_ = cms.string('HGCAL_noise_fC')
            ),
            thresholdFollowsMIP = cms.bool(True)
        ),
        digiCollection = cms.string('HGCDigisHEfront'),
        digitizationType = cms.uint32(0),
        digitizer = cms.string('HGCHEfrontDigitizer'),
        hitCollection = cms.string('HGCHitsHEfront'),
        hitsProducer = cms.string('g4SimHits'),
        makeDigiSimLinks = cms.bool(False),
        maxSimHitsAccTime = cms.uint32(100),
        premixStage1 = cms.bool(False),
        premixStage1MaxCharge = cms.double(1000000.0),
        premixStage1MinCharge = cms.double(0),
        tofDelay = cms.double(-11),
        useAllChannels = cms.bool(True),
        verbosity = cms.untracked.uint32(0)
    ),
    mergedtruth = cms.PSet(
        HepMCProductLabel = cms.InputTag("generatorSmeared"),
        accumulatorType = cms.string('TrackingTruthAccumulator'),
        allowDifferentSimHitProcesses = cms.bool(False),
        alwaysAddAncestors = cms.bool(True),
        createInitialVertexCollection = cms.bool(False),
        createMergedBremsstrahlung = cms.bool(True),
        createUnmergedCollection = cms.bool(True),
        genParticleCollection = cms.InputTag("genParticles"),
        ignoreTracksOutsideVolume = cms.bool(False),
        maximumPreviousBunchCrossing = cms.uint32(9999),
        maximumSubsequentBunchCrossing = cms.uint32(9999),
        removeDeadModules = cms.bool(False),
        select = cms.PSet(
            chargedOnlyTP = cms.bool(True),
            intimeOnlyTP = cms.bool(False),
            lipTP = cms.double(1000),
            maxRapidityTP = cms.double(5.0),
            minHitTP = cms.int32(0),
            minRapidityTP = cms.double(-5.0),
            pdgIdTP = cms.vint32(),
            ptMaxTP = cms.double(1e+100),
            ptMinTP = cms.double(0.1),
            signalOnlyTP = cms.bool(True),
            stableOnlyTP = cms.bool(False),
            tipTP = cms.double(1000)
        ),
        simHitCollections = cms.PSet(
            muon = cms.VInputTag(cms.InputTag("g4SimHits","MuonDTHits"), cms.InputTag("g4SimHits","MuonCSCHits"), cms.InputTag("g4SimHits","MuonRPCHits"), cms.InputTag("g4SimHits","MuonGEMHits")),
            pixel = cms.VInputTag(cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"), cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"), cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"), cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof")),
            tracker = cms.VInputTag()
        ),
        simTrackCollection = cms.InputTag("g4SimHits"),
        simVertexCollection = cms.InputTag("g4SimHits"),
        vertexDistanceCut = cms.double(0.003),
        volumeRadius = cms.double(120.0),
        volumeZ = cms.double(300.0)
    ),
    pixel = cms.PSet(
        AlgorithmCommon = cms.PSet(
            DeltaProductionCut = cms.double(0.03),
            makeDigiSimLinks = cms.untracked.bool(True)
        ),
        GeometryType = cms.string('idealForDigi'),
        PSPDigitizerAlgorithm = cms.PSet(
            AdcFullScale = cms.int32(255),
            AddInefficiency = cms.bool(False),
            AddNoise = cms.bool(True),
            AddNoisyPixels = cms.bool(True),
            AddThresholdSmearing = cms.bool(False),
            AddXTalk = cms.bool(True),
            Alpha2Order = cms.bool(True),
            BiasRailInefficiencyFlag = cms.int32(1),
            CellsToKill = cms.VPSet(),
            ClusterWidth = cms.double(3),
            DeadModules = cms.VPSet(),
            DeadModules_DB = cms.bool(False),
            EfficiencyFactors_Barrel = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999
            ),
            EfficiencyFactors_Endcap = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999
            ),
            ElectronPerAdc = cms.double(135.0),
            HIPThresholdInElectrons_Barrel = cms.double(10000000000.0),
            HIPThresholdInElectrons_Endcap = cms.double(10000000000.0),
            Inefficiency_DB = cms.bool(False),
            InterstripCoupling = cms.double(0.05),
            KillModules = cms.bool(False),
            LorentzAngle_DB = cms.bool(True),
            NoiseInElectrons = cms.double(200),
            Phase2ReadoutMode = cms.int32(0),
            ReadoutNoiseInElec = cms.double(-99.9),
            SigmaCoeff = cms.double(1.8),
            SigmaZero = cms.double(0.00037),
            TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
            TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
            ThresholdInElectrons_Barrel = cms.double(6300.0),
            ThresholdInElectrons_Endcap = cms.double(6300.0),
            ThresholdSmearing_Barrel = cms.double(630.0),
            ThresholdSmearing_Endcap = cms.double(630.0),
            TofLowerCut = cms.double(-12.5),
            TofUpperCut = cms.double(12.5),
            UseReweighting = cms.bool(False)
        ),
        PSSDigitizerAlgorithm = cms.PSet(
            AdcFullScale = cms.int32(255),
            AddInefficiency = cms.bool(False),
            AddNoise = cms.bool(True),
            AddNoisyPixels = cms.bool(True),
            AddThresholdSmearing = cms.bool(False),
            AddXTalk = cms.bool(True),
            Alpha2Order = cms.bool(True),
            CellsToKill = cms.VPSet(),
            ClusterWidth = cms.double(3),
            DeadModules = cms.VPSet(),
            DeadModules_DB = cms.bool(False),
            EfficiencyFactors_Barrel = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999
            ),
            EfficiencyFactors_Endcap = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999
            ),
            ElectronPerAdc = cms.double(135.0),
            HIPThresholdInElectrons_Barrel = cms.double(21000.0),
            HIPThresholdInElectrons_Endcap = cms.double(21000.0),
            Inefficiency_DB = cms.bool(False),
            InterstripCoupling = cms.double(0.05),
            KillModules = cms.bool(False),
            LorentzAngle_DB = cms.bool(True),
            NoiseInElectrons = cms.double(1010),
            Phase2ReadoutMode = cms.int32(0),
            ReadoutNoiseInElec = cms.double(-99.9),
            SigmaCoeff = cms.double(1.8),
            SigmaZero = cms.double(0.00037),
            TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
            TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
            ThresholdInElectrons_Barrel = cms.double(4800.0),
            ThresholdInElectrons_Endcap = cms.double(4800.0),
            ThresholdSmearing_Barrel = cms.double(480.0),
            ThresholdSmearing_Endcap = cms.double(480.0),
            TofLowerCut = cms.double(-12.5),
            TofUpperCut = cms.double(12.5),
            UseReweighting = cms.bool(False)
        ),
        Pixel3DDigitizerAlgorithm = cms.PSet(
            AdcFullScale = cms.int32(15),
            AddInefficiency = cms.bool(False),
            AddNoise = cms.bool(False),
            AddNoisyPixels = cms.bool(False),
            AddThresholdSmearing = cms.bool(False),
            AddXTalk = cms.bool(False),
            Alpha2Order = cms.bool(True),
            ApplyTimewalk = cms.bool(False),
            CellsToKill = cms.VPSet(),
            ClusterWidth = cms.double(3),
            DeadModules = cms.VPSet(),
            DeadModules_DB = cms.bool(False),
            EfficiencyFactors_Barrel = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999
            ),
            EfficiencyFactors_Endcap = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999
            ),
            ElectronPerAdc = cms.double(1500.0),
            Even_column_interchannelCoupling_next_column = cms.double(0.0),
            Even_row_interchannelCoupling_next_row = cms.double(0.0),
            HIPThresholdInElectrons_Barrel = cms.double(10000000000.0),
            HIPThresholdInElectrons_Endcap = cms.double(10000000000.0),
            Inefficiency_DB = cms.bool(False),
            InterstripCoupling = cms.double(0.0),
            KillModules = cms.bool(False),
            LorentzAngle_DB = cms.bool(True),
            NPColumnGap = cms.double(46.0),
            NPColumnRadius = cms.double(4.0),
            NoiseInElectrons = cms.double(0.0),
            Odd_column_interchannelCoupling_next_column = cms.double(0.0),
            Odd_row_interchannelCoupling_next_row = cms.double(0.2),
            OhmicColumnRadius = cms.double(4.0),
            Phase2ReadoutMode = cms.int32(3),
            ReadoutNoiseInElec = cms.double(-99.9),
            SigmaCoeff = cms.double(1.8),
            SigmaZero = cms.double(0.00037),
            TanLorentzAnglePerTesla_Barrel = cms.double(0.106),
            TanLorentzAnglePerTesla_Endcap = cms.double(0.106),
            ThresholdInElectrons_Barrel = cms.double(1000.0),
            ThresholdInElectrons_Endcap = cms.double(1000.0),
            ThresholdSmearing_Barrel = cms.double(0.0),
            ThresholdSmearing_Endcap = cms.double(0.0),
            TimewalkModel = cms.PSet(
                Curves = cms.VPSet(
                    cms.PSet(
                        charge = cms.vdouble(
                            1000, 1025, 1050, 1100, 1200,
                            1500, 2000, 6000, 10000, 15000,
                            20000, 30000
                        ),
                        delay = cms.vdouble(
                            26.8, 23.73, 21.92, 19.46, 16.52,
                            12.15, 8.88, 3.03, 1.69, 0.95,
                            0.56, 0.19
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            1200, 1225, 1250, 1500, 2000,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            26.28, 23.5, 21.79, 14.92, 10.27,
                            3.33, 1.86, 1.07, 0.66, 0.27
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            1500, 1525, 1550, 1600, 2000,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            25.36, 23.05, 21.6, 19.56, 12.94,
                            3.79, 2.14, 1.26, 0.81, 0.39
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            3000, 3025, 3050, 3100, 3500,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            25.63, 23.63, 22.35, 20.65, 14.92,
                            6.7, 3.68, 2.29, 1.62, 1.02
                        )
                    )
                ),
                ThresholdValues = cms.vdouble(1000, 1200, 1500, 3000)
            ),
            TofLowerCut = cms.double(-5.0),
            TofUpperCut = cms.double(20.0),
            UseReweighting = cms.bool(False)
        ),
        PixelDigitizerAlgorithm = cms.PSet(
            AdcFullScale = cms.int32(15),
            AddInefficiency = cms.bool(False),
            AddNoise = cms.bool(False),
            AddNoisyPixels = cms.bool(False),
            AddThresholdSmearing = cms.bool(False),
            AddXTalk = cms.bool(False),
            Alpha2Order = cms.bool(True),
            ApplyTimewalk = cms.bool(False),
            CellsToKill = cms.VPSet(),
            ClusterWidth = cms.double(3),
            DeadModules = cms.VPSet(),
            DeadModules_DB = cms.bool(False),
            EfficiencyFactors_Barrel = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999
            ),
            EfficiencyFactors_Endcap = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999
            ),
            ElectronPerAdc = cms.double(1500.0),
            Even_column_interchannelCoupling_next_column = cms.double(0.0),
            Even_row_interchannelCoupling_next_row = cms.double(0.0),
            HIPThresholdInElectrons_Barrel = cms.double(10000000000.0),
            HIPThresholdInElectrons_Endcap = cms.double(10000000000.0),
            Inefficiency_DB = cms.bool(False),
            InterstripCoupling = cms.double(0.0),
            KillModules = cms.bool(False),
            LorentzAngle_DB = cms.bool(True),
            NoiseInElectrons = cms.double(0.0),
            Odd_column_interchannelCoupling_next_column = cms.double(0.0),
            Odd_row_interchannelCoupling_next_row = cms.double(0.2),
            Phase2ReadoutMode = cms.int32(3),
            ReadoutNoiseInElec = cms.double(-99.9),
            SigmaCoeff = cms.double(0),
            SigmaZero = cms.double(0.00037),
            TanLorentzAnglePerTesla_Barrel = cms.double(0.106),
            TanLorentzAnglePerTesla_Endcap = cms.double(0.106),
            ThresholdInElectrons_Barrel = cms.double(1000.0),
            ThresholdInElectrons_Endcap = cms.double(1000.0),
            ThresholdSmearing_Barrel = cms.double(0.0),
            ThresholdSmearing_Endcap = cms.double(0.0),
            TimewalkModel = cms.PSet(
                Curves = cms.VPSet(
                    cms.PSet(
                        charge = cms.vdouble(
                            1000, 1025, 1050, 1100, 1200,
                            1500, 2000, 6000, 10000, 15000,
                            20000, 30000
                        ),
                        delay = cms.vdouble(
                            26.8, 23.73, 21.92, 19.46, 16.52,
                            12.15, 8.88, 3.03, 1.69, 0.95,
                            0.56, 0.19
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            1200, 1225, 1250, 1500, 2000,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            26.28, 23.5, 21.79, 14.92, 10.27,
                            3.33, 1.86, 1.07, 0.66, 0.27
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            1500, 1525, 1550, 1600, 2000,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            25.36, 23.05, 21.6, 19.56, 12.94,
                            3.79, 2.14, 1.26, 0.81, 0.39
                        )
                    ),
                    cms.PSet(
                        charge = cms.vdouble(
                            3000, 3025, 3050, 3100, 3500,
                            6000, 10000, 15000, 20000, 30000
                        ),
                        delay = cms.vdouble(
                            25.63, 23.63, 22.35, 20.65, 14.92,
                            6.7, 3.68, 2.29, 1.62, 1.02
                        )
                    )
                ),
                ThresholdValues = cms.vdouble(1000, 1200, 1500, 3000)
            ),
            TofLowerCut = cms.double(-5.0),
            TofUpperCut = cms.double(20.0),
            UseReweighting = cms.bool(False)
        ),
        ROUList = cms.vstring(
            'TrackerHitsPixelBarrelLowTof',
            'TrackerHitsPixelBarrelHighTof',
            'TrackerHitsPixelEndcapLowTof',
            'TrackerHitsPixelEndcapHighTof'
        ),
        SSDigitizerAlgorithm = cms.PSet(
            AdcFullScale = cms.int32(255),
            AddInefficiency = cms.bool(False),
            AddNoise = cms.bool(True),
            AddNoisyPixels = cms.bool(True),
            AddThresholdSmearing = cms.bool(False),
            AddXTalk = cms.bool(True),
            Alpha2Order = cms.bool(True),
            CBCDeadTime = cms.double(0.0),
            CellsToKill = cms.VPSet(),
            ClusterWidth = cms.double(3),
            DeadModules = cms.VPSet(),
            DeadModules_DB = cms.bool(False),
            EfficiencyFactors_Barrel = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999
            ),
            EfficiencyFactors_Endcap = cms.vdouble(
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999, 0.999, 0.999, 0.999, 0.999,
                0.999
            ),
            ElectronPerAdc = cms.double(135.0),
            HIPThresholdInElectrons_Barrel = cms.double(10000000000.0),
            HIPThresholdInElectrons_Endcap = cms.double(10000000000.0),
            HitDetectionMode = cms.int32(0),
            Inefficiency_DB = cms.bool(False),
            InterstripCoupling = cms.double(0.05),
            KillModules = cms.bool(False),
            LorentzAngle_DB = cms.bool(True),
            NoiseInElectrons = cms.double(1263),
            Phase2ReadoutMode = cms.int32(0),
            PulseShapeParameters = cms.vdouble(
                -3.0, 16.043703, 99.999857, 40.57165, 2.0,
                1.2459094
            ),
            ReadoutNoiseInElec = cms.double(-99.9),
            SigmaCoeff = cms.double(1.8),
            SigmaZero = cms.double(0.00037),
            TanLorentzAnglePerTesla_Barrel = cms.double(0.07),
            TanLorentzAnglePerTesla_Endcap = cms.double(0.07),
            ThresholdInElectrons_Barrel = cms.double(6000.0),
            ThresholdInElectrons_Endcap = cms.double(6000.0),
            ThresholdSmearing_Barrel = cms.double(600.0),
            ThresholdSmearing_Endcap = cms.double(600.0),
            TofLowerCut = cms.double(-12.5),
            TofUpperCut = cms.double(12.5),
            UseReweighting = cms.bool(False)
        ),
        accumulatorType = cms.string('Phase2TrackerDigitizer'),
        hitsProducer = cms.string('g4SimHits'),
        isOTreadoutAnalog = cms.bool(False),
        premixStage1 = cms.bool(False),
        usePseudoPixel3DAlgo = cms.bool(False)
    ),
    puVtx = cms.PSet(
        accumulatorType = cms.string('PileupVertexAccumulator'),
        hitsProducer = cms.string('generator'),
        makeDigiSimLinks = cms.untracked.bool(False),
        saveVtxTimes = cms.bool(True),
        vtxFallbackTag = cms.InputTag("generator"),
        vtxTag = cms.InputTag("generatorSmeared")
    )
)