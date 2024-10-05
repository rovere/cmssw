#!/usr/bin/env python3
import os
import argparse
import glob
from colorama import Fore, Style
from rich.console import Console

from DataFormats.FWLite import Events, Handle


### Event Analysis
def analyse_event(event, paths, verbosity=0):
    if verbosity >= 0:
        print("-" * 50)
        print("Run             =", event.eventAuxiliary().run())
        print("LuminosityBlock =", event.eventAuxiliary().luminosityBlock())
        print("Event           =", event.eventAuxiliary().event())
        print("-" * 50)

    TriggerResultsInputTagLabel = "TriggerResults"
    TriggerResultsHandle = Handle("edm::TriggerResults")
    event.getByLabel(TriggerResultsInputTagLabel, TriggerResultsHandle)
    TriggerResults = TriggerResultsHandle.product()
    triggerNames = event.object().triggerNames(TriggerResults).triggerNames()
    if verbosity >= 0:
        print(
            '\nTriggerResults = "'
            + TriggerResultsInputTagLabel
            + '"'
            + f" [{len(triggerNames)} Paths]"
        )
        print("")

    for triggerNameIdx, triggerName in enumerate(triggerNames):
        triggerNameIdxStr = "[" + str(triggerNameIdx) + "]"
        if verbosity >= 0:
            print(
                f"{triggerNameIdxStr:>6} | {int(TriggerResults.accept(triggerNameIdx))} | {str(triggerName)}"
            )
        if triggerName not in paths.keys():
            paths[str(triggerName)] = 0
        paths[str(triggerName)] += int(TriggerResults.accept(triggerNameIdx))

    ###MR    hltTriggerSummaryInputTagLabel = 'hltTriggerSummaryAOD'
    ###MR    hltTriggerSummaryHandle = Handle('trigger::TriggerEvent')
    ###MR    event.getByLabel(hltTriggerSummaryInputTagLabel, hltTriggerSummaryHandle)
    ###MR    hltTriggerSummary = hltTriggerSummaryHandle.product()
    ###MR
    ###MR    trigProcessName = hltTriggerSummary.usedProcessName()
    ###MR    trigSummarySizeFilters = hltTriggerSummary.sizeFilters()
    ###MR    trigObjects = hltTriggerSummary.getObjects()

    ###MR if verbosity >= 1:
    ###MR     print("\n" + "-" * 20)
    ###MR     print(
    ###MR         '\nTriggerEvent = "'
    ###MR         + hltTriggerSummaryInputTagLabel
    ###MR         + "::"
    ###MR         + trigProcessName
    ###MR         + '"'
    ###MR         + f" [{trigSummarySizeFilters} Filters, {len(trigObjects)} TriggerObjects]"
    ###MR     )

    ###MR     if verbosity >= 2:
    ###MR         print("")
    ###MR         for trigFilterIdx in range(trigSummarySizeFilters):
    ###MR             trigFilterTag = hltTriggerSummary.filterTag(trigFilterIdx)
    ###MR             trigFilterKeys = hltTriggerSummary.filterKeys(trigFilterIdx)
    ###MR             trigFilterIds = hltTriggerSummary.filterIds(trigFilterIdx)
    ###MR             trigFilterIdxStr = "[" + str(trigFilterIdx) + "]"

    ###MR             if verbosity < 3:
    ###MR                 print(
    ###MR                     f"{trigFilterIdxStr:>7} {trigFilterTag.encode()} ({len(trigFilterKeys)} TriggerObjects)"
    ###MR                 )
    ###MR             else:
    ###MR                 print("")
    ###MR                 print(
    ###MR                     f"{trigFilterIdxStr:>7} {trigFilterTag.encode()} ({len(trigFilterKeys)} TriggerObjects)"
    ###MR                 )
    ###MR                 for trigFilterKeyIdx, trigFilterKey in enumerate(trigFilterKeys):
    ###MR                     trigObj, trigObjId = (
    ###MR                         trigObjects[trigFilterKey],
    ###MR                         trigFilterIds[trigFilterKeyIdx],
    ###MR                     )
    ###MR                     trigObjStr = f"pt = {trigObj.pt(): >8.3f}, eta = {trigObj.eta(): >-6.3f}, phi = {trigObj.phi(): >-6.3f}, id = {trigObj.id(): >-3d}"
    ###MR                     trigFilterKeyStr = "[" + str(trigFilterKey) + "]"
    ###MR                     print(
    ###MR                         f"        {trigFilterKeyStr:>4} FilterId = {trigObjId: >-3d} : "
    ###MR                         + trigObjStr
    ###MR                     )


### main
if __name__ == "__main__":
    ### args
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--inputs",
        dest="inputs",
        required=True,
        nargs="+",
        default=None,
        help="path to input .root file(s)",
    )

    parser.add_argument(
        "-e",
        "--every",
        dest="every",
        action="store",
        type=int,
        default=1e2,
        help="show progress of processing every N events",
    )

    parser.add_argument(
        "-s",
        "--skipEvents",
        dest="skipEvents",
        action="store",
        type=int,
        default=0,
        help="index of first event to be processed (inclusive)",
    )

    parser.add_argument(
        "-n",
        "--maxEvents",
        dest="maxEvents",
        action="store",
        type=int,
        default=-1,
        help="maximum number of events to be processed (inclusive)",
    )

    parser.add_argument(
        "-v",
        "--verbosity",
        dest="verbosity",
        action="store",
        type=int,
        default=0,
        help="level of verbosity",
    )

    opts, opts_unknown = parser.parse_known_args()

    log_prx = os.path.basename(__file__) + " -- "

    paths_filtered_upon = [
        ###MR    "TripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9",
        ###MR    "TripleTkMuon_5_3p5_2p5_OS_Mass5to17",
        "pDoubleEGEle37_24",  # TAKE
        "pDoubleIsoTkPho22_12",  # TAKE
        "pDoublePuppiJet112_112",  # TAKE
        ###MR    "pDoublePuppiJet160_35_mass620",
        "pDoublePuppiTau52_52",  # TAKE
        "pDoubleTkEle25_12",  # TAKE
        ###MR    "pDoubleTkElePuppiHT_8_8_390",
        ###MR    "pDoubleTkMuPuppiHT_3_3_300",
        ###MR    "pDoubleTkMuPuppiJetPuppiMet_3_3_60_130",
        "pDoubleTkMuon15_7",  # TAKE
        ###MR    "pDoubleTkMuonTkEle5_5_9",
        ###MR    "pDoubleTkMuon_4_4_OS_Dr1p2",
        ###MR    "pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18",
        ###MR    "pDoubleTkMuon_OS_Er1p5_Dr1p4",
        "pIsoTkEleEGEle22_12",  # TAKE
        ###MR    "pNNPuppiTauPuppiMet_55_190",
        "pPuppiHT400",  # TAKE
        "pPuppiHT450",  # TAKE
        "pPuppiMET200",  # TAKE
        ###MR    "pPuppiMHT140",
        ###MR    "pPuppiTauTkIsoEle45_22",
        ###MR    "pPuppiTauTkMuon42_18",
        "pQuadJet70_55_40_40",  # TAKE
        "pSingleEGEle51",  # TAKE
        "pSingleIsoTkEle28",  # TAKE
        "pSingleIsoTkPho36",  # TAKE
        "pSinglePuppiJet230",  # TAKE
        "pSingleTkEle36",  # TAKE
        "pSingleTkMuon22",  # TAKE
        ###MR    "pTkEleIsoPuppiHT_26_190",
        ###MR    "pTkElePuppiJet_28_40_MinDR",
        ###MR    "pTkEleTkMuon10_20",
        ###MR    "pTkMuPuppiJetPuppiMet_3_110_120",
        ###MR    "pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax",
        ###MR    "pTkMuonDoubleTkEle6_17_17",
        ###MR    "pTkMuonPuppiHT6_320",
        ###MR    "pTkMuonTkEle7_23",
        ###MR    "pTkMuonTkIsoEle7_20",
        "pTripleTkMuon5_3_3",  # TAKE
    ]

    ### args validation
    INPUT_FILES = []
    for i_inpf in opts.inputs:
        i_inpf_ls = glob.glob(i_inpf)
        if len(i_inpf_ls) == 0:
            i_inpf_ls = [i_inpf]
        for i_inpf_2 in i_inpf_ls:
            i_inpfile = (
                os.path.abspath(os.path.realpath(i_inpf_2))
                if os.path.isfile(i_inpf_2)
                else i_inpf_2
            )
            INPUT_FILES += [i_inpfile]

    INPUT_FILES = sorted(list(set(INPUT_FILES)))

    if len(INPUT_FILES) == 0:
        raise RuntimeError(log_prx + "empty list of input files [-i]")

    SHOW_EVERY = opts.every
    if SHOW_EVERY <= 0:
        print(
            log_prx
            + 'invalid (non-positive) value for option "-e/--every" ('
            + str(SHOW_EVERY)
            + "), value will be changed to 100"
        )
        SHOW_EVERY = 1e2

    if len(opts_unknown) > 0:
        raise RuntimeError(
            log_prx + "unrecognized command-line arguments: " + str(opts_unknown)
        )

    ## histograms
    nEvtProcessed = 0

    paths = dict()
    for i_inpf in INPUT_FILES:
        if opts.verbosity >= 10:
            print(Fore.RED + os.path.relpath(i_inpf) + Fore.RESET)

        try:
            events = Events(i_inpf)
        except:
            print(
                log_prx
                + 'target TFile does not contain a TTree named "Events" (file will be ignored) [-t]: '
                + i_inpf
            )
            continue

        skipEvents = 0 if opts.skipEvents < 0 else opts.skipEvents

        eventIndex = 0
        for event in events:
            if (eventIndex < skipEvents) or (
                (opts.maxEvents >= 0) and (nEvtProcessed >= opts.maxEvents)
            ):
                continue

            if (not (eventIndex % SHOW_EVERY)) and (eventIndex > 0):
                print("-" * 10)
                print("Event #" + str(eventIndex))

            analyse_event(event, paths, opts.verbosity)
            nEvtProcessed += 1
            eventIndex += 1

    print("=" * 30)
    # Create a console that forces terminal color codes
    console = Console(force_terminal=True)
    console.print("[blue]Events processed =", nEvtProcessed)
    console.print("[red]Trigger Results Job Summary:")
    for path, success in paths.items():
        active = True if path in paths_filtered_upon else False
        bright = Style.BRIGHT if active else ""
        reset = Style.RESET_ALL
        #        path_color = (bright + Fore.BLUE) if active else Fore.BLUE
        #        success_color = (bright + Fore.GREEN) if active else Fore.RED
        #        processed_events_color = (bright + Fore.YELLOW) if active else Fore.RED
        #        path_color = (bright + Fore.BLUE) if active else Fore.BLUE
        path_color = "[bold blue]" if active else "[blue]"
        success_color = "[bold green]" if active else "[green]"
        processed_events_color = "[bold yellow]" if active else "[red]"
        path_width = max([len(p) for p in paths])
        console.print(
            path_color
            + path.rjust(path_width)
            + ":"
            + success_color
            + str(success).rjust(6)
            + " / "
            + processed_events_color
            + str(nEvtProcessed).ljust(6)
            + "[blue]{:4.2f}".format(100.0 * float(success / nEvtProcessed)).rjust(12)
        )
