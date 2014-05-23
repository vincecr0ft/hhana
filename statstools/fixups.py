from .histfactory import process_measurement
from rootpy.stats.histfactory import measurements_from_xml, write_measurement
from rootpy.io import MemFile
import os

from . import log; log = log[__name__]


def find_measurements(path):
    for dirpath, dirnames, filenames in os.walk(path):
        # does this directory contain a workspace?
        if not 'HistFactorySchema.dtd' in filenames:
            continue
        # find the top-level combination XMLs
        for filename in filenames:
            with open(os.path.join(dirpath, filename)) as f:
                if '<Combination' in f.read():
                    yield dirpath, filename


def fix(input, suffix='fixed', verbose=False, **kwargs):
    """
    Traverse all workspaces and apply HSG4 fixes
    """
    if not os.path.isdir(input):
        raise ValueError("input must be an existing directory")
    input = os.path.normpath(input)
    output = input + '_' + suffix
    for dirpath, measurement_file in find_measurements(input):
        path = os.path.join(dirpath, measurement_file)
        log.info("fixing {0} ...".format(path))
        measurements = measurements_from_xml(
            path,
            cd_parent=True,
            collect_histograms=True,
            silence=not verbose)
        for meas in measurements:
            with MemFile():
                fix_measurement(meas, **kwargs)
                write_measurement(meas,
                    output_path=os.path.join(
                        output, dirpath.replace(input, '', 1)[1:]),
                    write_workspaces=True,
                    silence=not verbose)


def fix_measurement(meas, fill_empties=False):
    """
    Apply the HSG4 fixes on a HistFactory::Measurement
    Changes are applied in-place
    """
    process_measurement(meas,
        split_norm_shape=True,
        drop_np_names=["*"],
        drop_np_types=['histosys'],
        drop_np_samples=['Others'],
        prune_overallsys=True,
        prune_overallsys_threshold=0.5, # percent
        uniform_binning=True,
        #fill_empties=True,
        #fill_empties_samples=['Fakes', 'Ztautau']
        )
    #--prune-histosys --prune-histosys-method chi2 --prune-histosys-threshold 0.99 \
    #--prune-histosys-blacklist QCDscale_ggH3in \
    #--rebin-channels channel_vbf_${mass} channel_boosted_${mass} \
    #--rebin 10 from Makefile

    process_measurement(meas,
        drop_np_names=["*TAU_ID*"],
        drop_np_types=["OverallSys"],
        drop_np_samples=['Ztautau'],
        #symmetrize_names=["*JVF*", "*TES_TRUE*", "*EMB_MFS*"],
        symmetrize_names=["*TES_TRUE_MODELING*"],
        symmetrize_types=["overallsys", "histosys"]
        )
    #--smooth-histosys --smooth-histosys-iterations 1 \
    #--smooth-histosys-samples Fakes Ztautau Others "Signal_VBF_*" "Signal_gg_*" \
    #--prune-histosys --prune-histosys-samples Fakes Others Ztautau \
    #--prune-histosys-method max --prune-histosys-threshold 0.1

    if fill_empties:
        process_measurement(meas,
            fill_empties=True,
            fill_empties_samples=['Fakes', 'Ztautau']
            )
