import rootpy
from rootpy.extern import argparse

from .categories import CATEGORIES
from .massregions import DEFAULT_LOW_MASS, DEFAULT_HIGH_MASS
from .variables import VARIABLES
from .regions import REGIONS
from .defaults import FAKES_REGION, TARGET_REGION


class formatter_class(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


def base_parser():
    return argparse.ArgumentParser(formatter_class=formatter_class)


def general_parser(parser=None, multi_years=False, multi_categories=False):
    if parser is None:
        parser = base_parser()
    parser.add_argument('--systematics', action='store_true', default=False,
            help="enable systematics")

    if multi_years:
        parser.add_argument('--years', type=int, default=[2012], choices=(2011, 2012),
                nargs='*',
                help='years')
    else:
        parser.add_argument('--year', type=int, default=2012, choices=(2011, 2012),
                help='the year')

    if multi_categories:
        parser.add_argument('--categories', default=['mva'],
                choices=CATEGORIES.keys(),
                nargs='*',
                help='category definitions')
    else:
        parser.add_argument('--categories', default='mva',
                choices=CATEGORIES.keys(),
                help='category definitions')

    parser.add_argument('--category-names', nargs="+", default=None,
            help='category names')
    parser.add_argument('--controls', default='mva_workspace_controls',
            help='control definitions')
    parser.add_argument('--unblind', action='store_true', default=False,
            help='plot the data in the signal region of the classifier output')
    parser.add_argument('--masses', default='125')
    parser.add_argument('--mixings', type=float, nargs="+", default=[0.0])
    parser.add_argument('--suffix', default=None, nargs='?',
            help='suffix to add to any output files or plots')
    parser.add_argument('--output-suffix', default=None, nargs='?',
            help='suffix to add to any output files or plots')
    parser.add_argument('--systematics-components', default=None,
            help='only include the following systematics in plots Example: '
                 'TES_TRUE_INSITU_UP,QCD_SHAPE_UP')
    return parser


def analysis_parser(parser=None):
    if parser is None:
        parser = base_parser()
    parser.add_argument('--random-mu', action='store_true', default=False,
            help='set mu (signal strength) to a random number')
    parser.add_argument('--mu', default=1., type=float,
            help='set mu (signal strength)')
    parser.add_argument('--no-embedding', action='store_false', default=True,
            dest='embedding',
            help='use ALPGEN Z->tau+tau instead of embedding')
    parser.add_argument('--fakes-region', choices=REGIONS.keys(),
            default=FAKES_REGION,
            help='fakes shape region')
    parser.add_argument('--target-region', choices=REGIONS.keys(),
            default=TARGET_REGION,
            help='target signal region')
    parser.add_argument('--constrain-norms',
            action='store_true', default=False)
    parser.add_argument('--decouple-qcd-shape',
            action='store_true', default=False)
    parser.add_argument('--no-qcd-shape-systematic',
            dest='qcd_shape_systematic',
            action='store_false', default=True)
    parser.add_argument('--no-ggf-weight',
            dest='ggf_weight',
            action='store_false', default=True)
    return parser


def mass_parser(parser=None):
    if parser is None:
        parser = base_parser()
    parser.add_argument('--low-mass-cut', type=int,
            default=DEFAULT_LOW_MASS,
            help='the low mass window cut. '
            'Norms of Z and QCD are fit below this and '
            'the signal region of the classifier output is above this')
    parser.add_argument('--high-mass-cut', type=int,
            default=DEFAULT_HIGH_MASS,
            help='the high mass window cut. '
            'Norms of Z and QCD are fit above this and '
            'the signal region of the classifier output is below this')
    parser.add_argument('--no-sideband-in-control',
            dest='high_sideband_in_control',
            action='store_false',
            default=True,
            help='Exclude the high mass sideband in the mass control and include '
            'it in the signal region')
    return parser


def plotting_parser(parser=None):
    if parser is None:
        parser = base_parser()
    parser.add_argument('--plots', nargs='*',
                        help='only draw these plots. see the keys in variables.py')
    parser.add_argument('--plot-cut', default=None, nargs='?',
                        help='extra cut to be applied on the plots, but excluded from the '
                        'QCD/Z normaliation and training and classifier output')
    parser.add_argument('--plot-expr', default=None, nargs='?',
                        help='expression to plot, instead of predefined ones in variables.py')
    parser.add_argument('--plot-name', default=None, nargs='?',
                        help='name of expr')
    parser.add_argument('--plot-min', type=float, default=0, nargs='?',
                        help='minimum of expr')
    parser.add_argument('--plot-max', type=float, default=1, nargs='?',
                        help='maximum of expr')
    parser.add_argument('--plot-bins', type=int, default=20, nargs='?',
                        help='number of bins to plot expr in')
    parser.add_argument('--no-weight', action='store_true', default=False,
                        help='do not apply correction weights')
    parser.add_argument('--output-formats', default=['png'], nargs='+',
                        choices=('png', 'eps', 'pdf'),
                        help='output formats')
    parser.add_argument('--no-data', action='store_true', default=False,
                        help='do not display data on the plot')
    parser.add_argument('--show-ratio', action='store_true', default=False,
                        help='Draw plot with a ratio plot below the main plot')
    return parser


def get_parser(actions=True):
    parser = general_parser()
    analysis_parser(parser)
    mass_parser(parser)
    plotting_parser(parser)
    if actions:
        parser.add_argument('actions', nargs='*',
            choices=[
                'stability',
                'validate',
                'weights',
                '2d',
                'plotevolving',
                'money',
                'scatter',
                'correlate',
                'evaluate',
                'workspace',
                'ntup',
                'ntuptruth',
                'top10',
                'overlap',
                'cuts_notmva',
                'massplot'],
            default=[],
            help='only perform these actions')
    return parser
