# stdlib imports
import os
import pickle
from math import pi, sqrt

import numpy as np

# rootpy imports
from rootpy.io import root_open
from rootpy.stats import histfactory
from root_numpy import fill_hist

# Higgs cross sections and branching ratios
import yellowhiggs

# local imports
from . import log
from .. import ETC_DIR, CACHE_DIR, DAT_DIR
from ..utils import uniform_hist
from .sample import MC, Signal


TAUTAUHADHADBR = 0.4197744 # = (1. - 0.3521) ** 2


class Higgs(MC, Signal):
    MASSES = range(100, 155, 5)
    MASS_POINTS = MASSES
    MIXINGS = [-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,
               -0.07,-0.05,-0.02, 0.00, 0.02, 0.05, 0.07, 
                0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
    MODES = ['Z', 'W', 'gg', 'VBF']
    MODES_COMBINED = [['Z', 'W'], ['gg'], ['VBF']]
    MODES_CP = [['Z', 'W', 'gg'], ['VBF']]
    MODES_DICT = {
        'gg': ('ggf', 'PowPyth_', 'PowPyth8_AU2CT10_'),
        'VBF': ('vbf', 'PowPyth_', 'PowPyth8_AU2CT10_'),
        'Z': ('zh', 'Pyth8_AU2CTEQ6L1_', 'Pyth8_AU2CTEQ6L1_'),
        'W': ('wh', 'Pyth8_AU2CTEQ6L1_', 'Pyth8_AU2CTEQ6L1_'),
    }
    MODES_WORKSPACE = {
        'gg': 'ggH',
        'VBF': 'VBF',
        'Z': 'ZH',
        'W': 'WH',
    }

    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HSG4Uncertainties
    # TODO: UPDATE
    QCD_SCALE = map(lambda token: token.strip().split(), '''\
    QCDscale_qqH     VBF    rest                   1.020/0.980
    QCDscale_qqH     VBF    boosted                1.014/0.986
    QCDscale_qqH     VBF    vbf                    1.020/0.980
    QCDscale_qqH     VBF    cuts_boosted_tight     1.020/0.980
    QCDscale_qqH     VBF    cuts_boosted_loose     1.020/0.980
    QCDscale_qqH     VBF    cuts_vbf_highdr_tight  1.020/0.980
    QCDscale_qqH     VBF    cuts_vbf_highdr_loose  1.020/0.980
    QCDscale_qqH     VBF    cuts_vbf_lowdr         1.020/0.980
    QCDscale_VH      VH     rest                   1.010/0.990
    QCDscale_VH      VH     boosted                1.041/0.960
    QCDscale_VH      VH     vbf                    1.010/0.990
    QCDscale_VH      VH     cuts_boosted_tight     1.040/0.960
    QCDscale_VH      VH     cuts_boosted_loose     1.040/0.960
    QCDscale_ggH     ggH    rest                   1.070/0.930
    QCDscale_ggH1in  ggH    boosted                1.320/0.760
    QCDscale_ggH1in  ggH    cuts_boosted_tight     1.280/0.780
    QCDscale_ggH1in  ggH    cuts_boosted_loose     1.280/0.780
    QCDscale_ggH2in  ggH    boosted                0.930/1.080
    QCDscale_ggH2in  ggH    vbf                    1.260/0.800
    QCDscale_ggH2in  ggH    cuts_boosted_tight     0.960/1.050
    QCDscale_ggH2in  ggH    cuts_boosted_loose     0.970/1.030
    QCDscale_ggH2in  ggH    cuts_vbf_highdr_tight  1.260/0.790
    QCDscale_ggH2in  ggH    cuts_vbf_highdr_loose  1.260/0.790
    QCDscale_ggH2in  ggH    cuts_vbf_lowdr         1.250/0.800'''.split('\n'))

    UE_UNCERT = map(lambda token: token.strip().split(), '''\
    ATLAS_UE_qq  VBF      vbf                      1.080/0.920
    ATLAS_UE_qq  VBF      boosted                  1.050/0.950
    ATLAS_UE_qq  VBF      cuts_vbf_lowdr           1.110/0.890
    ATLAS_UE_qq  VBF      cuts_vbf_highdr_tight    1.080/0.920
    ATLAS_UE_qq  VBF      cuts_vbf_highdr_loose    1.070/0.930
    ATLAS_UE_qq  VBF      cuts_boosted_tight       1.090/0.910
    ATLAS_UE_qq  VBF      cuts_boosted_loose       1.020/0.980
    ATLAS_UE_gg  ggH      vbf                      1.010/0.990
    ATLAS_UE_gg  ggH      boosted                  1.060/0.940
    ATLAS_UE_gg  ggH      cuts_vbf_lowdr           0.990/1.010
    ATLAS_UE_gg  ggH      cuts_vbf_highdr_tight    0.530/1.470
    ATLAS_UE_gg  ggH      cuts_vbf_highdr_loose    1.140/0.860
    ATLAS_UE_gg  ggH      cuts_boosted_tight       1.010/0.990
    ATLAS_UE_gg  ggH      cuts_boosted_loose       1.160/0.840'''.split('\n'))

    PDF_ACCEPT_NORM_UNCERT = map(lambda token: token.strip().split(), '''\
    pdf_Higgs_qq_ACCEPT  VBF      vbf                      1.010/0.990
    pdf_Higgs_qq_ACCEPT  VBF      boosted                  1.010/0.990
    pdf_Higgs_qq_ACCEPT  VBF      cuts_vbf_lowdr           1.010/0.990
    pdf_Higgs_qq_ACCEPT  VBF      cuts_vbf_highdr_tight    1.010/0.990
    pdf_Higgs_qq_ACCEPT  VBF      cuts_vbf_highdr_loose    1.020/0.980
    pdf_Higgs_qq_ACCEPT  VBF      cuts_boosted_tight       1.030/0.970
    pdf_Higgs_qq_ACCEPT  VBF      cuts_boosted_loose       1.010/0.990
    pdf_Higgs_qq_ACCEPT  VH       vbf                      1.010/0.990
    pdf_Higgs_qq_ACCEPT  VH       boosted                  1.010/0.990
    pdf_Higgs_qq_ACCEPT  VH       cuts_vbf_lowdr           1.010/0.990
    pdf_Higgs_qq_ACCEPT  VH       cuts_vbf_highdr_tight    1.010/0.990
    pdf_Higgs_qq_ACCEPT  VH       cuts_vbf_highdr_loose    1.020/0.980
    pdf_Higgs_qq_ACCEPT  VH       cuts_boosted_tight       1.030/0.970
    pdf_Higgs_qq_ACCEPT  VH       cuts_boosted_loose       1.010/0.990
    pdf_Higgs_gg_ACCEPT  ggH      vbf                      1.050/0.950
    pdf_Higgs_gg_ACCEPT  ggH      boosted                  1.060/0.940
    pdf_Higgs_gg_ACCEPT  ggH      cuts_vbf_lowdr           1.050/0.950
    pdf_Higgs_gg_ACCEPT  ggH      cuts_vbf_highdr_tight    1.050/0.950
    pdf_Higgs_gg_ACCEPT  ggH      cuts_vbf_highdr_loose    1.050/0.950
    pdf_Higgs_gg_ACCEPT  ggH      cuts_boosted_tight       1.060/0.940
    pdf_Higgs_gg_ACCEPT  ggH      cuts_boosted_loose       1.060/0.940'''.split('\n'))

    PDF_ACCEPT_SHAPE_UNCERT = map(lambda token: token.strip().split(), '''\
    pdf_Higgs_qq_ACCEPT  VBF  vbf      h_VBF_vbf_{0}TeV_Up/h_VBF_vbf_{0}TeV_Down
    pdf_Higgs_qq_ACCEPT  VBF  boosted  h_VBF_boosted_{0}TeV_Up/h_VBF_boosted_{0}TeV_Down
    pdf_Higgs_gg_ACCEPT  ggH  vbf      h_gg_vbf_{0}TeV_Up/h_gg_vbf_{0}TeV_Down
    pdf_Higgs_gg_ACCEPT  ggH  boosted  h_gg_boosted_{0}TeV_Up/h_gg_boosted_{0}TeV_Down'''.split('\n'))

    PDF_ACCEPT_file = root_open(
        os.path.join(DAT_DIR, 'ShapeUnc_PDF_hh.root'), 'read')

    NORM_BY_THEORY = True

    def __init__(self, year,
                 mode=None, modes=None,
                 mass=None, masses=None,
                 mixing=None, mixings=None,
                 sample_pattern=None, # i.e. PowhegJimmy_AUET2CT10_ggH{0:d}_tautauInclusive
                 ggf_weight=True,
                 vbf_weight=True,
                 suffix=None,
                 label=None,
                 inclusive_decays=False,
                 SM=False,
                 BSM=False,
                 **kwargs):
        self.inclusive_decays = inclusive_decays
        if masses is None:
            if mass is not None:
                assert mass in Higgs.MASSES
                masses = [mass]
            else:
                # default to 125
                masses = [125]
        else:
            assert len(masses) > 0
            for mass in masses:
                assert mass in Higgs.MASSES
            assert len(set(masses)) == len(masses)
        if mixings is None:
            if mixing is not None:
                assert mixing in Higgs.MIXINGS , [mixing,Higgs.MIXINGS]
                mixings = [mixing]
            else:
                # default to 125
                mixings = [0.0]
        else:
            assert len(mixings) > 0
            for mixing in mixings:
                assert mixing in Higgs.MIXINGS
            assert len(set(mixings)) == len(mixings)

        if modes is None:
            if mode is not None:
                assert mode in Higgs.MODES, mode
                modes = [mode]
            else:
                # default to all modes
                modes = Higgs.MODES
        else:
            assert len(modes) > 0
            for mode in modes:
                assert mode in Higgs.MODES
            assert len(set(modes)) == len(modes)

        self.SM=SM
        self.BSM=BSM

        name = 'Signal'
        isNonVBF=False
        str_mode = ''
        if len(modes) == 1:
            str_mode = modes[0]
            name += '_%s' % str_mode
        elif len(modes) == 2 and set(modes) == set(['W', 'Z']):
            str_mode = 'V'
            name += '_%s' % str_mode
        elif len(modes) == 3 and set(modes) == set(['W', 'Z', 'gg']):
            str_mode = 'gg+V'
            name += '_%s' % str_mode
            isNonVBF=True

        str_mixings = ''
        if len(mixings) == 1:
            str_mixings = "%0.2f" % mixings[0]
            name += "_%s" % str_mixings
            if mixings[0]==0.0 and BSM:
                name+='_copy'
        if label is None:
            label = '%s#font[52]{H}(%s)#rightarrow#tau#tau' % (
                str_mode, str_mixings)
            if isNonVBF:
                label = '%s#font[52]{H}(%s)#rightarrow#tau#tau' % (
                    str_mode, '125')



        if year == 2011:
            if suffix is None:
                suffix = '.mc11c'
            generator_index = 1
        elif year == 2012:
            if suffix is None:
                suffix = '.mc12a'
            generator_index = 2
        else:
            raise ValueError('No Higgs defined for year %d' % year)
        #####################
        self.samples = []
        self.masses = []
        self.mixings = []
        self.modes = []

        if sample_pattern is not None:
            assert len(modes) == 1
            for mass in masses:
                self.masses.append(mass)
                self.modes.append(modes[0])
                self.samples.append(sample_pattern.format(mass) + '.' + suffix)
        else:
            for mode in modes:
                generator = Higgs.MODES_DICT[mode][generator_index]
                for mass in masses:
                    self.samples.append('%s%sH%d_tautauhh%s' % (
                            generator, mode, mass, suffix))
                    self.masses.append(mass)
                    self.modes.append(mode)
        for mixing in mixings:
            self.mixings.append(mixing)

        if len(self.modes) == 1:
            self.mode = self.modes[0]
        else:
            self.mode = None
        if len(self.masses) == 1:
            self.mass = self.masses[0]
        else:
            self.mass = None
        if len(self.mixings) == 1:
            self.mixing = self.mixings[0]
        else:
            self.mixing = None

        #####################

        self.ggf_weight = ggf_weight
        self.ggf_weight_field = 'ggf_weight'
        self.vbf_weight = vbf_weight
        self.vbf_weight_field = 'vbf_weight'
        # use separate signal files by default
        kwargs.setdefault('student', 'hhskim_signal')
        super(Higgs, self).__init__(
            year=year, label=label, name=name, **kwargs)

    #def weight_systematics(self):
    #    systematics = super(Higgs, self).weight_systematics()
    #    if self.ggf_weight:
    #        systematics.update({
    #            'QCDscale_ggH1in'})
    #    return systematics

    def weight_fields(self):
        fields = super(Higgs, self).weight_fields()
        if self.ggf_weight:
            fields.append(self.ggf_weight_field)
        if self.vbf_weight:
            fields.append(self.vbf_weight_field)
        return fields

    def histfactory(self, sample, category, systematics=False,
                    rec=None, weights=None, mva=False,
                    uniform=False, nominal=None):
        if not systematics:
            return
        if len(self.modes) != 1:
            raise TypeError(
                'histfactory sample only valid for single production mode')
        if len(self.masses) != 1:
            raise TypeError(
                'histfactory sample only valid for single mass point')

        # isolation systematic
#        sample.AddOverallSys(
#            'ATLAS_ANA_HH_{0:d}_Isolation'.format(self.year),
#            1. - 0.06,
#            1. + 0.06)

        mode = self.modes[0]

        if mode in ('Z', 'W'):
            _uncert_mode = 'VH'
        else:
            _uncert_mode = self.MODES_WORKSPACE[mode]

        if self.year == 2011:
            energy = 7
        elif self.year == 2012:
            energy = 8
        else:
            raise ValueError(
                "collision energy is unknown for year {0:d}".format(self.year))

        # QCD_SCALE
        for qcd_scale_term, qcd_scale_mode, qcd_scale_category, values in self.QCD_SCALE:
            if qcd_scale_mode == _uncert_mode and qcd_scale_category == category.name:
                high, low = map(float, values.split('/'))
                sample.AddOverallSys(qcd_scale_term, low, high)

        # UE UNCERTAINTY
        for ue_term, ue_mode, ue_category, values in self.UE_UNCERT:
            if ue_mode == _uncert_mode and ue_category == category.name:
                high, low = map(float, values.split('/'))
                sample.AddOverallSys(ue_term, low, high)

        # PDF ACCEPTANCE UNCERTAINTY (OverallSys)
        for pdf_term, pdf_mode, pdf_category, values in self.PDF_ACCEPT_NORM_UNCERT:
            if pdf_mode == _uncert_mode and pdf_category == category.name:
                high, low = map(float, values.split('/'))
                sample.AddOverallSys(pdf_term, low, high)

        sample_nom = sample.hist

        # PDF ACCEPTANCE UNCERTAINTY (HistoSys) ONLY FOR MVA
        if mva:
            for pdf_term, pdf_mode, pdf_category, hist_names in self.PDF_ACCEPT_SHAPE_UNCERT:
                if pdf_mode == _uncert_mode and pdf_category == category.name:
                    high_name, low_name = hist_names.format(energy).split('/')
                    high_shape, low_shape = self.PDF_ACCEPT_file[high_name], self.PDF_ACCEPT_file[low_name]
                    if len(high_shape) != len(sample.hist):
                        log.warning("skipping pdf acceptance shape systematic "
                                    "since histograms are not compatible")
                        continue
                    high = sample_nom.Clone(shallow=True, name=sample_nom.name + '_{0}_UP'.format(pdf_term))
                    low = sample_nom.Clone(shallow=True, name=sample_nom.name + '_{0}_DOWN'.format(pdf_term))
                    high *= high_shape
                    low *= low_shape
                    histsys = histfactory.HistoSys(
                        pdf_term, low=low, high=high)
                    sample.AddHistoSys(histsys)

        #mixing Norms
        if self.SM:
            print 'adding norm factor'
            sample.AddNormFactor('ATLAS_epsilon', 1., -200., 200., False)
        elif self.BSM:
            print 'adding norm factor'
            sample.AddNormFactor('ATLAS_epsilon_rejected', 1., -200., 200., False)
        else:
            print 'no norms for ',self.name



        # BR_tautau
        _, (br_up, br_down) = yellowhiggs.br(
            self.mass, 'tautau', error_type='factor')
        sample.AddOverallSys('ATLAS_BR_tautau', br_down, br_up)

        # <NormFactor Name="mu_BR_tautau" Val="1" Low="0" High="200" />
        sample.AddNormFactor('mu_BR_tautau', 1., 0., 200., True)

        #mu_XS[energy]_[mode]
        #_, (xs_up, xs_down) = yellowhiggs.xs(
        #    energy, self.mass, self.MODES_DICT[self.mode][0],
        #    error_type='factor')
        #sample.AddOverallSys(
        #    'mu_XS{0:d}_{1}'.format(energy, self.MODES_WORKSPACE[self.mode]),
        #    xs_down, xs_up)
        sample.AddNormFactor(
            'mu_XS{0:d}_{1}'.format(energy, self.MODES_WORKSPACE[self.mode]),
            1., 0., 200., True)

        # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HSG4Uncertainties
        # pdf uncertainty
        if mode == 'gg':
            if energy == 8:
                sample.AddOverallSys('pdf_Higgs_gg', 0.93, 1.08)
            else: # 7 TeV
                sample.AddOverallSys('pdf_Higgs_gg', 0.92, 1.08)
        else:
            if energy == 8:
                sample.AddOverallSys('pdf_Higgs_qq', 0.97, 1.03)
            else: # 7 TeV
                sample.AddOverallSys('pdf_Higgs_qq', 0.98, 1.03)

        # EWK NLO CORRECTION FOR VBF ONLY
        if mode == 'VBF':
            sample.AddOverallSys('NLO_EW_Higgs', 0.98, 1.02)

        # QCDscale_ggH3in HistoSys ONLY FOR MVA
        # also see ggH3in script
        if mva and mode == 'gg' and category.name == 'vbf':
            Rel_Error_2j = 0.215
            Error_exc = 0.08613046469238815 # Abs error on the exclusive xsec
            xsec_exc = 0.114866523583739 # Exclusive Xsec
            Error_3j = sqrt(Error_exc**2 - (Rel_Error_2j*xsec_exc)**2)
            rel_error = Error_3j / xsec_exc

            dphi = rec['true_dphi_jj_higgs_no_overlap']
            scores = rec['classifier']

            idx_2j = ((pi - dphi) < 0.2) & (dphi >= 0)
            idx_3j = ((pi - dphi) >= 0.2) & (dphi >= 0)

            # get normalization factor
            dphi_2j = weights[idx_2j].sum()
            dphi_3j = weights[idx_3j].sum()

            weight_up = np.ones(len(weights))
            weight_dn = np.ones(len(weights))

            weight_up[idx_2j] -= (dphi_3j / dphi_2j) * rel_error
            weight_dn[idx_2j] += (dphi_3j / dphi_2j) * rel_error

            weight_up[idx_3j] += rel_error
            weight_dn[idx_3j] -= rel_error

            weight_up *= weights
            weight_dn *= weights

            up_hist = nominal.clone(shallow=True, name=sample_nom.name + '_QCDscale_ggH3in_UP')
            up_hist.Reset()
            dn_hist = nominal.clone(shallow=True, name=sample_nom.name + '_QCDscale_ggH3in_DOWN')
            dn_hist.Reset()

            fill_hist(up_hist, scores, weight_up)
            fill_hist(dn_hist, scores, weight_dn)

            if uniform:
                up_hist = uniform_hist(up_hist)
                dn_hist = uniform_hist(dn_hist)

            shape = histfactory.HistoSys('QCDscale_ggH3in',
                low=dn_hist,
                high=up_hist)
            norm, shape = histfactory.split_norm_shape(shape, sample_nom)
            sample.AddHistoSys(shape)

    def xsec_kfact_effic(self, isample):
        # use yellowhiggs for cross sections
        xs, _ = yellowhiggs.xsbr(
            self.energy, self.masses[isample],
            Higgs.MODES_DICT[self.modes[isample]][0], 'tautau')
        log.debug("{0} {1} {2} {3} {4} {5}".format(
            self.samples[isample],
            self.masses[isample],
            self.modes[isample],
            Higgs.MODES_DICT[self.modes[isample]][0],
            self.energy,
            xs))
        if not self.inclusive_decays:
            xs *= TAUTAUHADHADBR
        kfact = 1.
        effic = 1.
        return xs, kfact, effic


class InclusiveHiggs(MC, Signal):
    # for overlap study
    SAMPLES = {
        'ggf': 'PowPyth8_AU2CT10_ggH125p5_inclusive.mc12b',
        'vbf': 'PowPyth8_AU2CT10_VBFH125p5_inclusive.mc12b',
        'zh': 'Pyth8_AU2CTEQ6L1_ZH125p5_inclusive.mc12b',
        'wh': 'Pyth8_AU2CTEQ6L1_WH125p5_inclusive.mc12b',
        'tth': 'Pyth8_AU2CTEQ6L1_ttH125p5_inclusive.mc12b',
    }

    def __init__(self, mode=None, **kwargs):
        self.energy = 8
        if mode is not None:
            self.modes = [mode]
            self.masses = [125]
            self.samples = [self.SAMPLES[mode]]
        else:
            self.masses = [125] * 5
            self.modes = ['ggf', 'vbf', 'zh', 'wh', 'tth']
            self.samples = [
                'PowPyth8_AU2CT10_ggH125p5_inclusive.mc12b',
                'PowPyth8_AU2CT10_VBFH125p5_inclusive.mc12b',
                'Pyth8_AU2CTEQ6L1_ZH125p5_inclusive.mc12b',
                'Pyth8_AU2CTEQ6L1_WH125p5_inclusive.mc12b',
                'Pyth8_AU2CTEQ6L1_ttH125p5_inclusive.mc12b',
            ]
        super(InclusiveHiggs, self).__init__(
            year=2012, name='Signal', label='Signal',
            ntuple_path='ntuples/prod_v29',
            student='hhskim_overlap',
            **kwargs)

    def xsec_kfact_effic(self, isample):
        # use yellowhiggs for cross sections
        xs, _ = yellowhiggs.xs(
            self.energy, self.masses[isample], self.modes[isample])
        log.debug("{0} {1} {2} {3} {4}".format(
            self.samples[isample],
            self.masses[isample],
            self.modes[isample],
            self.energy,
            xs))
        kfact = 1.
        effic = 1.
        return xs, kfact, effic
