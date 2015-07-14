# std lib imports
import os
import sys
import pickle
from operator import add, itemgetter
from collections import namedtuple

# numpy imports
import numpy as np
from numpy.lib import recfunctions

# rootpy imports
import ROOT
from rootpy.plotting import Hist, Hist2D, Canvas, HistStack
from rootpy.plotting.hist import _HistBase
from rootpy.tree import Tree, Cut
from rootpy.stats import histfactory
from rootpy import asrootpy

# root_numpy imports
from root_numpy import rec2array, stack, fill_hist

# higgstautau imports
from higgstautau import samples as samples_db

# local imports
from . import log; log = log[__name__]
from .. import variables
from .. import NTUPLE_PATH, DEFAULT_STUDENT, ETC_DIR, CACHE_DIR
from ..utils import print_hist, ravel_hist, uniform_hist
from ..classify import histogram_scores, Classifier
from ..regions import REGIONS
from ..systematics import (
    get_systematics, SYSTEMATICS_BY_WEIGHT,
    iter_systematics, systematic_name)
from ..lumi import LUMI, get_lumi_uncert
from .db import DB, TEMPFILE, get_file
from ..cachedtable import CachedTable
from ..variables import get_binning, get_scale
from .calcOptObs import observize,new_observize,dphize,weight_vize,new_weight_vize

BCH_UNCERT = pickle.load(open(os.path.join(CACHE_DIR, 'bch_cleaning.cache')))


class Dataset(namedtuple('Dataset',
                         ('ds', 'tables', 'events',
                          'xs', 'kfact', 'effic'))):
    @property
    def name(self):
        return self.ds.name


def get_workspace_np_name(sample, syst, year):
    """
    HSG4 naming convention for NPs in the workspaces
    """
    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsPropertiesNuisanceParameterNames
    npname = 'ATLAS_{0}_{1:d}'.format(syst, year)
    # special cases
    npname = npname.replace('JES_Detector_{0}'.format(year),
                            'JES_{0}_Detector1'.format(year))
    npname = npname.replace('JES_EtaMethod_{0}'.format(year),
                            'JES_{0}_Eta_StatMethod'.format(year))
    npname = npname.replace('JES_EtaModelling_{0}'.format(year),
                            'JES_Eta_Modelling')
    npname = npname.replace('JES_FlavComp_TAU_G_{0}'.format(year),
                            'JES_FlavComp_TAU_G')
    npname = npname.replace('JES_FlavComp_TAU_Q_{0}'.format(year),
                            'JES_FlavComp_TAU_Q')
    npname = npname.replace('JES_FlavResp_{0}'.format(year),
                            'JES_FlavResp')
    npname = npname.replace('JES_Modelling_{0}'.format(year),
                            'JES_{0}_Modelling1'.format(year))
    npname = npname.replace('JES_PURho_TAU_GG_{0}'.format(year),
                            'JES_{0}_PileRho_TAU_GG'.format(year))
    npname = npname.replace('JES_PURho_TAU_QG_{0}'.format(year),
                            'JES_{0}_PileRho_TAU_QG'.format(year))
    npname = npname.replace('JES_PURho_TAU_QQ_{0}'.format(year),
                            'JES_{0}_PileRho_TAU_QQ'.format(year))
    npname = npname.replace('FAKERATE', 'TAU_JFAKE')
    npname = npname.replace('MET_RESOSOFTTERMS',
                            'MET_RESOSOFT')
    npname = npname.replace('MET_SCALESOFTTERMS',
                            'MET_SCALESOFT')
    from .ztautau import Embedded_Ztautau
    if isinstance(sample, Embedded_Ztautau):
        npname = npname.replace('TRIGGER', 'TRIGGER_EMB_HH')
    else:
        npname = npname.replace('TRIGGER', 'TRIGGER_HH')
    # Decorrelate embedding NPs
    # * Decorrelate the NP between 2011 and 2012 for MFS because the cell
    #   subtraction yield was changed from 30% in 2011 to 20% in 2012.
    npname = npname.replace('ISOL', 'ANA_EMB_ISOL')
    npname = npname.replace('MFS', 'ANA_EMB_MFS')
    # correlate across years:
    npname = npname.replace('JER_{0}'.format(year), 'JER')
    return npname


class Sample(object):

    def weight_fields(self):
        return []

    def corrections(self, rec):
        return []

    def __init__(self, year, scale=1., cuts=None,
                 ntuple_path=NTUPLE_PATH,
                 student=DEFAULT_STUDENT,
                 force_reopen=False,
                 trigger=True,
                 name='Sample',
                 label='Sample',
                 **hist_decor):
        self.year = year
        if year == 2011:
            self.energy = 7
        else:
            self.energy = 8
        self.scale = scale
        if cuts is None:
            self._cuts = Cut()
        else:
            self._cuts = cuts
        self.ntuple_path = ntuple_path
        self.student = student
        self.force_reopen = force_reopen
        self.name = name
        self.label = label
        self.hist_decor = hist_decor
        if 'fillstyle' not in hist_decor:
            self.hist_decor['fillstyle'] = 'solid'
        self.trigger = trigger

    def optimate(self, rec):
        pt1 = rec2array(rec[['jet1_pt']])
        pt2 = rec2array(rec[['jet2_pt']])
        pt3 = rec2array(rec[['jet3_pt']])

        eta1 = rec2array(rec[['jet1_eta']])
        eta2 = rec2array(rec[['jet2_eta']])
        eta3 = rec2array(rec[['jet3_eta']])

        phi1 = rec2array(rec[['jet1_phi']])
        phi2 = rec2array(rec[['jet2_phi']])
        phi3 = rec2array(rec[['jet3_phi']])

        tpt1 = rec2array(rec[['mmc1_resonance_pt']])
        teta1 = rec2array(rec[['mmc1_resonance_eta']])
        tphi1 = rec2array(rec[['mmc1_resonance_phi']])
        tm1 = rec2array(rec[['mmc1_resonance_m']])
        numJets = rec2array(rec[['numJets']])

        if phi1.size==phi2.size and phi1.size>0:
            Opt_vec = new_observize(pt1,pt2,pt3,eta1,eta2,eta3,phi1,phi2,phi3,tpt1,teta1,tphi1,tm1,numJets)
        else:
            Opt_vec = np.zeros(pt1.size)
        return Opt_vec

    def dphate(self, rec):
        pt1 = rec2array(rec[['jet1_pt']])
        pt2 = rec2array(rec[['jet2_pt']])
        pt3 = rec2array(rec[['jet3_pt']])

        eta1 = rec2array(rec[['jet1_eta']])
        eta2 = rec2array(rec[['jet2_eta']])
        eta3 = rec2array(rec[['jet3_eta']])

        phi1 = rec2array(rec[['jet1_phi']])
        phi2 = rec2array(rec[['jet2_phi']])
        phi3 = rec2array(rec[['jet3_phi']])

        numJets = rec2array(rec[['numJets']])

        if phi1.size==phi2.size and phi1.size>0:
            dphi_vec = dphize(pt1,pt2,pt3,eta1,eta2,eta3,phi1,phi2,phi3,numJets)
        else:
            dphi_vec = np.zeros(pt1.size)
        return dphi_vec

    def optiweight(self, d_tilde, rec, weight):
        pt1 = rec2array(rec[['jet1_pt']])
        pt2 = rec2array(rec[['jet2_pt']])
        pt3 = rec2array(rec[['jet3_pt']])

        eta1 = rec2array(rec[['jet1_eta']])
        eta2 = rec2array(rec[['jet2_eta']])
        eta3 = rec2array(rec[['jet3_eta']])

        phi1 = rec2array(rec[['jet1_phi']])
        phi2 = rec2array(rec[['jet2_phi']])
        phi3 = rec2array(rec[['jet3_phi']])

        m1 = rec2array(rec[['jet1_m']])
        m2 = rec2array(rec[['jet2_m']])
        m3 = rec2array(rec[['jet3_m']])

        tpt1 = rec2array(rec[['mmc1_resonance_pt']])
        teta1 = rec2array(rec[['mmc1_resonance_eta']])
        tphi1 = rec2array(rec[['mmc1_resonance_phi']])
        tm1 = rec2array(rec[['mmc1_resonance_m']])

        numJets = rec2array(rec[['numJets']])

        Opt_weight = weight_vize(d_tilde,pt1,pt2,pt3,eta1,eta2,eta3,phi1,phi2,phi3,m1,m2,m3,tpt1,teta1,tphi1,tm1,numJets,weight)

        if np.isnan(np.sum(Opt_weight)):
            print "nooooo weights returned"
            Opt_weight=np.zeros(pt1.size)
        return Opt_weight

    def new_optiweight(self, d_tilde, rec, weight):
        pt1 = rec2array(rec[['parton1_pt']])
        pt2 = rec2array(rec[['parton2_pt']])
        pt3 = rec2array(rec[['parton3_pt']])

        eta1 = rec2array(rec[['parton1_eta']])
        eta2 = rec2array(rec[['parton2_eta']])
        eta3 = rec2array(rec[['parton3_eta']])

        phi1 = rec2array(rec[['parton1_phi']])
        phi2 = rec2array(rec[['parton2_phi']])
        phi3 = rec2array(rec[['parton3_phi']])

        tpt1 = rec2array(rec[['higgs_pt']])
        teta1 = rec2array(rec[['higgs_eta']])
        tphi1 = rec2array(rec[['higgs_phi']])
        tm1 = rec2array(rec[['higgs_m']])

        partin1 = rec2array(rec[['partin1_pdgId']])
        partin2 = rec2array(rec[['partin2_pdgId']])
        parton1 = rec2array(rec[['parton1_pdgId']])
        parton2 = rec2array(rec[['parton2_pdgId']])
        parton3 = rec2array(rec[['parton3_pdgId']])

        x1_true = rec2array(rec[['mcevent_pdf_x1_0']])
        x2_true = rec2array(rec[['mcevent_pdf_x2_0']])

        numJets = rec2array(rec[['numJets']])

        Opt_weight = new_weight_vize(d_tilde,pt1,pt2,pt3,eta1,eta2,eta3,phi1,phi2,phi3,tpt1,teta1,tphi1,tm1,partin1,partin2,parton1,parton2,parton3,numJets,weight)

        if np.isnan(np.sum(Opt_weight)):
            print "nooooo weights returned"
            Opt_weight=np.zeros(pt1.size)
        return Opt_weight


    def decorate(self, name=None, label=None, **hist_decor):
        if name is not None:
            self.name = name
        if label is not None:
            self.label = label
        if hist_decor:
            self.hist_decor.update(hist_decor)
        return self

    def get_field_hist(self, vars, category, templates=None):
        """
        retrieve a dictionnary of histograms for the requested
        variables in the given category
        ------
        Parameters:
        - vars: dictionnary of variables (see variables.py)
        - category: Analysis category see categories/*
        - template: dictionnary of Histograms. If specified used those
        histograms as templates to retrieve the various fields. If not
        specified, take the default binning specified in variables.py
        """

        field_hist = {}
        field_scale = {}
        for field, var_info in vars.items():
            if templates is not None and field in templates:
                field_hist[field] = templates[field].Clone(
                    title=self.label, **self.hist_decor)
                continue
            elif isinstance(var_info, _HistBase):
                field_hist[field] = var_info.Clone(
                    title=self.label, **self.hist_decor)
                continue
            bins = get_binning(field, category, self.year)
            #bins = var_info['binning']
            if isinstance(bins, list):
                hist = Hist(bins,
                    title=self.label,
                    type='D',
                    **self.hist_decor)
            elif isinstance(bins, tuple):
                #_range = var_info['range']
                # if isinstance(_range, dict):
                #    _range = _range.get(category.name.upper(), _range[None])
                if len(bins) == 3:
                    bins, min, max = bins
                else:
                    min, max = bins
                hist = Hist(bins, min, max,
                    title=self.label,
                    type='D',
                    **self.hist_decor)
            field_hist[field] = hist
            if 'scale' in var_info and var_info['scale'] != 1.:
                field_scale[field] = var_info['scale']
        return field_hist, field_scale

    def get_hist_array(self,
                       field_hist_template,
                       category, region,
                       cuts=None,
                       clf=None,
                       scores=None,
                       min_score=None,
                       max_score=None,
                       systematics=False,
                       systematics_components=None,
                       suffix=None,
                       field_scale=None,
                       weight_hist=None,
                       weighted=True,
                       bootstrap_data=False):

        do_systematics = (isinstance(self, SystematicsSample)
                          and self.systematics
                          and systematics)
        if do_systematics and systematics_components is None:
            systematics_components = self.systematics_components()

        histname = 'hh_category_{0}_{1}'.format(category.name, self.name)
        if suffix is not None:
            histname += suffix

        field_hist = {}
        for field, hist in field_hist_template.items():
            if isinstance(field, basestring):
                new_hist = hist.Clone(name=histname + '_{0}'.format(field))
            else:
                new_hist = hist.Clone(name=histname)
            new_hist.Reset()
            new_hist.decorate(**self.hist_decor)
            new_hist.title = self.label
            field_hist[field] = new_hist

        if (scores is not None and len(field_hist) == 1
                and 'classifier' in field_hist):
            # special case where we are only filling scores
            # and the scores have been supplied
            from ..classify import histogram_scores
            hist = field_hist['classifier']
            histogram_scores(hist, scores,
                             min_score=min_score,
                             max_score=max_score,
                             inplace=True)
            return field_hist
        rec, weights = self.draw_array(
            field_hist, category, region,
            cuts=cuts,
            weighted=weighted,
            field_scale=field_scale,
            weight_hist=weight_hist,
            clf=clf,
            scores=scores,
            min_score=min_score,
            max_score=max_score,
            systematics=do_systematics,
            systematics_components=systematics_components,
            bootstrap_data=bootstrap_data)

        return field_hist, rec, weights

    def get_histfactory_sample(self, hist_template, expr_or_clf,
                               category, region, **kwargs):
        return self.get_histfactory_sample_array(
            {expr_or_clf: hist_template}, category, region,
            **kwargs)[expr_or_clf]

    def get_histfactory_sample_array(self,
                                     field_hist_template,
                                     category, region,
                                     cuts=None,
                                     clf=None,
                                     scores=None,
                                     min_score=None,
                                     max_score=None,
                                     systematics=False,
                                     suffix=None,
                                     field_scale=None,
                                     weight_hist=None,
                                     weighted=True,
                                     no_signal_fixes=False,
                                     bootstrap_data=False,
                                     ravel=True,
                                     uniform=False,
                                     mva=False):

        from .data import Data
        from .qcd import QCD
        from .others import Others
        log.info("creating histfactory samples for {0}".format(self.name))
        field_hist, rec, weights = self.get_hist_array(
            field_hist_template,
            category, region,
            cuts=cuts,
            clf=clf,
            scores=scores,
            min_score=min_score,
            max_score=max_score,
            systematics=systematics,
            suffix=suffix,
            field_scale=field_scale,
            weight_hist=weight_hist,
            weighted=weighted,
            bootstrap_data=bootstrap_data)
        do_systematics = (not isinstance(self, Data)
                          and self.systematics
                          and systematics)
        if do_systematics and isinstance(self, QCD) and self.shape_systematic:
            qcd_shapes = self.get_shape_systematic_array(
                field_hist,
                category, region,
                cuts=cuts,
                clf=clf,
                min_score=min_score,
                max_score=max_score,
                suffix=suffix,
                field_scale=field_scale,
                weight_hist=weight_hist,
                weighted=weighted)
        if isinstance(self, Data):
            sample_cls = histfactory.Data
        else:
            sample_cls = histfactory.Sample
        field_samples = {}
        for field, hist in field_hist.items():
            sample = sample_cls(self.name)
            nominal_hist = hist.Clone()
            if ravel:
                # convert to 1D if 2D (also handles systematics if present)
                hist = ravel_hist(hist)
            if uniform:
                hist = uniform_hist(hist)
            #print_hist(hist)
            # set the nominal histogram
            sample.hist = hist
            # add systematics samples
            if do_systematics:
                SYSTEMATICS = get_systematics(self.year)
                for sys_component in self.systematics_components():
                    terms = SYSTEMATICS[sys_component]
                    if len(terms) == 1:
                        up_term = terms[0]
                        hist_up = hist.systematics[up_term]
                        # use nominal hist for "down" side
                        hist_down = hist.Clone(name=hist.name + '_copy')
                    else:
                        up_term, down_term = terms
                        hist_up = hist.systematics[up_term]
                        hist_down = hist.systematics[down_term]
                    if sys_component == 'JES_FlavComp':
                        if ((isinstance(self, Signal) and self.mode == 'gg') or
                             isinstance(self, Others)):
                            sys_component += '_TAU_G'
                        else:
                            sys_component += '_TAU_Q'
                    elif sys_component == 'JES_PURho':
                        if isinstance(self, Others):
                            sys_component += '_TAU_QG'
                        elif isinstance(self, Signal):
                            if self.mode == 'gg':
                                sys_component += '_TAU_GG'
                            else:
                                sys_component += '_TAU_QQ'
                    npname = get_workspace_np_name(self, sys_component, self.year)
                    histsys = histfactory.HistoSys(
                        npname,
                        low=hist_down,
                        high=hist_up)
                    sample.AddHistoSys(histsys)
                if isinstance(self, QCD) and self.shape_systematic:
                    high, low = qcd_shapes[field]
                    if ravel:
                        low = ravel_hist(low)
                        high = ravel_hist(high)
                    if uniform:
                        low = uniform_hist(low)
                        high = uniform_hist(high)
                    npname = 'ATLAS_ANA_HH_{0:d}_QCD'.format(self.year)
                    if self.decouple_shape or not self.coherent_shape:
                        npname += '_{0}'.format(category.name)
                    if self.coherent_shape:
                        histsys = histfactory.HistoSys(
                            npname, low=low, high=high)
                        sample.AddHistoSys(histsys)
                    else:
                        shapesys = histfactory.ShapeSys(npname)
                        high -= hist
                        high /= hist
                        for bin in high.bins():
                            bin.value = abs(bin.value)
                        shapesys.hist = high
                        sample.AddShapeSys(shapesys)
            if isinstance(self, Signal):
                log.info('got signal sample')
#                sample.AddNormFactor('SigXsecOverSM', 0., -200., 200., False)
                if self.SM:
                    log.info('it is standard model')
                    sample.AddNormFactor('ATLAS_epsilon', 1., -200., 200., False)               
                elif self.BSM:
                    log.info('it is NOT standard model')
                    sample.AddNormFactor('ATLAS_epsilon_rejected', 1., -200., 200., False)
                else:
                    log.info('it is NOT either standard model or bsm')
            if isinstance(self, Background): #formerly elif
                # only activate stat error on background samples
                sample.ActivateStatError()
            if not isinstance(self, Data):
                norm_by_theory = getattr(self, 'NORM_BY_THEORY', True)
                sample.SetNormalizeByTheory(norm_by_theory)
                if norm_by_theory:
                    lumi_uncert = get_lumi_uncert(self.year)
                    lumi_sys = histfactory.OverallSys(
                        'ATLAS_LUMI_{0:d}'.format(self.year),
                        high=1. + lumi_uncert,
                        low=1. - lumi_uncert)
                    sample.AddOverallSys(lumi_sys)
                    if (self.year == 2012 and
                            do_systematics and category.name in BCH_UNCERT):
                        bch_uncert = BCH_UNCERT[category.name]
                        bch_sys = histfactory.OverallSys(
                            'ATLAS_BCH_Cleaning',
                            high=1. + bch_uncert,
                            low=1. - bch_uncert)
                        sample.AddOverallSys(bch_sys)
            # HACK: disable calling this on signal for now since while plotting
            # we only want to show the combined signal but in the histfactory
            # method we require only a single mode
            if hasattr(self, 'histfactory') and not (
                    isinstance(self, Signal) and no_signal_fixes):
                # perform sample-specific items
                self.histfactory(sample, category,
                                 systematics=do_systematics,
                                 rec=rec, weights=weights,
                                 mva=mva, uniform=uniform,
                                 nominal=nominal_hist)
            field_samples[field] = sample
        return field_samples

    def partitioned_records(self,
              category=None,
              region=None,
              fields=None,
              cuts=None,
              include_weight=True,
              systematic='NOMINAL',
              key=None,
              num_partitions=2,
              return_idx=False):
        """
        Partition sample into num_partitions chunks of roughly equal size
        assuming no correlation between record index and field values.
        """
        partitions = []
        for start in range(num_partitions):
            if key is None:
                # split by index
                log.info("splitting records by index parity")
                recs = self.records(
                    category=category,
                    region=region,
                    fields=fields,
                    include_weight=include_weight,
                    cuts=cuts,
                    systematic=systematic,
                    return_idx=return_idx,
                    start=start,
                    step=num_partitions)
            else:
                # split by field values modulo the number of partitions
                partition_cut = Cut('((abs({0})%{1})>={2})&&((abs({0})%{1})<{3})'.format(
                    key, num_partitions, start, start + 1))
                log.info(
                    "splitting records by key parity: {0}".format(
                        partition_cut))
                recs = self.records(
                    category=category,
                    region=region,
                    fields=fields,
                    include_weight=include_weight,
                    cuts=partition_cut & cuts,
                    systematic=systematic,
                    return_idx=return_idx)
            if return_idx:
                partitions.append(recs)
            else:
                partitions.append(np.hstack(recs))
        return partitions

    def merged_records(self,
                       category=None,
                       region=None,
                       fields=None,
                       cuts=None,
                       clf=None,
                       clf_name='classifier',
                       scores=None,
                       include_weight=True,
                       systematic='NOMINAL'):
        recs = self.records(
            category=category,
            region=region,
            fields=fields,
            include_weight=include_weight,
            cuts=cuts,
            systematic=systematic)
        if include_weight and fields is not None:
            if 'weight' not in fields:
                fields = list(fields) + ['weight']
            if 'o1' not in fields:
                fields = list(fields) + ['o1']
            if 'dphi' not in fields:
                fields = list(fields) + ['dphi']

        rec = stack(recs, fields=fields)
        if clf is not None or scores is not None:
            if scores is None:
                scores, _ = clf.classify(
                    self, category, region,
                    cuts=cuts, systematic=systematic)
            elif isinstance(scores, dict):
                # ignore weights
                scores = scores[systematic][0]
            elif isinstance(scores, tuple):
                # ignore weights
                scores = scores[0]
            rec = recfunctions.rec_append_fields(rec,
                names=clf_name,
                data=scores,
                dtypes='f4')
        return rec

    def array(self, *args, **kwargs):
        return rec2array(self.merged_records(*args, **kwargs))

    def weights(self, systematic='NOMINAL'):
        weight_fields = self.weight_fields()
        if isinstance(self, SystematicsSample):
            systerm, variation = \
                SystematicsSample.get_sys_term_variation(systematic)
            for term, variations in self.weight_systematics().items():
                # handle cases like TAU_ID and TAU_ID_STAT
                if (systerm is not None
                    and systerm.startswith(term)
                    and systematic[0][len(term) + 1:] in variations):
                    weight_fields += variations[systematic[0][len(term) + 1:]]
                elif term == systerm:
                    weight_fields += variations[variation]
                else:
                    weight_fields += variations['NOMINAL']
        # HACK
        if not self.trigger and 'tau1_trigger_sf' in weight_fields:
            log.info("replacing trigger_sf with trigger_eff")
            weight_fields.remove('tau1_trigger_sf')
            weight_fields.remove('tau2_trigger_sf')
            weight_fields.extend(['tau1_trigger_eff', 'tau2_trigger_eff'])
        return weight_fields

    def cuts(self, category=None, region=None, systematic='NOMINAL', **kwargs):
        cuts = Cut(self._cuts)
        if category is not None:
            cuts &= category.get_cuts(self.year, **kwargs)
        if region is not None:
            cuts &= REGIONS[region]
        if self.trigger:
            cuts &= Cut('trigger')
        if isinstance(self, SystematicsSample):
            systerm, variation = SystematicsSample.get_sys_term_variation(
                systematic)
            for term, variations in self.cut_systematics().items():
                if term == systerm:
                    cuts &= variations[variation]
                else:
                    cuts &= variations['NOMINAL']
        return cuts

    def draw_array_helper(self, field_hist, category, region,
                          cuts=None,
                          weighted=True,
                          field_scale=None,
                          weight_hist=None,
                          field_weight_hist=None,
                          clf=None,
                          scores=None,
                          min_score=None,
                          max_score=None,
                          systematic='NOMINAL',
                          scale=1.,
                          bootstrap_data=False):
        from .data import Data, DataInfo
        all_fields = []
        classifiers = []
        for f in field_hist.iterkeys():
            if isinstance(f, basestring):
                all_fields.append(f)
            elif isinstance(f, Classifier):
                classifiers.append(f)
            elif f is not None:
                all_fields.extend(list(f))
        # hack
        if isinstance(self, Signal) and len(self.modes) == 1 and self.modes[0] == 'gg':
            # for ggH3in systematic
            all_fields.append('true_dphi_jj_higgs_no_overlap')
        if len(classifiers) > 1:
            raise RuntimeError(
                "more than one classifier in fields is not supported")
        elif len(classifiers) == 1:
            classifier = classifiers[0]
        else:
            classifier = None
        if isinstance(self, Data) and bootstrap_data:
            log.info("using bootstrapped data")
            analysis = bootstrap_data
            recs = []
            scores = []
            for s in analysis.backgrounds:
                rec = s.merged_records(category, region,
                    fields=all_fields, cuts=cuts,
                    include_weight=True,
                    clf=clf,
                    systematic=systematic)
                recs.append(rec)
            b_rec = stack(recs, fields=all_fields + ['classifier', 'weight', 'o1','dphi'])
            s_rec = analysis.higgs_125.merged_records(category, region,
                fields=all_fields, cuts=cuts,
                include_weight=True,
                clf=clf,
                systematic=systematic)

            # handle negative weights separately
            b_neg = b_rec[b_rec['weight'] < 0]
            b_pos = b_rec[b_rec['weight'] >= 0]

            def bootstrap(rec):
                prob = np.abs(rec['weight'])
                prob = prob / prob.sum()
                # random sample without replacement
                log.warning(str(int(round(abs(rec['weight'].sum())))))
                sample_idx = np.random.choice(
                    rec.shape[0], size=int(round(abs(rec['weight'].sum()))),
                    replace=False, p=prob)
                return rec[sample_idx]

            rec = stack([
                bootstrap(b_neg),
                bootstrap(b_pos),
                bootstrap(s_rec)],
                fields=all_fields + ['classifier', 'weight'])

            rec['weight'][:] = 1.
            scores = rec['classifier']
        elif all_fields or scores is None:
            # TODO: only get unblinded vars
            rec = self.merged_records(category, region,
                fields=all_fields, cuts=cuts,
                include_weight=True,
                clf=classifier,
                #scores=scores,
                scores=scores[0] if isinstance(scores, tuple) else scores,
                systematic=systematic)
            if scores is None and 'classifier' in rec.dtype.names:
                scores = rec['classifier']
        else:
            # use scores only
            if isinstance(self, Data) and not isinstance(scores, tuple):
                weights = np.ones(len(scores))
            else:
                scores, weights = scores
            rec = np.core.records.fromarrays([scores, weights],
                                             names=['classifier', 'weight'])

        if isinstance(scores, tuple):
            # sanity
            #assert (scores[1] == rec['weight']).all()
            # ignore the score weights since they should be the same as the rec
            # weights
            scores = scores[0]

        weights = rec['weight']

        if min_score is not None:
            if scores is None:
                raise RuntimeError("min_score specified when scores is None")
            # cut below a minimum classifier score
            idx = scores > min_score
            rec = rec[idx]
            weights = weights[idx]
            scores = scores[idx]
        if max_score is not None:
            if scores is None:
                raise RuntimeError("max_score specified when scores is None")
            # cut above a maximum classifier score
            idx = scores < max_score
            rec = rec[idx]
            weights = weights[idx]
            scores = scores[idx]

        def get_weight(array, hist):
            edges = np.array(list(hist.xedges()))
            # handle overflow
            edges[0] -= 1E100
            edges[-1] += 1E100
            return np.array(list(hist.y())).take(
                edges.searchsorted(array) - 1)

        if weight_hist is not None and scores is not None:
            # apply weight according to the classifier score
            log.warning("applying a score weight histogram")
            weights *= get_weight(scores, weight_hist)

        if field_weight_hist is not None:
            # apply weight corrections according to certain fields
            for field, hist in field_weight_hist.items():
                log.warning(
                    "applying a weight histogram with field {0}".format(field))
                if field not in field_hist:
                    raise ValueError(
                        "attempting to apply a weight histogram using "
                        "field {0} but that field is not present in the "
                        "requested array")
                weights *= get_weight(rec[field], hist)

        if scale != 1.:
            weights *= scale

        for fields, hist in field_hist.items():
            if isinstance(fields, Classifier) or fields is None:
                fields = ['classifier']
            # fields can be a single field or list of fields
            elif not isinstance(fields, (list, tuple)):
                fields = [fields]
            if hist is None:
                # this var might be blinded
                continue
            # defensive copy
            if isinstance(fields, tuple):
                # select columns in numpy recarray with a list
                fields = list(fields)
            arr = np.copy(rec[fields])
            if field_scale is not None:
                for field in fields:
                    if field in field_scale:
                        arr[field] *= field_scale[field]
            # convert to array
            arr = rec2array(arr, fields=fields)
            # include the scores if the histogram dimensionality allows
            if scores is not None and hist.GetDimension() == len(fields) + 1:
                arr = np.c_[arr, scores]
            elif hist.GetDimension() != len(fields):
                raise TypeError(
                    'histogram dimensionality does not match '
                    'number of fields: %s' % (', '.join(fields)))
            fill_hist(hist, arr, weights)
            if isinstance(self, Data):
                if hasattr(hist, 'datainfo'):
                    hist.datainfo += self.info
                else:
                    hist.datainfo = DataInfo(self.info.lumi, self.info.energies)

        if scores is not None and 'classifier' not in rec.dtype.names:
            rec = recfunctions.rec_append_fields(rec,
                names='classifier',
                data=scores,
                dtypes='f4')
        return rec, weights

    def events(self, category=None, region=None,
               cuts=None, systematic='NOMINAL', hist=None,
               weighted=True, scale=1.):
        """
        This method returns the number of events selected.  The selection is
        specified by the different arguments.  By default, the output is a
        one-bin histogram with number of event as content.

        Parameters
        ----------
        category :
            A given analysis category. See categories/__init__.py for the list
        region :
            A given analyis regions based on the sign and isolation of the
            taus. The signal region is 'OS'
        cuts :
            In addition to the category (where cuts are specified), extra
            cuts can be added See categories/common.py for a list of possible
            cuts
        systematic :
            By default look at the nominal tree but could also do it
            on specified syst.
        weighted :
            if True, return the weighted number of events
        hist :
            if specified, fill this histogram. if not create a new one an
            return it.
        scale :
            if specified, multiply the number of events by the given
            scale.
        """
        if hist is None:
            hist = Hist(1, -100, 100)
        rec = self.merged_records(category=category, region=region,
                                  cuts=cuts, systematic=systematic)
        if weighted:
            if scale != 1:
                rec['weight'] *= scale
            fill_hist(hist, np.ones(len(rec)), rec['weight'])
        else:
            fill_hist(hist, np.ones(len(rec)))
        return hist

    def events_root(self, category=None, region=None, cuts=None, hist=None,
                    weighted=True, scale=1.):
        """
        QUICK FIX.
        DO NOT USE UNLESS YOU KNOW WHAT YOU'RE DOING
        """
        rfile = get_file(self.ntuple_path, self.student, force_reopen=self.force_reopen)
        if hist is None:
            hist = Hist(1, -100, 100)
        for ds in self.datasets:
            treename = ds.name
            treename = treename.replace('.', '_')
            tree = rfile[treename]
            events = ds.events['NOMINAL']
            weight = LUMI[self.year] * scale * ds.xs * ds.kfact * ds.effic / events
            selection =  (self.cuts(category, region) & cuts)
            if weighted:
                weight_branches = self.weights()
                selection *= Cut(str(weight * scale))
                selection *= Cut('*'.join(weight_branches))
            log.debug("requesing number of events from %s using cuts: %s"
                % (tree.GetName(), selection))
            tree.Draw('1', selection, hist=hist)
        return hist


class Signal(object):
    # mixin
    pass


class Background(object):
    # mixin
    pass


class SystematicsSample(Sample):

    @classmethod
    def get_sys_term_variation(cls, systematic):
        if systematic == 'NOMINAL':
            systerm = None
            variation = 'NOMINAL'
        elif len(systematic) > 1:
            # no support for this yet...
            systerm = None
            variation = 'NOMINAL'
        else:
            systerm, variation = systematic[0].rsplit('_', 1)
        return systerm, variation

    def systematics_components(self):
        common = [
            'MET_RESOSOFTTERMS',
            'MET_SCALESOFTTERMS',
            'TAU_ID',
            'TRIGGER',
        ]
        # No FAKERATE for embedding since fakes are data
        # so don't include FAKERATE here
        if self.year == 2011:
            return common + [
                'TES_TRUE_FINAL',
                'TES_FAKE_FINAL',
            ]
        else:
            return common + [
                'TAU_ID_STAT',
                'TES_TRUE_INSITUINTERPOL',
                'TES_TRUE_SINGLEPARTICLEINTERPOL',
                'TES_TRUE_MODELING',
                'TES_FAKE_TOTAL',
                'TRIGGER_STAT_PERIODA',
                'TRIGGER_STAT_PERIODBD_BARREL',
                'TRIGGER_STAT_PERIODBD_ENDCAP',
                'TRIGGER_STAT_PERIODEM_BARREL',
                'TRIGGER_STAT_PERIODEM_ENDCAP',
            ]

    def weight_systematics(self):
        systematics = {}
        if self.tau_id_sf:
            if self.year == 2011:
                tauid = {
                    'TAU_ID': {
                        'UP': [
                            'tau1_id_sf_high',
                            'tau2_id_sf_high'],
                        'DOWN': [
                            'tau1_id_sf_low',
                            'tau2_id_sf_low'],
                        'NOMINAL': [
                            'tau1_id_sf',
                            'tau2_id_sf']}
                    }
            else:
                tauid = {
                    'TAU_ID': {
                        'STAT_UP': [
                            'tau1_id_sf_stat_high',
                            'tau2_id_sf_stat_high'],
                        'STAT_DOWN': [
                            'tau1_id_sf_stat_low',
                            'tau2_id_sf_stat_low'],
                        'UP': [
                            'tau1_id_sf_sys_high',
                            'tau2_id_sf_sys_high'],
                        'DOWN': [
                            'tau1_id_sf_sys_low',
                            'tau2_id_sf_sys_low'],
                        'NOMINAL': [
                            'tau1_id_sf',
                            'tau2_id_sf']},
                    }
            systematics.update(tauid)
        return systematics

    def cut_systematics(self):
        return {}

    def __init__(self, year, db=DB, systematics=False,
                 tau_id_sf=True, **kwargs):

        if isinstance(self, Background):
            sample_key = self.__class__.__name__.lower()
            sample_info = samples_db.get_sample(
                'hadhad', year, 'background', sample_key)
            kwargs.setdefault('name', sample_info['name'])
            kwargs.setdefault('label', sample_info['root'])
            if 'color' in sample_info and 'color' not in kwargs:
                kwargs['color'] = sample_info['color']
            self.samples = sample_info['samples']

        elif isinstance(self, Signal):
            # samples already defined in Signal subclass
            # see Higgs class below
            assert len(self.samples) > 0

        else:
            raise TypeError(
                'MC sample %s does not inherit from Signal or Background' %
                self.__class__.__name__)

        super(SystematicsSample, self).__init__(year=year, **kwargs)

        self.db = db
        self.datasets = []
        self.systematics = systematics
        self.tau_id_sf = tau_id_sf
        self.norms = {}
        rfile = get_file(self.ntuple_path, self.student, force_reopen=self.force_reopen)
        h5file = get_file(self.ntuple_path, self.student, hdf=True, force_reopen=self.force_reopen)

        from .ztautau import Embedded_Ztautau

        for i, name in enumerate(self.samples):

            ds = self.db[name]
            treename = name.replace('.', '_')
            treename = treename.replace('-', '_')

            tables = {}
            events = {}

            if isinstance(self, Embedded_Ztautau):
                events_bin = 1
            else:
                # use mc_weighted second bin
                events_bin = 2
            events_hist_suffix = '_cutflow'

            tables['NOMINAL'] =  CachedTable.hook(getattr(
                h5file.root, treename))
            cutflow_hist = rfile[treename + events_hist_suffix]
            events['NOMINAL'] = cutflow_hist[events_bin].value
            del cutflow_hist

            if self.systematics:

                systematics_terms, systematics_samples = \
                    samples_db.get_systematics('hadhad', self.year, name)

                if systematics_terms:
                    for sys_term in systematics_terms:
                        sys_name = treename + '_' + '_'.join(sys_term)
                        tables[sys_term] = CachedTable.hook(getattr(
                            h5file.root, sys_name))
                        cutflow_hist = rfile[sys_name + events_hist_suffix]
                        events[sys_term] = cutflow_hist[events_bin].value
                        del cutflow_hist

                if systematics_samples:
                    for sample_name, sys_term in systematics_samples.items():
                        log.info("%s -> %s %s" % (name, sample_name, sys_term))
                        sys_term = tuple(sys_term.split(','))
                        sys_ds = self.db[sample_name]
                        sample_name = sample_name.replace('.', '_')
                        sample_name = sample_name.replace('-', '_')
                        tables[sys_term] = CachedTable.hook(getattr(
                            h5file.root, sample_name))
                        cutflow_hist = rfile[sample_name + events_hist_suffix]
                        events[sys_term] = cutflow_hist[events_bin].value
                        del cutflow_hist

            if hasattr(self, 'xsec_kfact_effic'):
                xs, kfact, effic = self.xsec_kfact_effic(i)
            else:
                xs, kfact, effic = ds.xsec_kfact_effic
            log.debug(
                "dataset: {0}  cross section: {1} [pb] "
                "k-factor: {2} "
                "filtering efficiency: {3} "
                "events {4}".format(
                    ds.name, xs, kfact, effic, events['NOMINAL']))
            dataset = Dataset(ds=ds, tables=tables, events=events,
                              xs=xs, kfact=kfact, effic=effic)
            self.datasets.append(dataset)

    def draw(self, field, hist, category=None, region=None, **kwargs):
        field_scale = {field: get_scale(field)}
        return self.draw_array({field: hist},
                               category=category,
                               region=region,
                               field_scale=field_scale,
                               **kwargs)

    def draw_array(self, field_hist,
                   category=None, region=None,
                   cuts=None,
                   weighted=True,
                   field_scale=None,
                   weight_hist=None,
                   field_weight_hist=None,
                   clf=None,
                   scores=None,
                   min_score=None,
                   max_score=None,
                   systematics=False,
                   systematics_components=None,
                   scale=1.,
                   bootstrap_data=False):

        do_systematics = self.systematics and systematics
        if scores is None and clf is not None:
            scores = self.scores(
                clf, category, region, cuts=cuts,
                systematics=systematics,
                systematics_components=systematics_components)

        rec, weights = self.draw_array_helper(field_hist, category, region,
            cuts=cuts,
            weighted=weighted,
            field_scale=field_scale,
            weight_hist=weight_hist,
            field_weight_hist=field_weight_hist,
            scores=scores['NOMINAL'] if scores else None,
            min_score=min_score,
            max_score=max_score,
            systematic='NOMINAL',
            scale=scale)

        if not do_systematics:
            return rec, weights

        all_sys_hists = {}

        for field, hist in field_hist.items():
            if not hasattr(hist, 'systematics'):
                hist.systematics = {}
            all_sys_hists[field] = hist.systematics

        for systematic in iter_systematics(False,
                year=self.year,
                components=systematics_components):

            sys_field_hist = {}
            for field, hist in field_hist.items():
                if systematic in all_sys_hists[field]:
                    sys_hist = all_sys_hists[field][systematic]
                else:
                    sys_hist = hist.Clone(
                        name=hist.name + '_' + systematic_name(systematic))
                    sys_hist.Reset()
                    all_sys_hists[field][systematic] = sys_hist
                sys_field_hist[field] = sys_hist

            self.draw_array_helper(sys_field_hist, category, region,
                cuts=cuts,
                weighted=weighted,
                field_scale=field_scale,
                weight_hist=weight_hist,
                field_weight_hist=field_weight_hist,
                scores=scores[systematic] if scores else None,
                min_score=min_score,
                max_score=max_score,
                systematic=systematic,
                scale=scale)

        return rec, weights

    def scores(self, clf, category, region,
               cuts=None, scores_dict=None,
               systematics=False,
               systematics_components=None,
               scale=1.):

        # TODO check that weight systematics are included
        do_systematics = self.systematics and systematics
        if scores_dict is None:
            scores_dict = {}
        for systematic in iter_systematics(True,
                year=self.year,
                components=systematics_components):
            if not do_systematics and systematic != 'NOMINAL':
                continue
            scores, weights = clf.classify(self,
                category=category,
                region=region,
                cuts=cuts,
                systematic=systematic)
            weights *= scale
            if systematic not in scores_dict:
                scores_dict[systematic] = (scores, weights)
            else:
                prev_scores, prev_weights = scores_dict[systematic]
                scores_dict[systematic] = (
                    np.concatenate((prev_scores, scores)),
                    np.concatenate((prev_weights, weights)))
        return scores_dict

    def records(self,
                category=None,
                region=None,
                fields=None,
                cuts=None,
                include_weight=True,
                systematic='NOMINAL',
                scale=1.,
                return_idx=False,
                **kwargs):

        from .ztautau import Ztautau
        if include_weight and fields is not None:
            if 'weight' not in fields:
                fields = list(fields) + ['weight']
            if 'o1' not in fields:
                fields = list(fields) + ['o1']
            if 'dphi' not in fields:
                fields = list(fields) + ['dphi']

        selection = self.cuts(category, region, systematic) & cuts
        table_selection = selection.where()
        if systematic == 'NOMINAL':
            log.info("requesting table from %s" %
                     (self.__class__.__name__))
        else:
            log.info("requesting table from %s for systematic %s " %
                     (self.__class__.__name__, systematic_name(systematic)))
        log.debug("using selection: %s" % selection)
        weight_branches = self.weights(systematic)
        if systematic in SYSTEMATICS_BY_WEIGHT:
            systematic = 'NOMINAL'
        recs = []
        if return_idx:
            idxs = []
        for ds in self.datasets:
            try:
                table = ds.tables[systematic]
                events = ds.events[systematic]
            except KeyError:
                log.warning(
                    "table for %s not present for %s "
                    "using NOMINAL" % (systematic, ds.name))
                table = ds.tables['NOMINAL']
                events = ds.events['NOMINAL']
            log.debug(
                "\ndataset: {0}"
                "\ntable: {1}"
                "\ncross section: {2} [pb]"
                "\nk-factor: {3}"
                "\nfiltering efficiency: {4}"
                "\nevents {5}".format(
                    ds.name, table.name, ds.xs, ds.kfact, ds.effic, events))
            actual_scale = self.scale
            if isinstance(self, Ztautau):
                if systematic == ('ZFIT_UP',):
                    log.debug("scaling up for ZFIT_UP")
                    actual_scale += self.scale_error
                elif systematic == ('ZFIT_DOWN',):
                    log.debug("scaling down for ZFIT_DOWN")
                    actual_scale -= self.scale_error
            weight = (
                scale * actual_scale *
                LUMI[self.year] *
                ds.xs * ds.kfact * ds.effic / events)
            if systematic in self.norms:
                weight *= self.norms[systematic]
            # read the table with a selection
            try:
                if table_selection:
                    rec = table.read_where(table_selection, **kwargs)
                else:
                    rec = table.read(**kwargs)
            except Exception as e:
                print table
                print e
                continue
                #raise
            if return_idx:
                # only valid if table_selection is non-empty
                idx = table.get_where_list(table_selection, **kwargs)
                idxs.append(idx)
            # add weight field

            # add optimal observable
            o1 = self.optimate(rec)
            rec = recfunctions.rec_append_fields(
                rec, names='o1',data=o1,dtypes='f8')
            dphi = self.dphate(rec)
            rec = recfunctions.rec_append_fields(
                rec, names='dphi',data=dphi,dtypes='f8')


            if include_weight:
                weights = np.empty(rec.shape[0], dtype='f8')
                weights.fill(weight)
                # merge the weight fields
                weights *= reduce(np.multiply,
                    [rec[br] for br in weight_branches])
                correction_weights = self.corrections(rec)
                if correction_weights:
                    weights *= reduce(np.multiply, correction_weights)
                if isinstance(self,Signal):
                    preHawkScale = np.sum(weights)
                    weights = self.new_optiweight(self.mixing,rec,weights)
                    postHawkScale = np.sum(weights)
                    sf = preHawkScale/postHawkScale
                    weights *= sf
                # add the combined weight
                rec = recfunctions.rec_append_fields(rec,
                    names='weight',
                    data=weights,
                    dtypes='f8')
                if rec['weight'].shape[0] > 1 and rec['weight'].sum() == 0:
                    log.warning("{0}: weights sum to zero!".format(table.name))
            if fields is not None:
                try:
                    rec = rec[fields]
                except Exception as e:
                    print table
                    print rec.shape
                    print rec.dtype
                    print e
                    raise
            recs.append(rec)
        if return_idx:
            return zip(recs, idxs)
        return recs


class MC(SystematicsSample):

    def __init__(self, *args, **kwargs):
        self.pileup_weight = kwargs.pop('pileup_weight', True)
        super(MC, self).__init__(*args, **kwargs)

    def systematics_components(self):
        components = super(MC, self).systematics_components()
        components = components + [
            'JES_Modelling',
            'JES_Detector',
            'JES_EtaModelling',
            'JES_EtaMethod',
            'JES_PURho',
            'JES_FlavComp',
            'JES_FlavResp',
            'JER',
            'FAKERATE',
            'PU_RESCALE',
        ]
        if self.year == 2012:
            components += [
                'JVF',
            ]
        return components

    def weight_fields(self):
        return super(MC, self).weight_fields() + [
            'mc_weight',
            # uncertainty on these are small and are ignored:
            'tau1_fakerate_sf_reco',
            'tau2_fakerate_sf_reco',
        ]

    def weight_systematics(self):
        systematics = super(MC, self).weight_systematics()
        systematics.update({
            'FAKERATE': {
                'UP': [
                    'tau1_fakerate_sf_high',
                    'tau2_fakerate_sf_high'],
                'DOWN': [
                    'tau1_fakerate_sf_low',
                    'tau2_fakerate_sf_low'],
                'NOMINAL': [
                    'tau1_fakerate_sf',
                    'tau2_fakerate_sf']},
            })
        if self.pileup_weight:
            systematics.update({
                'PU_RESCALE': {
                    'UP': ['pileup_weight_high'],
                    'DOWN': ['pileup_weight_low'],
                    'NOMINAL': ['pileup_weight']},
                })
        if self.year == 2011:
            systematics.update({
                'TRIGGER': {
                    'UP': [
                        'tau1_trigger_sf_high',
                        'tau2_trigger_sf_high'],
                    'DOWN': [
                        'tau1_trigger_sf_low',
                        'tau2_trigger_sf_low'],
                    'NOMINAL': [
                        'tau1_trigger_sf',
                        'tau2_trigger_sf']}})
        else:
            systematics.update({
                'TRIGGER': {
                    'UP': [
                        'tau1_trigger_sf_sys_high',
                        'tau2_trigger_sf_sys_high'],
                    'DOWN': [
                        'tau1_trigger_sf_sys_low',
                        'tau2_trigger_sf_sys_low'],
                    'NOMINAL': [
                        'tau1_trigger_sf',
                        'tau2_trigger_sf']},
                'TRIGGER_STAT': {
                    'PERIODA_UP': [
                        'tau1_trigger_sf_stat_scale_PeriodA_high',
                        'tau2_trigger_sf_stat_scale_PeriodA_high'],
                    'PERIODA_DOWN': [
                        'tau1_trigger_sf_stat_scale_PeriodA_low',
                        'tau2_trigger_sf_stat_scale_PeriodA_low'],
                    'PERIODBD_BARREL_UP': [
                        'tau1_trigger_sf_stat_scale_PeriodBD_Barrel_high',
                        'tau2_trigger_sf_stat_scale_PeriodBD_Barrel_high'],
                    'PERIODBD_BARREL_DOWN': [
                        'tau1_trigger_sf_stat_scale_PeriodBD_Barrel_low',
                        'tau2_trigger_sf_stat_scale_PeriodBD_Barrel_low'],
                    'PERIODBD_ENDCAP_UP': [
                        'tau1_trigger_sf_stat_scale_PeriodBD_EndCap_high',
                        'tau2_trigger_sf_stat_scale_PeriodBD_EndCap_high'],
                    'PERIODBD_ENDCAP_DOWN': [
                        'tau1_trigger_sf_stat_scale_PeriodBD_EndCap_low',
                        'tau2_trigger_sf_stat_scale_PeriodBD_EndCap_low'],
                    'PERIODEM_BARREL_UP': [
                        'tau1_trigger_sf_stat_scale_PeriodEM_Barrel_high',
                        'tau2_trigger_sf_stat_scale_PeriodEM_Barrel_high'],
                    'PERIODEM_BARREL_DOWN': [
                        'tau1_trigger_sf_stat_scale_PeriodEM_Barrel_low',
                        'tau2_trigger_sf_stat_scale_PeriodEM_Barrel_low'],
                    'PERIODEM_ENDCAP_UP': [
                        'tau1_trigger_sf_stat_scale_PeriodEM_EndCap_high',
                        'tau2_trigger_sf_stat_scale_PeriodEM_EndCap_high'],
                    'PERIODEM_ENDCAP_DOWN': [
                        'tau1_trigger_sf_stat_scale_PeriodEM_EndCap_low',
                        'tau2_trigger_sf_stat_scale_PeriodEM_EndCap_low'],
                    'NOMINAL': []}})
        return systematics


class CompositeSample(object):
    """
    This class adds together the events from a list of samples
    and also return the summed histograms of all of those samples
    for the requested fields
    TODO: Implement a naming from the components.
    """
    def __init__(self, samples_list, name='Sample', label='Sample'):
        if not isinstance( samples_list, (list,tuple)):
            samples_list = [samples_list]
        if not isinstance (samples_list[0], Sample):
            raise ValueError( "samples_list must be filled with Samples")
        self.samples_list = samples_list
        self.name = name
        self.label = label

    def events(self, *args, **kwargs ):
        """
        Return a one-bin histogram with the total sum of events
        of all the samples
        Parameters:
        - See the events() method in the Sample class
        """
        return sum([s.events(*args, **kwargs) for s in self.samples_list])

    def draw_array(self, field_hist_tot, category, region,
                   systematics=False, **kwargs):
        """
        Construct histograms of the sum of all the samples.
        Parameters:
        - field_hist_tot: dictionnary of Histograms that constain the structure we want to retrieve
        - category: the analysis category
        - region: the analysis region (for example 'OS')
        - systematics: boolean flag
        """
        field_hists_list = []
        # -------- Retrieve the histograms dictionnary from each sample and store it into a list
        for s in self.samples_list:
            # field_hists_temp = s.get_hist_array( field_hist_tot, category, region, systematics=systematics,**kwargs)
            field_hists_temp = {}
            for field,hist in field_hist_tot.items():
                field_hists_temp[field] = hist.Clone()
                field_hists_temp[field].Reset()
            s.draw_array(field_hists_temp, category, region, systematics=systematics,**kwargs)
            field_hists_list.append( field_hists_temp )

        # -------- Reset the output histograms
        for field, hist in field_hists_list[0].items():
            hist_tot = hist.Clone()
            hist_tot.Reset()
            field_hist_tot[field] = hist_tot

        # -------- Add the nominal histograms
        for field_hist in field_hists_list:
            for field, hist in field_hist.items():
                field_hist_tot[field].Add( hist )

        # --- Systematic Uncertainties block
        if systematics:
            #--- loop over the dictionnary of the summed histograms
            for field,hist in field_hist_tot.items():
                # --- Add a dictionary to the nominal summed histogram
                if not hasattr( hist,'systematics'):
                    hist.systematics = {}
                # --- loop over the systematic uncercainties
                for sys in iter_systematics(self.samples_list[0].year):
                    if sys is 'NOMINAL':
                        continue
                    log.info ( "Fill the %s syst for the field %s" % (sys,field) )
                    # -- Create an histogram for each systematic uncertainty
                    hist.systematics[sys] =  hist.Clone()
                    hist.systematics[sys].Reset()
                    # -- loop over the samples and sum-up the syst-applied histograms
                    for field_hist_sample in field_hists_list:
                        field_hist_syst = field_hist_sample[field].systematics
                        hist.systematics[sys].Add( field_hist_syst[sys] )
        return
