# rootpy imports
from rootpy import asrootpy
from rootpy.plotting import Hist

# numpy imports
import numpy as np

# local imports
from .sample import Sample, Background
from . import log; log = log[__name__]
from ..systematics import systematic_name
from ..regions import REGION_SYSTEMATICS
from ..defaults import FAKES_REGION


class QCD(Sample, Background):
    # don't include MC systematics in workspace for QCD
    WORKSPACE_SYSTEMATICS = [] #MC.WORKSPACE_SYSTEMATICS
    NORM_BY_THEORY = False

    def systematics_components(self):
        return []

    def histfactory(self, sample, category, systematics=True, **kwargs):
        if self.workspace_norm is False:
            return
        elif self.workspace_norm is not None:
            sample.AddNormFactor(
                'ATLAS_norm_HH_{0:d}_QCD'.format(self.year),
                self.workspace_norm,
                self.workspace_norm,
                self.workspace_norm,
                True) # const
        elif self.constrain_norm:
            # overallsys
            error = self.scale_error / self.scale
            sample.AddOverallSys(
                'ATLAS_norm_HH_{0:d}_QCD'.format(self.year),
                1. - error,
                1. + error)
        else:
            sample.AddNormFactor(
                'ATLAS_norm_HH_{0:d}_QCD'.format(self.year),
                1., 0., 50., False) # floating

    @staticmethod
    def sample_compatibility(data, mc):
        if not isinstance(mc, (list, tuple)):
            raise TypeError("mc must be a list or tuple of MC samples")
        if not mc:
            raise ValueError("mc must contain at least one MC sample")
        systematics = mc[0].systematics
        for m in mc:
            if data.year != m.year:
                raise ValueError("MC and Data years do not match")
            if m.systematics != systematics:
                raise ValueError(
                    "two MC samples with inconsistent systematics setting")

    def __init__(self, data, mc,
                 scale=1.,
                 scale_error=0.,
                 data_scale=1.,
                 mc_scales=None,
                 shape_region=FAKES_REGION,
                 decouple_shape=False,
                 coherent_shape=True,
                 workspace_norm=None,
                 constrain_norm=False,
                 shape_systematic=True,
                 cuts=None,
                 name='Fakes',
                 label='Fakes',
                 **kwargs):
        QCD.sample_compatibility(data, mc)
        super(QCD, self).__init__(
            year=data.year,
            scale=scale,
            name=name,
            label=label,
            **kwargs)
        self.data = data
        self.mc = mc
        self.scale = 1.
        self.data_scale = data_scale
        if mc_scales is not None:
            if len(mc_scales) != len(mc):
                raise ValueError("length of MC scales must match number of MC")
            self.mc_scales = mc_scales
        else:
            # default scales to 1.
            self.mc_scales = [1. for m in self.mc]
        self.scale_error = scale_error
        self.shape_region = shape_region
        self.decouple_shape = decouple_shape
        self.coherent_shape = coherent_shape
        self.workspace_norm = workspace_norm
        self.constrain_norm = constrain_norm
        self.shape_systematic = shape_systematic
        self.systematics = mc[0].systematics

    def draw_array(self, field_hist, category, region,
                   cuts=None,
                   weighted=True,
                   field_scale=None,
                   weight_hist=None,
                   field_weight_hist=None,
                   clf=None,
                   scores=None,
                   min_score=None,
                   max_score=None,
                   systematics=True,
                   systematics_components=None,
                   bootstrap_data=False):

        if scores is not None:
            log.warning(
                "scores is not None in QCD.draw_array() and will be ignored")
            scores = None

        # TODO: support for field_weight_hist
        do_systematics = self.systematics and systematics

        field_hist_MC_bkg = dict([(expr, hist.Clone())
            for expr, hist in field_hist.items()])

        for mc_scale, mc in zip(self.mc_scales, self.mc):
            mc.draw_array(field_hist_MC_bkg, category, self.shape_region,
                cuts=cuts,
                weighted=weighted,
                field_scale=field_scale,
                weight_hist=weight_hist,
                clf=clf,
                #scores=scores,
                min_score=min_score,
                max_score=max_score,
                systematics=systematics,
                systematics_components=systematics_components,
                scale=mc_scale)

        field_hist_data = dict([(expr, hist.Clone())
            for expr, hist in field_hist.items()])

        self.data.draw_array(field_hist_data,
            category, self.shape_region,
            cuts=cuts,
            weighted=weighted,
            field_scale=field_scale,
            weight_hist=weight_hist,
            clf=clf,
            #scores=scores,
            min_score=min_score,
            max_score=max_score)

        for expr, h in field_hist.items():
            mc_h = field_hist_MC_bkg[expr]
            d_h = field_hist_data[expr]
            h += (d_h * self.data_scale - mc_h) * self.scale
            h.SetTitle(self.label)

        if do_systematics and systematics_components is None:
            # get shape systematic
            field_shape_sys = self.get_shape_systematic_array(
                # use the nominal hist here
                field_hist,
                category, region,
                cuts=cuts,
                clf=clf,
                min_score=min_score,
                max_score=max_score,
                field_scale=field_scale,
                weight_hist=weight_hist)
        else:
            field_shape_sys = None

        for expr, h in field_hist.items():
            mc_h = field_hist_MC_bkg[expr]
            d_h = field_hist_data[expr]

            if not do_systematics or not hasattr(mc_h, 'systematics'):
                continue
            if not hasattr(h, 'systematics'):
                h.systematics = {}

            if field_shape_sys is not None:
                # add shape systematics
                high, low = field_shape_sys[expr]
                h.systematics[('QCDSHAPE_DOWN',)] = low
                h.systematics[('QCDSHAPE_UP',)] = high

            for sys_term, sys_hist in mc_h.systematics.items():
                if sys_term in (('QCDSHAPE_UP',), ('QCDSHAPE_DOWN',)):
                    continue
                scale = self.scale
                if sys_term == ('QCDFIT_UP',):
                    scale = self.scale + self.scale_error
                elif sys_term == ('QCDFIT_DOWN',):
                    scale = self.scale - self.scale_error
                qcd_hist = (d_h * self.data_scale - sys_hist) * scale
                qcd_hist.name = h.name + '_' + systematic_name(sys_term)
                if sys_term not in h.systematics:
                    h.systematics[sys_term] = qcd_hist
                else:
                    h.systematics[sys_term] += qcd_hist
        # hack: no rec or weights
        return None, None

    def scores(self, clf, category, region,
               cuts=None,
               systematics=True,
               systematics_components=None,
               **kwargs):
        # data
        data_scores, data_weights = self.data.scores(
            clf,
            category,
            region=self.shape_region,
            cuts=cuts,
            **kwargs)

        scores_dict = {}
        # subtract MC
        for mc_scale, mc in zip(self.mc_scales, self.mc):
            mc.scores(
                clf,
                category,
                region=self.shape_region,
                cuts=cuts,
                scores_dict=scores_dict,
                systematics=systematics,
                systematics_components=systematics_components,
                scale=mc_scale,
                **kwargs)

        for sys_term in scores_dict.keys()[:]:
            sys_scores, sys_weights = scores_dict[sys_term]
            scale = self.scale
            if sys_term == ('QCDFIT_UP',):
                scale += self.scale_error
            elif sys_term == ('QCDFIT_DOWN',):
                scale -= self.scale_error
            # subtract MC
            sys_weights *= -1 * scale
            # add data
            # same order as in records()
            sys_scores = np.concatenate(
                (np.copy(data_scores), sys_scores))
            sys_weights = np.concatenate(
                (data_weights * scale, sys_weights))
            scores_dict[sys_term] = (sys_scores, sys_weights)

        return scores_dict

    def records(self,
                category=None,
                region=None,
                fields=None,
                cuts=None,
                include_weight=True,
                systematic='NOMINAL',
                return_idx=False,
                **kwargs):
        assert include_weight == True
        data_records = self.data.records(
            category=category,
            region=self.shape_region,
            fields=fields,
            cuts=cuts,
            include_weight=include_weight,
            systematic='NOMINAL',
            return_idx=return_idx,
            **kwargs)

        if return_idx:
            arrays = [(d.copy(), idx) for d, idx in data_records]
        else:
            arrays = [d.copy() for d in data_records]

        for mc_scale, mc in zip(self.mc_scales, self.mc):
            _arrays = []
            _arrs = mc.records(
                category=category,
                region=self.shape_region,
                fields=fields,
                cuts=cuts,
                include_weight=include_weight,
                systematic=systematic,
                scale=mc_scale,
                return_idx=return_idx,
                **kwargs)
            # FIX: weight may not be present if include_weight=False
            if return_idx:
                for partition, idx in _arrs:
                    partition = partition.copy()
                    partition['weight'] *= -1
                    _arrays.append((partition, idx))
            else:
                for partition in _arrs:
                    partition = partition.copy()
                    partition['weight'] *= -1
                    _arrays.append(partition)
            arrays.extend(_arrays)

        scale = self.scale
        if systematic == ('QCDFIT_UP',):
            scale += self.scale_error
        elif systematic == ('QCDFIT_DOWN',):
            scale -= self.scale_error

        # FIX: weight may not be present if include_weight=False
        if return_idx:
            for partition, idx in arrays:
                partition['weight'] *= scale
        else:
            for partition in arrays:
                partition['weight'] *= scale

        return arrays

    def get_shape_systematic(self, nominal_hist, expr_or_clf,
                             category, region, **kwargs):
        return self.get_shape_systematic_array(
            {expr_or_clf: nominal_hist}, category, region,
            **kwargs)[expr_or_clf]

    def get_shape_systematic_array(self, field_nominal_hist,
                                   category, region,
                                   cuts=None,
                                   clf=None,
                                   min_score=None,
                                   max_score=None,
                                   suffix=None,
                                   field_scale=None,
                                   weight_hist=None,
                                   weighted=True):

        log.info("creating QCD shape systematic")

        curr_model = self.shape_region
        if curr_model not in REGION_SYSTEMATICS:
            raise ValueError(
                "no QCD shape systematic defined for nominal {0}".format(
                    curr_model))
        shape_models = REGION_SYSTEMATICS[curr_model]
        if not isinstance(shape_models, (tuple, list)):
            shape_models = [shape_models]
        elif len(shape_models) > 2:
            raise ValueError("a maximum of two shape models is supported")
        # reflect single shape about the nominal
        reflect = len(shape_models) == 1

        # use preselection as reference in which all models should have the
        # same expected number of QCD events
        # get number of events at preselection for nominal model
        from ..categories import Category_Preselection
        nominal_events = self.events(Category_Preselection, None)[1].value

        field_hist_template = {}
        for field, hist in field_nominal_hist.items():
            new_hist = hist.Clone()
            new_hist.Reset()
            field_hist_template[field] = new_hist

        # get QCD shape systematic
        model_events = {}
        model_field_shape_sys = {}
        for model in shape_models:
            self.shape_region = model
            log.info("getting QCD shape for {0}".format(model))
            model_field_shape_sys[model] = self.get_hist_array(
                field_hist_template,
                category, region,
                cuts=cuts,
                clf=clf,
                scores=None,
                min_score=min_score,
                max_score=max_score,
                systematics=False,
                suffix=(suffix or '') + '_{0}'.format(model),
                field_scale=field_scale,
                weight_hist=weight_hist,
                weighted=weighted)[0]
            model_events[model] = self.events(
                Category_Preselection, None)[1].value
        # restore previous shape model
        self.shape_region = curr_model

        def fix_empty_shape(nominal, shape):
            """
            If the shape uncertainty has zero/negative events in a given bin
            where the nominal model has events, then assign at least 100%
            uncertainty in the same direction as the variation in the previous
            bin.
            """
            for bin_nom, bin_shape in zip(nominal.bins(), shape.bins()):
                if bin_shape.value <= 0 and bin_nom.value > 0:
                    idx = bin_shape.idx
                    if shape[idx - 1].value - nominal[idx - 1].value >= 0:
                        # at least 100% uncertainty
                        bin_shape.value = 2 * bin_nom.value - bin_shape.value
                    # fill with average of variation to the left and right
                    #bin_shape.value = bin_nom.value + (
                    #    (shape[idx - 1].value - nominal[idx - 1].value) +
                    #    (shape[idx + 1].value - nominal[idx + 1].value)) / 2.

        if reflect:
            # single shape
            model = shape_models[0]
            field_shape_sys_reflect = {}
            for field, shape_sys in model_field_shape_sys[model].items():
                nominal_hist = field_nominal_hist[field]
                # normalize shape_sys such that it would have the same number of
                # events as the nominal at preselection
                shape_sys *= nominal_events / float(model_events[model])
                fix_empty_shape(nominal_hist, shape_sys)
                # reflect shape about the nominal to get high and low variations
                shape_sys_reflect = nominal_hist + (nominal_hist - shape_sys)
                shape_sys_reflect.name = shape_sys.name + '_reflected'
                field_shape_sys_reflect[field] = (shape_sys, shape_sys_reflect)

            return field_shape_sys_reflect

        model_high, model_low = shape_models
        field_shape_sys = {}
        for field, shape_sys_low in model_field_shape_sys[model_low].items():
            shape_sys_high = model_field_shape_sys[model_high][field]
            nominal_hist = field_nominal_hist[field]
            # normalize shape_sys such that it would have the same number of
            # events as the nominal at preselection
            shape_sys_low *= nominal_events / float(model_events[model_low])
            shape_sys_high *= nominal_events / float(model_events[model_high])
            fix_empty_shape(nominal_hist, shape_sys_low)
            fix_empty_shape(nominal_hist, shape_sys_high)
            field_shape_sys[field] = (shape_sys_high, shape_sys_low)

        return field_shape_sys
