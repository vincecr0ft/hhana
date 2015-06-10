from . import log
from .draw import draw
from rootpy.stats.histfactory import HistoSys, split_norm_shape


def draw_channel(channel, fit=None, no_data=False,
                 ypadding=None, log_ypadding=None, **kwargs):
    """
    Draw a HistFactory::Channel only include OverallSys systematics
    in resulting band as an illustration of the level of uncertainty
    since correlations of the NPs are not known and it is not
    possible to draw the statistically correct error band.
    """
    if fit is not None:
        log.warning("applying snapshot on channel {0}".format(channel.name))
        channel = channel.apply_snapshot(fit)
    if channel.data and channel.data.hist and not no_data:
        data_hist = channel.data.hist
    else:
        data_hist = None
    model_hists = []
    signal_hists = []
    systematics_terms = {}
    for sample in channel.samples:
        nominal_hist = sample.hist.Clone(shallow=True)
        _systematics = {}
        for sys_name, osys, hsys in sample.iter_sys():
            systematics_terms[sys_name] = (
                sys_name + '_UP',
                sys_name + '_DOWN')
            if hsys is not None:
                # include only overallsys component
                norm, shape = split_norm_shape(hsys, nominal_hist)
                if osys is not None:
                    osys.low *= norm.low
                    osys.high *= norm.high
                else:
                    osys = norm
            _systematics[sys_name + '_DOWN'] = nominal_hist * osys.low
            _systematics[sys_name + '_UP'] = nominal_hist * osys.high
            log.debug("sample: {0} overallsys: {1} high: {2} low: {3}".format(
                sample.name, sys_name, osys.high, osys.low))
        nominal_hist.systematics = _systematics
        if sample.GetNormFactor('SigXsecOverSM') is not None:
            XsecNF=5.8777/nominal_hist.Integral() #SM integral                                          
            print 'scaling by ',XsecNF
            nominal_hist.Scale(XsecNF)
            signal_hists.append(nominal_hist)
        else:
            model_hists.append(nominal_hist)
    if 'systematics' in kwargs:
        del kwargs['systematics']
    figs = []
    for logy in (False, True):
        figs.append(draw(
            data=data_hist,
            model=model_hists or None,
            signal=signal_hists or None,
            systematics=systematics_terms,
            logy=logy,
            ypadding=(log_ypadding or ypadding) if logy else ypadding,
            **kwargs))
    return figs


def draw_channel_array(
        analysis,
        vars,
        category,
        region,
        cuts=None,
        mass=125,
        mode=None,
        scale_125=False,
        ravel=False,
        weighted=True,
        weight_hist=None,
        clf=None,
        min_score=None,
        max_score=None,
        templates=None,
        plots=None,
        output_suffix='',
        unblind=False,
        bootstrap_data=False,
        **kwargs):
    # filter out plots that will not be made
    used_vars = {}
    field_scale = {}
    if plots is not None:
        for plot in plots:
            if plot in vars:
                var_info = vars[plot]
                if (var_info.get('cats', None) is not None and
                    category.name.upper() not in var_info['cats']):
                    raise ValueError(
                        "variable %s is not valid in the category %s" %
                        (plot, category.name.upper()))
                used_vars[plot] = var_info
                if 'scale' in var_info:
                    field_scale[plot] = var_info['scale']
            else:
                raise ValueError(
                    "variable %s is not defined in mva/variables.py" % plot)
    else:
        for expr, var_info in vars.items():
            if (var_info.get('cats', None) is not None and
                category.name.upper() not in var_info['cats']):
                continue
            used_vars[expr] = var_info
            if 'scale' in var_info:
                field_scale[expr] = var_info['scale']
    vars = used_vars
    if not vars:
        raise RuntimeError("no variables selected")

    field_channel = analysis.get_channel_array(vars,
        category, region, cuts,
        include_signal=True,
        mass=mass,
        mode=mode,
        scale_125=scale_125,
        clf=clf,
        min_score=min_score,
        max_score=max_score,
        weighted=weighted,
        templates=templates,
        field_scale=field_scale,
        weight_hist=weight_hist,
        no_signal_fixes=True,
        bootstrap_data=bootstrap_data)

    figs = {}
    for field, var_info in vars.items():
        if unblind:
            blind = False
        else:
            blind = var_info.get('blind', False)
        output_name = var_info['filename'] + output_suffix
        if cuts:
            output_name += '_' + cuts.safe()
        ypadding = kwargs.pop('ypadding', var_info.get('ypadding', None))
        log_ypadding = kwargs.pop('log_ypadding', var_info.get('log_ypadding', None))
        legend_position = kwargs.pop('legend_position', var_info.get('legend', 'right'))
        fig = draw_channel(field_channel[field],
                           data_info=str(analysis.data.info),
                           category=category,
                           name=var_info['root'],
                           units=var_info.get('units', None),
                           output_name=output_name,
                           blind=blind,
                           integer=var_info.get('integer', False),
                           ypadding=ypadding,
                           log_ypadding=log_ypadding,
                           legend_position=legend_position,
                           **kwargs)
        figs[field] = fig
    return field_channel, figs
