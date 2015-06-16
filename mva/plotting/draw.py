# stdlib imports
import math

# ROOT/rootpy imports
import ROOT
from rootpy.plotting import Legend, Hist, HistStack
from rootpy.plotting.shapes import Arrow
import rootpy.plotting.utils as rootpy_utils

# local imports
from .. import PLOTS_DIR, save_canvas
from .templates import RatioPlot, SimplePlot
from ..utils import fold_overflow
from .utils import label_plot, legend_params, set_colors
from . import log


def uncertainty_band(model, systematics): #, systematics_components):
    # TODO determine systematics from model itself
    if not isinstance(model, (list, tuple)):
        model = [model]
    # add separate variations in quadrature
    # also include stat error in quadrature
    total_model = sum(model)
    var_high = []
    var_low = []
    for term, variations in systematics.items():
        if len(variations) == 2:
            high, low = variations
        elif len(variations) == 1:
            high = variations[0]
            low = 'NOMINAL'
        else:
            raise ValueError(
                "only one or two variations "
                "per term are allowed: {0}".format(variations))

        """
        if systematics_components is not None:
            if high not in systematics_components:
                log.warning("filtering out {0}".format(high))
                high = 'NOMINAL'
            if low not in systematics_components:
                log.warning("filtering out {0}".format(low))
                low = 'NOMINAL'
        """

        if high == 'NOMINAL' and low == 'NOMINAL':
            continue

        total_high = model[0].Clone()
        total_high.Reset()
        total_low = total_high.Clone()
        total_max = total_high.Clone()
        total_min = total_high.Clone()
        for m in model:
            if high == 'NOMINAL' or high not in m.systematics:
                total_high += m.Clone()
            else:
                #print m.title, high, list(m.systematics[high])
                total_high += m.systematics[high]

            if low == 'NOMINAL' or low not in m.systematics:
                total_low += m.Clone()
            else:
                #print m.title, low, list(m.systematics[low])
                total_low += m.systematics[low]

        if total_low.Integral() <= 0:
            log.warning("{0}_DOWN is non-positive".format(term))
        if total_high.Integral() <= 0:
            log.warning("{0}_UP is non-positive".format(term))

        for i in total_high.bins_range(overflow=True):
            total_max[i].value = max(total_high[i].value, total_low[i].value, total_model[i].value)
            total_min[i].value = min(total_high[i].value, total_low[i].value, total_model[i].value)

        if total_min.Integral() <= 0:
            log.warning("{0}: lower bound is non-positive".format(term))
        if total_max.Integral() <= 0:
            log.warning("{0}: upper bound is non-positive".format(term))

        var_high.append(total_max)
        var_low.append(total_min)

        log.debug("{0} {1}".format(str(term), str(variations)))
        log.debug("{0} {1} {2}".format(
            total_max.integral(),
            total_model.integral(),
            total_min.integral()))

    #log.debug(str(systematics_components))
    # include stat error variation
    total_model_stat_high = total_model.Clone()
    total_model_stat_low = total_model.Clone()
    for i in total_model.bins_range(overflow=True):
        total_model_stat_high[i].value += total_model.yerrh(i)
        total_model_stat_low[i].value -= total_model.yerrl(i)
    var_high.append(total_model_stat_high)
    var_low.append(total_model_stat_low)

    # sum variations in quadrature bin-by-bin
    high_band = total_model.Clone()
    high_band.Reset()
    low_band = high_band.Clone()
    for i in high_band.bins_range(overflow=True):
        sum_high = math.sqrt(
            sum([(v[i].value - total_model[i].value)**2 for v in var_high]))
        sum_low = math.sqrt(
            sum([(v[i].value - total_model[i].value)**2 for v in var_low]))
        high_band[i].value = sum_high
        low_band[i].value = sum_low
    return total_model, high_band, low_band


def draw(name,
         category,
         data=None,
         data_info=None,
         model=None,
         model_colors=None,
         signal=None,
         signal_scale=1.,
         signal_on_top=False,
         signal_linestyles=None,
         signal_colors=None,
         show_signal_error=False,
         fill_signal=False,
         stack_signal=False,
         units=None,
         plot_label=None,
         ylabel='Events',
         blind=False,
         show_ratio=False,
         ratio_range=(0, 2),
         ratio_height=0.15,
         ratio_margin=0.06,
         output_formats=None,
         systematics=None,
         #systematics_components=None,
         integer=False,
         textsize=22,
         logy=False,
         logy_min=None,
         separate_legends=False,
         ypadding=None,
         legend_position='right',
         range=None,
         output_name=None,
         output_dir=PLOTS_DIR,
         arrow_values=None,
         overflow=False,
         show_pvalue=False,
         top_label=None,
         poisson_errors=True):

    if model is None and data is None and signal is None:
        # insufficient input
        raise ValueError(
            "at least one of model, data, "
            "or signal must be specified")

    if model is not None:
        if not isinstance(model, (list, tuple)):
            model = [model]
        if overflow:
            for hist in model:
                fold_overflow(hist)
    if signal is not None:
        if not isinstance(signal, (list, tuple)):
            signal = [signal]
        if overflow:
            for hist in signal:
                fold_overflow(hist)
    if data is not None and overflow:
        fold_overflow(data)

    # objects will be populated with all histograms in the main pad
    objects = []
    legends = []

    if show_ratio and (data is None or model is None):
        # cannot show the ratio if data or model was not specified
        show_ratio=False

    if ypadding is None:
        # select good defaults for log or linear scales
        if logy:
            ypadding = (.6, .0)
        else:
            ypadding = (.35, .0)

    template = data or model[0]
    xdivisions = min(template.nbins(), 7) if integer else 507

    if show_ratio:
        fig = RatioPlot(
            logy=logy,
            ratio_title='Data / Model',
            ratio_limits=(0, 2),
            offset=-72,
            ratio_margin=22,
            prune_ratio_ticks=True)
    else:
        fig = SimplePlot(logy=logy)

    if signal is not None:
        if signal_scale != 1.:
            scaled_signal = []
            for sig in signal:
                scaled_h = sig * signal_scale
                scaled_h.SetTitle(r'%g #times %s' % (
                    signal_scale,
                    sig.GetTitle()))
                scaled_signal.append(scaled_h)
        else:
            scaled_signal = signal
        if signal_colors is not None:
            set_colors(scaled_signal, signal_colors)
        for i, s in enumerate(scaled_signal):
            s.drawstyle = 'HIST'
            if fill_signal:
                s.fillstyle = 'solid'
                s.fillcolor = s.linecolor
                s.linewidth = 0
                s.linestyle = 'solid'
                alpha = .75
            else:
                s.fillstyle = 'hollow'
                s.linewidth = 3
                if signal_linestyles is not None:
                    s.linestyle = signal_linestyles[i]
                else:
                    s.linestyle = 'solid'
                alpha = 1.

    if model is not None:
        if model_colors is not None:
            set_colors(model, model_colors)
        # create the model stack
        model_stack = HistStack()
        for hist in model:
            hist.SetLineWidth(0)
            hist.drawstyle = 'hist'
            model_stack.Add(hist)
        if signal is not None and signal_on_top:
            if stack_signal:
                for s in scaled_signal:
                    model_stack.Add(s)                
        objects.append(model_stack)

    if model is not None:
        if model_colors is not None:
            set_colors(model, model_colors)
        model_stack = HistStack()
        for hist in model:
            hist.SetLineWidth(0)
            hist.drawstyle = 'hist'
            model_stack.Add(hist)
        if signal is not None and signal_on_top:
            if not stack_signal:
                for s in scaled_signal:
                    this_stack=model_stack
                    this_stack.Add(s)                
                    objects.extend(this_stack)


    if signal is not None and not signal_on_top:
        if stack_signal:
            # create the signal stack
            signal_stack = HistStack()
            for hist in scaled_signal:
                signal_stack.Add(hist)
            objects.append(signal_stack)
        else:
            objects.extend(scaled_signal)

    if model is not None:
        # draw uncertainty band
        total_model, high_band_model, low_band_model = uncertainty_band(
            model, systematics) #, systematics_components)
        high = total_model + high_band_model
        low = total_model - low_band_model
        error_band_model = rootpy_utils.get_band(
            low, high,
            middle_hist=total_model)
        error_band_model.fillstyle = '/'
        error_band_model.fillcolor = 13
        error_band_model.linecolor = 10
        error_band_model.markersize = 0
        error_band_model.markercolor = 10
        error_band_model.drawstyle = 'e2'
        objects.append(error_band_model)

    if signal is not None and show_signal_error:
        total_signal, high_band_signal, low_band_signal = uncertainty_band(
            signal, systematics) #, systematics_components)
        high = (total_signal + high_band_signal) * signal_scale
        low = (total_signal - low_band_signal) * signal_scale
        if signal_on_top:
            high += total_model
            low += total_model
        error_band_signal = rootpy_utils.get_band(
            low, high,
            middle_hist=total_signal * signal_scale)
        error_band_signal.fillstyle = '\\'
        error_band_signal.fillcolor = 13
        error_band_signal.linecolor = 10
        error_band_signal.markersize = 0
        error_band_signal.markercolor = 10
        error_band_signal.drawstyle = 'e2'
        objects.append(error_band_signal)

    if data is not None and blind is not True:
        # create the data histogram
        if isinstance(blind, tuple):
            low, high = blind
            # zero out bins in blind region
            for bin in data.bins():
                if (low < bin.x.high <= high or low <= bin.x.low < high):
                    data[bin.idx] = (0., 0.)

        if poisson_errors:
            # convert data to TGraphAsymmErrors with Poisson errors
            data_poisson = data.poisson_errors()
            data_poisson.markersize = 1.2
            data_poisson.drawstyle = 'PZ'
            objects.append(data_poisson)
        else:
            # Gaussian errors
            data.drawstyle = 'E0'
            objects.append(data)

        # draw ratio plot
        if model is not None and show_ratio:
            fig.cd('ratio')
            total_model = sum(model)
            ratio_hist = Hist.divide(data, total_model)
            # remove bins where data is zero
            max_dev = 0
            for bin in data.bins():
                if bin.value <= 0:
                    ratio_hist[bin.idx] = (-100, 0)
                else:
                    ratio_value = ratio_hist[bin.idx].value
                    dev =  abs(ratio_value - 1)
                    if dev > max_dev:
                        max_dev = dev

            if max_dev < 0.2:
                ratio_range = (0.8, 1.2)
            elif max_dev < 0.4:
                ratio_range = (0.6, 1.4)
            ruler_high = (ratio_range[1] + 1.) / 2.
            ruler_low = (ratio_range[0] + 1.) / 2.

            ratio_hist.linecolor = 'black'
            ratio_hist.linewidth = 2
            ratio_hist.fillstyle = 'hollow'
            ratio_hist.drawstyle = 'E0'

            """
            # draw empty copy of ratio_hist first so lines will show
            ratio_hist_tmp = ratio_hist.Clone()
            ratio_hist_tmp.Reset()
            ratio_hist_tmp.Draw()
            ratio_hist_tmp.yaxis.SetLimits(*ratio_range)
            ratio_hist_tmp.yaxis.SetRangeUser(*ratio_range)
            ratio_hist_tmp.yaxis.SetTitle('Data / Model')
            ratio_hist_tmp.yaxis.SetNdivisions(4)
            # not certain why the following is needed
            ratio_hist_tmp.yaxis.SetTitleOffset(style.GetTitleYOffset())

            ratio_xrange = range or ratio_hist.bounds()

            ratio_hist_tmp.xaxis.SetLimits(*ratio_xrange)
            #ratio_hist_tmp.xaxis.SetRangeUser(*ratio_xrange)

            ratio_hist_tmp.xaxis.SetTickLength(
                ratio_hist_tmp.xaxis.GetTickLength() * 2)

            # draw ratio=1 line
            line = Line(ratio_xrange[0], 1,
                        ratio_xrange[1], 1)
            line.linestyle = 'dashed'
            line.linewidth = 2
            line.Draw()

            # draw high ratio line
            line_up = Line(ratio_xrange[0], ruler_high,
                           ratio_xrange[1], ruler_high)
            line_up.linestyle = 'dashed'
            line_up.linewidth = 2
            line_up.Draw()

            # draw low ratio line
            line_dn = Line(ratio_xrange[0], ruler_low,
                           ratio_xrange[1], ruler_low)
            line_dn.linestyle = 'dashed'
            line_dn.linewidth = 2
            line_dn.Draw()
            """

            # draw band below points on ratio plot
            ratio_hist_high = Hist.divide(
                total_model + high_band_model, total_model)
            ratio_hist_low = Hist.divide(
                total_model - low_band_model, total_model)
            fig.cd('ratio')
            error_band = rootpy_utils.get_band(
                ratio_hist_high, ratio_hist_low)
            error_band.fillstyle = '/'
            error_band.fillcolor = '#858585'
            error_band.drawstyle = 'E2'
            fig.draw('ratio', [error_band, ratio_hist], xdivisions=xdivisions)

    if separate_legends:
        fig.cd('main')
        right_legend = Legend(len(signal) + 1 if signal is not None else 1,
            pad=fig.pad('main'),
            **legend_params('right', textsize))
        right_legend.AddEntry(data, style='lep')
        if signal is not None:
            for s in reversed(scaled_signal):
                right_legend.AddEntry(s, style='F' if fill_signal else 'L')
        legends.append(right_legend)
        if model is not None:
            n_entries = len(model)
            if systematics:
                n_entries += 1
            model_legend = Legend(n_entries,
                pad=fig.pad('main'),
                **legend_params('left', textsize))
            for hist in reversed(model):
                model_legend.AddEntry(hist, style='F')
            if systematics:
                model_err_band = error_band_model.Clone()
                model_err_band.linewidth = 0
                model_err_band.linecolor = 'white'
                model_err_band.fillcolor = '#858585'
                model_err_band.title = 'Uncert.'
                model_legend.AddEntry(model_err_band, style='F')
            legends.append(model_legend)
    else:
        n_entries = 1
        if signal is not None:
            n_entries += len(scaled_signal)
        if model is not None:
            n_entries += len(model)
            if systematics:
                n_entries += 1
        fig.cd('main')
        legend = Legend(
            n_entries,
            pad=fig.pad('main'),
            **legend_params(legend_position, textsize))
        if data is not None:
            legend.AddEntry(data, style='lep')
        if signal is not None:
            for s in reversed(scaled_signal):
                legend.AddEntry(s, style='F' if fill_signal else 'L')
        if model:
            for hist in reversed(model):
                legend.AddEntry(hist, style='F')
            model_err_band = error_band_model.Clone()
            model_err_band.linewidth = 0
            model_err_band.linecolor = 'white'
            model_err_band.fillcolor = '#858585'
            model_err_band.title = 'Uncert.'
            legend.AddEntry(model_err_band, style='F')
        legends.append(legend)

    # draw the objects
    bounds = fig.draw('main', objects, ypadding=ypadding,
                      logy_crop_value=1E-1,
                      xdivisions=xdivisions)
    xaxis, yaxis = fig.axes('main')

    base_xaxis = xaxis
    base_xaxis.range_user = template.bounds()
    base_xaxis.limits = template.bounds()
    xmin, xmax, ymin, ymax = bounds

    if show_ratio:
        base_xaxis = fig.axes('ratio')[0]
        base_xaxis.range_user = template.bounds()
        base_xaxis.limits = template.bounds()

    # draw the legends
    fig.cd('main')
    for legend in legends:
        legend.Draw()

    label_plot(fig.pad('main'), template=template,
               xaxis=base_xaxis, yaxis=yaxis,
               xlabel=name, ylabel=ylabel, units=units,
               category_label=category.label,
               extra_label=plot_label,
               extra_label_position='right' if legend_position == 'left' else 'left',
               data_info=data_info)

    if logy and logy_min is not None:
        yaxis.min = logy_min
        ymin = logy_min

    # draw arrows
    if arrow_values is not None:
        arrow_top = ymin + (ymax - ymin) / 2.
        fig.cd('main')
        for value in arrow_values:
            arrow = Arrow(value, arrow_top, value, ymin, 0.05, '|>')
            arrow.SetAngle(30)
            arrow.SetLineWidth(2)
            arrow.Draw()

    if show_pvalue and data is not None and model:
        fig.cd('main')
        total_model = sum(model)
        # show p-value and chi^2
        pvalue = total_model.Chi2Test(data, 'WW')
        pvalue_label = ROOT.TLatex(
            0.6, 0.97,
            "p-value={0:.2f}".format(pvalue))
        pvalue_label.SetNDC(True)
        pvalue_label.SetTextFont(43)
        pvalue_label.SetTextSize(16)
        pvalue_label.Draw()
        chi2 = total_model.Chi2Test(data, 'WW CHI2/NDF')
        chi2_label = ROOT.TLatex(
            0.78, 0.97,
            "#chi^{{2}}/ndf={0:.2f}".format(chi2))
        chi2_label.SetNDC(True)
        chi2_label.SetTextFont(43)
        chi2_label.SetTextSize(16)
        chi2_label.Draw()

    if top_label is not None:
        fig.cd('main')
        label = ROOT.TLatex(
            fig.pad('main').GetLeftMargin() + 0.08, 0.97,
            top_label)
        label.SetNDC(True)
        label.SetTextFont(43)
        label.SetTextSize(16)
        label.Draw()

    if output_name is not None:
        # create the output filename
        filename = 'var_{0}_{1}'.format(
            category.name,
            output_name.lower().replace(' ', '_'))
        if logy:
            filename += '_logy'
        filename += '_root'
        if output_formats is None:
            output_formats = ('png',)
        # save the figure
        save_canvas(fig, output_dir, filename, formats=output_formats)

    return fig
