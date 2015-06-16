import ROOT
from itertools import izip
from matplotlib import cm
from rootpy.plotting.style.atlas.labels import ATLAS_label
from rootpy.memory.keepalive import keepalive
from .. import ATLAS_LABEL


def set_colors(hists, colors='jet'):
    if isinstance(colors, basestring):
        colors = cm.get_cmap(colors, len(hists))
    if hasattr(colors, '__call__'):
        for i, h in enumerate(hists):
            color = colors((i + 1) / float(len(hists) + 1))
            h.SetColor(color)
    else:
        for h, color in izip(hists, colors):
            h.SetColor(color)


def category_lumi_atlas(pad, category_label=None,
                        data_info=None, atlas_label=None,
                        textsize=22):

    left, right, bottom, top = pad.margin_pixels
    height = float(pad.height_pixels)

    # draw the category label
    if category_label:
        label = ROOT.TLatex(
            1. - pad.GetRightMargin(),
            1. - (textsize - 2) / height,
            category_label)
        label.SetNDC()
        label.SetTextFont(43)
        label.SetTextSize(textsize)
        label.SetTextAlign(31)
        with pad:
            label.Draw()
        keepalive(pad, label)

    # draw the luminosity label
    if data_info is not None:
        plabel = ROOT.TLatex(
            1. - pad.GetRightMargin() - 0.03,
            1. - (top + textsize + 15) / height,
            str(data_info))
        plabel.SetNDC()
        plabel.SetTextFont(43)
        plabel.SetTextSize(textsize)
        plabel.SetTextAlign(31)
        with pad:
            plabel.Draw()
        keepalive(pad, plabel)

    # draw the ATLAS label
    if atlas_label is not False:
        label = atlas_label or ATLAS_LABEL
        ATLAS_label(pad.GetLeftMargin() + 0.03,
                    1. - (top + textsize + 15) / height,
                    sep=0.132, pad=pad, sqrts=None,
                    text=label,
                    textsize=textsize)
    pad.Update()
    pad.Modified()


def label_plot(pad, template, xaxis, yaxis,
               ylabel='Events', xlabel=None,
               units=None, data_info=None,
               category_label=None,
               atlas_label=None,
               extra_label=None,
               extra_label_position='left',
               textsize=22):

    # set the axis labels
    binw = list(template.xwidth())
    binwidths = list(set(['%.2g' % w for w in binw]))
    if units is not None:
        if xlabel is not None:
            xlabel = '%s [%s]' % (xlabel, units)
        if ylabel and len(binwidths) == 1 and binwidths[0] != '1':
            # constant width bins
            ylabel = '%s / %s %s' % (ylabel, binwidths[0], units)
    elif ylabel and len(binwidths) == 1 and binwidths[0] != '1':
        ylabel = '%s / %s' % (ylabel, binwidths[0])

    if ylabel:
        yaxis.SetTitle(ylabel)
    if xlabel:
        xaxis.SetTitle(xlabel)

    left, right, bottom, top = pad.margin_pixels
    height = float(pad.height_pixels)

    category_lumi_atlas(pad, category_label, data_info, atlas_label)

    # draw the extra label
    if extra_label is not None:
        if extra_label_position == 'left':
            label = ROOT.TLatex(pad.GetLeftMargin() + 0.03,
                                1. - (top + 2 * (textsize + 15)) / height,
                                extra_label)
        else: # right
            label = ROOT.TLatex(1. - pad.GetRightMargin() - 0.03,
                                1. - (top + 2 * (textsize + 15)) / height,
                                extra_label)
            label.SetTextAlign(31)
        label.SetNDC()
        label.SetTextFont(43)
        label.SetTextSize(textsize)
        with pad:
            label.Draw()
        keepalive(pad, label)

    pad.Update()
    pad.Modified()


def legend_params(position, textsize):
    location = 'upper {0}'.format(position)
    if location == 'upper left':
        print "UPPER LEFT"
        print "text size ",textsize
        return dict(
            leftmargin=0.05,
            rightmargin=0.5,
            topmargin=0.05,
            # margin=0.25,
            entrysep=0.005,
            #entryheight=textsize + 0.04,
            textsize=textsize*0.6)
    elif location == 'upper right':
        print "UPPER RIGHT"
        print "text size ",textsize
        return dict(
            leftmargin=0.05,
            rightmargin=0.5,
            topmargin=0.05,
            # margin=0.25,
            entrysep=0.005,
            #entryheight=textsize + 0.4,
            textsize=textsize*0.6)
    else:
        raise RuntimeError('wrong position argument -- legend_params should be updated')
