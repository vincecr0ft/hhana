
import os
import math

from rootpy.plotting import Hist, Hist2D
from rootpy.io import root_open
from rootpy.stats import histfactory
from rootpy.utils.path import mkdir_p

from statstools.histfactory import (
    to_uniform_binning, apply_remove_window, is_signal)

from . import log; log = log[__name__]
from . import CONST_PARAMS, CACHE_DIR, MMC_MASS, POI
from .categories import CATEGORIES
from .plotting import hist_scores

import pickle
import os


def write_workspaces(path, prefix, year_mass_category_channel,
                     controls=None, shapeControls=None,
                     silence=False):
    log.info("writing workspaces ... for {0}".format(str(year_mass_category_channel)))
    if controls is None:
        controls = []
    if shapeControls is None:
        shapeControls = []

    if not os.path.exists(path):
        mkdir_p(path)
    for year, mass_category_channel in year_mass_category_channel.items():
        # write workspaces for each year
        for mass, category_channel in mass_category_channel.items():
            if isinstance(controls, dict):
                log.info("controls are dicks")
                if isinstance(controls[year], dict):
                    log.info("controls by year are dicks")
                    mass_controls = controls[year][mass].values()
                else:
                    log.info("controls years are dicks")
                    mass_controls = controls[year]
            else:
                mass_controls = controls
            if isinstance(shapeControls, dict):
                log.info("shapeControls are dicks")
                if isinstance(shapeControls[year], dict):
                    log.info("shapeControls by year are dicks")
                    shape_controls = shapeControls[year][mass].values()
                else:
                    log.info("shapeControls years are dicks")
                    shape_controls = shapeControls[year]
            else:
                shape_controls = shapeControls

            log.info("comparing {0} with {1}".format(str(mass_controls),str(shape_controls)))
            channels = []
            # make workspace for each category
            # include the control region in each
            nCategories=0
            for category, channel in category_channel.items():
                nCategories+=1
                name = "{0}_{1}_{2}_{3}".format(
                    prefix, year % 1000, category, mass)
                log.info("writing {0} ...".format(name))
                if mass<0.:
                    parity='m'
                else:
                    parity='p'
                newname = "AllSys_cp_{}_0_{:02.0f}_{}".format(parity, abs(mass*100), prefix)

                # print newname
                measurement = histfactory.make_measurement(
                    newname, [channel] + mass_controls + shape_controls,
                    POI=POI,
                    const_params=CONST_PARAMS)
                workspace = histfactory.make_workspace(measurement, name=newname,
                                                       silence=silence)
                with root_open(os.path.join(path, '{0}.root'.format(newname)),
                               'recreate') as workspace_file:
                    workspace.Write()
                    # mu=1 for Asimov data
#                    measurement.SetParamValue('ATLAS_epsilon', 1)
                    histfactory.write_measurement(measurement,
                        root_file=workspace_file,
                        xml_path=os.path.join(path, newname),
                        silence=silence)
                channels.append(channel)            
            log.info("length of channels is {0}".format(str(len(channels))))
            # make combined workspace
"""
            name = "{0}_{1}_combination_{2}".format(prefix, year % 1000, mass)
            log.info("writing {0} ...".format(name))
            measurement = histfactory.make_measurement(
                name, channels + mass_controls,
                POI=POI,
                const_params=CONST_PARAMS)
            workspace = histfactory.make_workspace(measurement, name=name,
                                                   silence=silence)
            with root_open(os.path.join(path, '{0}.root'.format(name)),
                           'recreate') as workspace_file:
                workspace.Write()
                # mu=1 for Asimov data
                #measurement.SetParamValue('SigXsecOverSM', 1)
                histfactory.write_measurement(measurement,
                    root_file=workspace_file,
                    xml_path=os.path.join(path, name),
                    silence=silence)
    # write combined workspaces over all years
    years = year_mass_category_channel.keys()
    if len(years) == 1:
        return
    masses = year_mass_category_channel[years[0]].keys()
    categories = year_mass_category_channel[years[0]][masses[0]].keys()
    for mass in masses:
        if isinstance(controls, dict):
            if isinstance(controls[year], dict):
                mass_controls = [control for year in years
                                 for control in controls[year][mass].values()]
            else:
                mass_controls = [control for year in years
                                 for control in controls[year]]
        else:
            mass_controls = controls
        channels = []
        # make workspace for each category
        # include the control region in each
        # TODO: categories might be different across years

        for category in categories:
            cat_channels = [year_mass_category_channel[year][mass][category]
                            for year in years]
            name = "{0}_full_{1}_{2}".format(
                prefix, category, mass)
            log.info("writing {0} ...".format(name))
            # make workspace
            measurement = histfactory.make_measurement(
                name, cat_channels + mass_controls,
                POI=POI,
                const_params=CONST_PARAMS)
            workspace = histfactory.make_workspace(measurement, name=name,
                                                   silence=silence)
            with root_open(os.path.join(path, '{0}.root'.format(name)),
                           'recreate') as workspace_file:
                workspace.Write()
                # mu=1 for Asimov data
                #measurement.SetParamValue('SigXsecOverSM', 1)
                histfactory.write_measurement(measurement,
                    root_file=workspace_file,
                    xml_path=os.path.join(path, name),
                    silence=silence)
            channels.extend(cat_channels)

        channels = [chan for year in years
                    for chan in year_mass_category_channel[year][mass].values()]
        # make combined workspace
        name = "{0}_full_combination_{1}".format(prefix, mass)
        log.info("writing {0} ...".format(name))
        measurement = histfactory.make_measurement(
            name, channels + mass_controls,
            POI=POI,
            const_params=CONST_PARAMS)
        workspace = histfactory.make_workspace(measurement, name=name,
                                               silence=silence)
        with root_open(os.path.join(path, '{0}.root'.format(name)),
                       'recreate') as workspace_file:
            workspace.Write()
            # mu=1 for Asimov data
            #measurement.SetParamValue('SigXsecOverSM', 1)
            histfactory.write_measurement(measurement,
                root_file=workspace_file,
                xml_path=os.path.join(path, name),
                silence=silence)
"""

def mva_workspace(analysis, categories, masses,
                  clf_mass=None,
                  clf_bins='optimal',
                  clf_swap=False,
                  unblind=False,
                  systematics=False,
                  cuts=None):
    hist_template = Hist(5, 0, 1.5, type='D')
    controls = analysis.make_var_channels(
        hist_template, 'dEta_tau1_tau2',
        CATEGORIES['mva_workspace_controls'],
        analysis.target_region,
        include_signal=True, masses=masses,
        systematics=systematics)
    mass_category_channel = {}
    for category in analysis.iter_categories(categories):
        for mass in masses:
            clf = analysis.get_clf(category, load=True,
                                   mass=clf_mass or mass,
                                   transform=True,
                                   swap=clf_swap)
            if isinstance(clf_bins, basestring):
                if clf_bins == 'optimal':
                    # get the binning (see the optimize-binning script)
                    bins = clf.binning(analysis.year, overflow=1E5)
                    log.info("binning: {0}".format(str(bins)))
                else:
                    bins = int(clf_bins)
            else:
                bins = clf_bins
            # construct a "channel" for each mass point
            scores, channel = analysis.clf_channels(
                clf, category,
                region=analysis.target_region,
                bins=bins,
                mass=mass,
                mode='workspace',
                systematics=systematics,
                cuts=cuts,
                unblind=unblind or 2,
                hybrid_data=not unblind,
                uniform=True, mva=True)
            if mass not in mass_category_channel:
                mass_category_channel[mass] = {}
            mass_category_channel[mass][category.name] = channel
    return mass_category_channel, controls


def cuts_workspace(analysis, categories, masses,
                   unblind=False,
                   systematics=False,
                   cuts=None,
                   sideband=False):
    hybrid_data = None if unblind else {MMC_MASS:(100., 150.)}
    channels = {}
    for category in analysis.iter_categories(categories):
        binning = category.limitbins
        if isinstance(binning, dict):
            binning = binning[analysis.year]
        hist_template = Hist(binning, type='D')
        for mass in masses:
            channel = analysis.get_channel_array(
                {MMC_MASS: hist_template},
                category=category,
                region=analysis.target_region,
                cuts=cuts,
                include_signal=True,
                mass=mass,
                mode='workspace',
                systematics=systematics,
                hybrid_data=hybrid_data,
                uniform=False)[MMC_MASS]
            if sideband:
                # remove the signal window
                remove_window = (100, 150)
                channel.data.hist = apply_remove_window(
                    channel.data.hist, remove_window)
                for s in channel.samples:
                    s.hist = apply_remove_window(
                        s.hist, remove_window)
                    for histosys in s.histo_sys:
                        histosys.high = apply_remove_window(
                            histosys.high, remove_window)
                        histosys.low = apply_remove_window(
                            histosys.low, remove_window)
            # convert to uniform binning
            channel.data.hist = to_uniform_binning(channel.data.hist)
            for s in channel.samples:
                s.hist = to_uniform_binning(s.hist)
                for histosys in s.histo_sys:
                    histosys.high = to_uniform_binning(histosys.high)
                    histosys.low = to_uniform_binning(histosys.low)
            if mass not in channels:
                channels[mass] = {}
            channels[mass][category.name] = channel
    return channels, []


def feature_workspace(field, template,
                      analysis, categories, masses,
                      systematics=False,
                      cuts=None):
    channels = {}
    for category in analysis.iter_categories(categories):
        for mass in masses:
            channel = analysis.get_channel_array(
                {field: template},
                category=category,
                region=analysis.target_region,
                cuts=cuts,
                include_signal=True,
                mass=mass,
                mode='workspace',
                systematics=systematics,
                uniform=False)[field]
            if mass not in channels:
                channels[mass] = {}
            channels[mass][category.name] = channel
    return channels, []


def weighted_mass_workspace(analysis, categories, masses,
                            systematics=False,
                            cuts=None):
    hist_template = Hist(20, 50, 250, type='D')
    channels = {}
    for category in analysis.iter_categories(categories):
        clf = analysis.get_clf(category, load=True, mass=125)
        clf_bins = clf.binning(analysis.year, overflow=1E5)
        scores = analysis.get_scores(
            clf, category, analysis.target_region,
            masses=[125], mode='combined',
            systematics=False,
            unblind=True)
        bkg_scores = scores.bkg_scores
        sig_scores = scores.all_sig_scores[125]
        min_score = scores.min_score
        max_score = scores.max_score
        bkg_score_hist = Hist(clf_bins, type='D')
        sig_score_hist = bkg_score_hist.Clone()
        hist_scores(bkg_score_hist, bkg_scores)
        _bkg = bkg_score_hist.Clone()
        hist_scores(sig_score_hist, sig_scores)
        _sig = sig_score_hist.Clone()
        sob_hist = (1 + _sig / _bkg)
        _log = math.log
        for bin in sob_hist.bins(overflow=True):
            bin.value = _log(bin.value)
        log.info(str(list(sob_hist.y())))
        for mass in masses:
            channel = analysis.get_channel_array(
                {MMC_MASS: hist_template},
                category=category,
                region=analysis.target_region,
                include_signal=True,
                weight_hist=sob_hist,
                clf=clf,
                cuts=cuts,
                mass=mass,
                mode='workspace',
                systematics=systematics)[MMC_MASS]
            if mass not in channels:
                channels[mass] = {}
            channels[mass][category.name] = channel
    return channels, []

def weighted_mass_cba_workspace(analysis, categories, masses,
                                systematics=False,
                                cuts=None):
    hist_template = Hist(20, 50, 250, type='D')
    channels = {}

    def scaled(hist, factor):
        new_hist = hist * factor
        new_hist.name = hist.name + '_scaled'
        return new_hist

    for category in analysis.iter_categories(categories):
        for mass in masses:
            channel = analysis.get_channel_array(
                {MMC_MASS: hist_template},
                category=category,
                region=analysis.target_region,
                include_signal=True,
                cuts=cuts,
                mass=mass,
                mode='workspace',
                systematics=systematics)[MMC_MASS]
            # weight by ln(1 + s / b)
            total_s = hist_template.Clone()
            total_s.Reset()
            total_b = total_s.Clone()
            for sample in channel.samples:
                if is_signal(sample):
                    total_s += sample.hist
                else:
                    total_b += sample.hist
            sob = math.log(1 + total_s.integral() / total_b.integral())
            channel.data.hist = scaled(channel.data.hist, sob)
            for sample in channel.samples:
                sample.hist = scaled(sample.hist, sob)
                for hsys in sample.histo_sys:
                    hsys.high = scaled(hsys.high, sob)
                    hsys.low = scaled(hsys.low, sob)
            if mass not in channels:
                channels[mass] = {}
            channels[mass][category.name] = channel
    return channels, []


def weighted_mixing_workspace(analysis, categories, masses, mixings,
                            systematics=False,
                            cuts=None):
    hist_template = Hist(6, -10, 10, type='D')
    channels = {}
    for category in analysis.iter_categories(categories):
        clf = analysis.get_clf(category, load=True, mass=125)
        clf_bins = clf.binning(analysis.year, overflow=1E5)
        scores = analysis.get_scores(
            clf, category, analysis.target_region,
            masses=[125], mode='CP',
            systematics=False,
            unblind=True)
        bkg_scores = scores.bkg_scores
        sig_scores = scores.all_sig_scores[125]
        min_score = scores.min_score
        max_score = scores.max_score
        bkg_score_hist = Hist(clf_bins, type='D')
        sig_score_hist = bkg_score_hist.Clone()
        hist_scores(bkg_score_hist, bkg_scores)
        _bkg = bkg_score_hist.Clone()
        hist_scores(sig_score_hist, sig_scores)
        _sig = sig_score_hist.Clone()
        sob_hist = (1 + _sig / _bkg)
        _log = math.log
        for bin in sob_hist.bins(overflow=True):
            bin.value = _log(bin.value)
        log.info(str(list(sob_hist.y())))

        log.info("getting workspaces for {0}".format(str(mixings)))
        for mixing in mixings:
            if not mixing==0.0:
                comparison=[0.0,mixing]
            else:
                comparison=[0.0]
            channel = analysis.get_channel_array(
                {'o1': hist_template},
                category=category,
                region=analysis.target_region,
                include_signal=True,
                weight_hist=sob_hist,
                clf=clf,
                cuts=cuts,
                mass=125,
                mixings=comparison,
                mode='CPworkspace',
                systematics=systematics)['o1']
            if mixing not in channels:
                channels[mixing] = {}
            channels[mixing][category.name] = channel
    return channels, []

def mixing_workspace(analysis, categories, masses, mixings=None,
                            systematics=False,
                            cuts=None):
    hist_template = Hist([-15,-3.5,-1.5,0,1.5,3.5,15], type='D')
    channels = {}
    print "mixings now",mixings
    log.info("mixings have become {0}".format(mixings))
    for category in analysis.iter_categories(categories):
        log.info("mixing workspace for category {0}".format(category.name))
        clf = analysis.get_clf(category, load=True, mass=125)
        clf_bins = clf.binning(analysis.year, overflow=1E5)
        log.info("getting workspaces for {0}".format(str(mixings)))
        controls = {}
        shapeControls = {}
        for mixing in mixings:
            if not mixing==0.0:
                comparison=[0.0,mixing]
            else:
                comparison=[0.0]
            log.info("this comparison is {0}".format(str(comparison)))
            channel = analysis.get_channel_array(
                {'o1': hist_template},
                category=category,
                region=analysis.target_region,
                include_signal=True,
                clf=clf,
                cuts=cuts,
                mass=125,
                mixings=comparison,
                min_score=0.815967620756,#0.567454611796
                mode='CPworkspace',
                systematics=systematics)['o1']
            if mixing not in channels:
                channels[mixing] = {}                
            channels[mixing][category.name] = channel
            log.info("got SR moving in for CRs")

            #get crs
            cr_template = Hist(5, 0, 1.5, type='D')
            control = analysis.make_var_channels(
                cr_template, 'dEta_tau1_tau2',
                CATEGORIES['mva_workspace_controls'],
                analysis.target_region,
                mixings=comparison,
                mode='CPworkspace',
                include_signal=True, masses=masses,
                systematics=systematics)
            if mixing not in controls:
                controls[mixing] = {}
            controls[mixing] = control[125]

            bins = 10
            scr_template = Hist(8, -1, 0.6, type='D')
            shapeControl = analysis.make_var_channels(
                scr_template, clf,
                CATEGORIES['mva_vbf'],
                analysis.target_region,
                mixings=comparison,
                mode='CPworkspace',
                include_signal=True, masses=masses,
                systematics=systematics)
            if mixing not in shapeControls:
                shapeControls[mixing] = {}
            shapeControls[mixing] = shapeControl[125]

    return channels, controls, shapeControls

def mass2d_workspace(analysis, categories, masses,
                     systematics=False):
    hist_template = Hist2D(250, 0, 250, 200, -1, 1, type='D')
    channels = {}
    for category in analysis.iter_categories(categories):
        clf = analysis.get_clf(category, load=True)
        for mass in masses:
            channel = analysis.get_channel_array(
                {MMC_MASS: hist_template},
                category=category,
                region=analysis.target_region,
                clf=clf,
                include_signal=True,
                mass=mass,
                mode='workspace',
                systematics=systematics,
                ravel=False)[MMC_MASS]
            if mass not in channels:
                channels[mass] = {}
            channels[mass][category.name] = channel
    return channels, []
