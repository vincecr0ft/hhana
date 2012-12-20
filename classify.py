from logger import log; log = log[__name__]

import numpy as np
# for reproducibilty
# especially for test/train set selection
np.random.seed(1987) # my birth year ;)

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter

from rootpy.plotting import Hist
from rootpy.io import open as ropen

from samples import *
from utils import draw


def search_flat_bins(bkg_scores, min_score, max_score, bins):

    scores = []
    weights = []
    for bkg, scores_dict in bkg_scores:
        s, w = scores_dict['NOMINAL']
        scores.append(s)
        weights.append(w)
    scores = np.concatenate(scores)
    weights = np.concatenate(weights)

    selection = (min_score <= scores) & (scores < max_score)
    scores = scores[selection]
    weights = weights[selection]

    sort_idx = np.argsort(scores)
    scores = scores[sort_idx]
    weights = weights[sort_idx]

    total_weight = weights.sum()
    bin_width = total_weight / bins

    # inefficient linear search for now
    weights_cumsum = np.cumsum(weights)
    boundaries = [min_score]
    curr_total = bin_width
    for i, cs in enumerate(weights_cumsum):
        if cs >= curr_total:
            boundaries.append((scores[i] + scores[i+1])/2)
            curr_total += bin_width
        if len(boundaries) == bins:
            break
    boundaries.append(max_score)
    return boundaries


def plot_grid_scores(
        grid_scores, best_point, params, name,
        label_all_bins=False,
        label_all_ticks=False,
        n_ticks=10,
        title=None,
        format='png'):

    param_names = sorted(grid_scores[0][0].keys())
    param_values = dict([(pname, []) for pname in param_names])
    for pvalues, score, cv_scores in grid_scores:
        for pname in param_names:
            param_values[pname].append(pvalues[pname])

    # remove duplicates
    for pname in param_names:
        param_values[pname] = np.unique(param_values[pname]).tolist()

    scores = np.empty(shape=[len(param_values[pname]) for pname in param_names])

    for pvalues, score, cv_scores in grid_scores:
        index = []
        for pname in param_names:
            index.append(param_values[pname].index(pvalues[pname]))
        scores.itemset(tuple(index), score)

    fig = plt.figure(figsize=(7, 5), dpi=100)
    ax = plt.axes([.12, .15, .8, .75])
    cmap = cm.get_cmap('jet', 100)
    img = ax.imshow(scores, interpolation="nearest", cmap=cmap,
            aspect='auto',
            origin='lower')

    if label_all_ticks:
        plt.xticks(range(len(param_values[param_names[1]])),
                param_values[param_names[1]])
        plt.yticks(range(len(param_values[param_names[0]])),
                param_values[param_names[0]])
    else:
        trees = param_values[param_names[1]]
        def tree_formatter(x, pos):
            if x < 0 or x >= len(trees):
                return ''
            return str(trees[int(x)])

        leaves = param_values[param_names[0]]
        def leaf_formatter(x, pos):
            if x < 0 or x >= len(leaves):
                return ''
            return str(leaves[int(x)])

        ax.xaxis.set_major_formatter(FuncFormatter(tree_formatter))
        ax.yaxis.set_major_formatter(FuncFormatter(leaf_formatter))
        ax.xaxis.set_major_locator(MaxNLocator(n_ticks, integer=True,
            prune='lower', steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_locator(MaxNLocator(n_ticks, integer=True,
            steps=[1, 2, 5, 10]))
        xlabels = ax.get_xticklabels()
        for label in xlabels:
            label.set_rotation(45)

    ax.set_xlabel(params[param_names[1]], fontsize=12,
            position=(1., 0.), ha='right')
    ax.set_ylabel(params[param_names[0]], fontsize=12,
            position=(0., 1.), va='top')

    ax.set_frame_on(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    for row in range(scores.shape[0]):
        for col in range(scores.shape[1]):
            decor={}
            if ((param_values[param_names[0]].index(best_point[param_names[0]])
                 == row) and
                (param_values[param_names[1]].index(best_point[param_names[1]])
                 == col)):
                decor = dict(weight='bold',
                             bbox=dict(boxstyle="round,pad=0.5",
                                       ec='black',
                                       fill=False))
            if label_all_bins or decor:
                plt.text(col, row, "%.3f" % (scores[row][col]), ha='center',
                         va='center', **decor)
    if title:
        plt.suptitle(title)

    plt.colorbar(img, fraction=.06, pad=0.03)
    plt.axis("tight")
    plt.savefig("grid_scores_%s.%s" % (name, format), bbox_inches='tight')
    plt.clf()


def plot_clf(
        background_scores,
        category,
        category_name,
        signal_scores=None,
        signal_scale=1.,
        data_scores=None,
        name=None,
        draw_histograms=True,
        draw_data=False,
        save_histograms=False,
        bins=10,
        min_score=0,
        max_score=1,
        signal_on_top=False,
        signal_colour_map=cm.spring,
        plot_signal_significance=True,
        fill_signal=True,
        systematics=None,
        **kwargs):

    if hasattr(bins, '__iter__'):
        # variable width bins
        hist_template = Hist(bins)
        min_score = min(bins)
        max_score = max(bins)
    else:
        hist_template = Hist(bins, min_score, max_score)

    bkg_hists = []
    for bkg, scores_dict in background_scores:
        hist = hist_template.Clone(title=bkg.label)
        scores, weight = scores_dict['NOMINAL']
        hist.fill_array(scores, weight)
        hist.decorate(**bkg.hist_decor)
        hist.systematics = {}
        for sys_term in scores_dict.keys():
            if sys_term == 'NOMINAL':
                continue
            sys_hist = hist_template.Clone()
            scores, weight = scores_dict[sys_term]
            sys_hist.fill_array(scores, weight)
            hist.systematics[sys_term] = sys_hist
        bkg_hists.append(hist)

    if signal_scores is not None:
        sig_hists = []
        for sig, scores_dict in signal_scores:
            sig_hist = hist_template.Clone(title=sig.label)
            scores, weight = scores_dict['NOMINAL']
            sig_hist.fill_array(scores, weight)
            sig_hist.decorate(**sig.hist_decor)
            sig_hist.systematics = {}
            for sys_term in scores_dict.keys():
                if sys_term == 'NOMINAL':
                    continue
                sys_hist = hist_template.Clone()
                scores, weight = scores_dict[sys_term]
                sys_hist.fill_array(scores, weight)
                sig_hist.systematics[sys_term] = sys_hist
            sig_hists.append(sig_hist)
    else:
        sig_hists = None

    if data_scores is not None and draw_data:
        data, data_scores = data_scores
        data_hist = hist_template.Clone(title=data.label)
        data_hist.decorate(**data.hist_decor)
        data_hist.fill_array(data_scores)
        log.info("Data events: %d" % sum(data_hist))
        log.info("Model events: %f" % sum(sum(bkg_hists)))
        for hist in bkg_hists:
            log.info("{0} {1}".format(hist.GetTitle(), sum(hist)))
        log.info("Data / Model: %f" % (sum(data_hist) / sum(sum(bkg_hists))))
    else:
        data_hist = None

    if draw_histograms:
        output_name = 'event_bdt_score'
        if name is not None:
            output_name += '_' + name
        draw(data=data_hist,
             model=bkg_hists,
             signal=sig_hists,
             signal_scale=signal_scale,
             plot_signal_significance=plot_signal_significance,
             category=category,
             category_name=category_name,
             name="BDT Score",
             output_name=output_name,
             range=(min_score, max_score),
             show_ratio=data_hist is not None,
             model_colour_map=None,
             signal_colour_map=signal_colour_map,
             signal_on_top=signal_on_top,
             fill_signal=fill_signal,
             systematics=systematics,
             **kwargs)
    return bkg_hists, sig_hists, data_hist


def std(X):

    return (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)


def make_classification(
        signals,
        backgrounds,
        category,
        region,
        branches,
        train_fraction,
        cuts=None,
        max_sig_train=None,
        max_bkg_train=None,
        max_sig_test=None,
        max_bkg_test=None,
        norm_sig_to_bkg_train=True,
        norm_sig_to_bkg_test=False,
        same_size_train=True,
        same_size_test=False,
        standardize=False,
        remove_negative_train_weights=False,
        systematic='NOMINAL'):

    signal_train_arrs = []
    signal_weight_train_arrs = []
    signal_test_arrs = []
    signal_weight_test_arrs = []

    for signal in signals:
        train, test = signal.train_test(
            category=category,
            region=region,
            branches=branches,
            train_fraction=train_fraction,
            cuts=cuts,
            systematic=systematic)
        signal_weight_train_arrs.append(train['weight'])
        signal_weight_test_arrs.append(test['weight'])

        signal_train_arrs.append(
            np.vstack(train[branch] for branch in branches).T)
        signal_test_arrs.append(
            np.vstack(test[branch] for branch in branches).T)

    background_train_arrs = []
    background_weight_train_arrs = []
    background_test_arrs = []
    background_weight_test_arrs = []

    for background in backgrounds:
        train, test = background.train_test(
            category=category,
            region=region,
            branches=branches,
            train_fraction=train_fraction,
            cuts=cuts,
            systematic=systematic)
        background_weight_train_arrs.append(train['weight'])
        background_weight_test_arrs.append(test['weight'])

        background_train_arrs.append(
            np.vstack(train[branch] for branch in branches).T)
        background_test_arrs.append(
            np.vstack(test[branch] for branch in branches).T)

    signal_train = np.concatenate(signal_train_arrs)
    signal_weight_train = np.concatenate(signal_weight_train_arrs)
    signal_test = np.concatenate(signal_test_arrs)
    signal_weight_test = np.concatenate(signal_weight_test_arrs)

    background_train = np.concatenate(background_train_arrs)
    background_weight_train = np.concatenate(background_weight_train_arrs)
    background_test = np.concatenate(background_test_arrs)
    background_weight_test = np.concatenate(background_weight_test_arrs)

    if remove_negative_train_weights:
        # remove samples from the training sample with a negative weight
        signal_train = signal_train[signal_weight_train >= 0]
        background_train = background_train[background_weight_train >= 0]

        signal_weight_train = signal_weight_train[signal_weight_train >= 0]
        background_weight_train = background_weight_train[background_weight_train >= 0]

    if max_sig_train is not None and max_sig_train < len(signal_train):
        subsample = np.random.permutation(len(signal_train))[:max_sig_train]
        signal_train = signal_train[subsample]
        signal_weight_train = signal_weight_train[subsample]

    if max_bkg_train is not None and max_bkg_train < len(background_train):
        subsample = np.random.permutation(len(background_train))[:max_bkg_train]
        background_train = background_train[subsample]
        background_weight_train = background_weight_train[subsample]

    if max_sig_test is not None and max_sig_test < len(signal_test):
        subsample = np.random.permutation(len(signal_test))[:max_sig_test]
        signal_test = signal_test[subsample]
        signal_weight_test = signal_weight_test[subsample]

    if max_bkg_test is not None and max_bkg_test < len(background_test):
        subsample = np.random.permutation(len(background_test))[:max_bkg_test]
        background_test = background_test[subsample]
        background_weight_test = background_weight_test[subsample]

    if same_size_train:
        if len(background_train) > len(signal_train):
            # random subsample of background so it's the same size as signal
            subsample = np.random.permutation(
                len(background_train))[:len(signal_train)]
            background_train = background_train[subsample]
            background_weight_train = background_weight_train[subsample]
        elif len(background_train) < len(signal_train):
            # random subsample of signal so it's the same size as background
            subsample = np.random.permutation(
                len(signal_train))[:len(background_train)]
            signal_train = signal_train[subsample]
            signal_weight_train = signal_weight_train[subsample]

    if same_size_test:
        if len(background_test) > len(signal_test):
            # random subsample of background so it's the same size as signal
            subsample = np.random.permutation(
                len(background_test))[:len(signal_test)]
            background_test = background_test[subsample]
            background_weight_test = background_weight_test[subsample]
        elif len(background_test) < len(signal_test):
            # random subsample of signal so it's the same size as background
            subsample = np.random.permutation(
                len(signal_test))[:len(background_test)]
            signal_test = signal_test[subsample]
            signal_weight_test = signal_weight_test[subsample]

    if norm_sig_to_bkg_train:
        # normalize signal to background
        signal_weight_train *= (
            background_weight_train.sum() / signal_weight_train.sum())

    if norm_sig_to_bkg_test:
        # normalize signal to background
        signal_weight_test *= (
            background_weight_test.sum() / signal_weight_test.sum())

    log.info("Training Samples:")
    log.info("Signal: %d events, %s features" % signal_train.shape)
    log.info("Sum(signal weights): %f" % signal_weight_train.sum())
    log.info("Background: %d events, %s features" % background_train.shape)
    log.info("Sum(background weight): %f" % background_weight_train.sum())
    log.info("")
    log.info("Test Samples:")
    log.info("Signal: %d events, %s features" % signal_test.shape)
    log.info("Sum(signal weights): %f" % signal_weight_test.sum())
    log.info("Background: %d events, %s features" % background_test.shape)
    log.info("Sum(background weight): %f" % background_weight_test.sum())

    # create training/testing samples
    sample_train = np.concatenate((background_train, signal_train))
    sample_test = np.concatenate((background_test, signal_test))

    sample_weight_train = np.concatenate(
        (background_weight_train, signal_weight_train))
    sample_weight_test = np.concatenate(
        (background_weight_test, signal_weight_test))

    if standardize:
        sample_train = std(sample_train)
        sample_test = std(sample_test)

    labels_train = np.concatenate(
        (np.zeros(len(background_train)), np.ones(len(signal_train))))
    labels_test = np.concatenate(
        (np.zeros(len(background_test)), np.ones(len(signal_test))))

    # random permutation of training sample
    perm = np.random.permutation(len(labels_train))
    sample_train = sample_train[perm]
    sample_weight_train = sample_weight_train[perm]
    labels_train = labels_train[perm]

    return sample_train, sample_test,\
        sample_weight_train, sample_weight_test,\
        labels_train, labels_test


def staged_score(self, X, y, sample_weight, n_estimators=-1):
    """
    calculate maximum signal significance
    """
    bins = 50
    for p in self.staged_predict_proba(X, n_estimators=n_estimators):

        scores = p[:,-1]

        # weighted mean accuracy
        y_pred = scores >= .5
        acc = np.average((y_pred == y), weights=sample_weight)

        min_score, max_score = scores.min(), scores.max()
        b_hist = Hist(bins, min_score, max_score + 0.0001)
        s_hist = b_hist.Clone()

        scores_s, w_s = scores[y==1], sample_weight[y==1]
        scores_b, w_b = scores[y==0], sample_weight[y==0]

        # fill the histograms
        s_hist.fill_array(scores_s, w_s)
        b_hist.fill_array(scores_b, w_b)

        # reverse cumsum
        #bins = list(b_hist.xedges())[:-1]
        s_counts = np.array(s_hist)
        b_counts = np.array(b_hist)
        S = s_counts[::-1].cumsum()[::-1]
        B = b_counts[::-1].cumsum()[::-1]

        # S / sqrt(S + B)
        s_sig = np.divide(list(S), np.sqrt(list(S + B)))

        #max_bin = np.argmax(np.ma.masked_invalid(significance)) #+ 1
        #max_sig = significance[max_bin]
        #max_cut = bins[max_bin]

        s_sig_max = np.max(np.ma.masked_invalid(s_sig))
        yield s_sig_max * acc


def write_score_hists(f, mass, scores_list, hist_template, no_neg_bins=True):

    sys_hists = {}
    for samp, scores_dict in scores_list:
        for sys_term, (scores, weights) in scores_dict.items():
            if sys_term == 'NOMINAL':
                suffix = ''
            else:
                suffix = '_' + '_'.join(sys_term)
            hist = hist_template.Clone(
                    name=samp.name + ('_%d' % mass) + suffix)
            hist.fill_array(scores, weights)
            if sys_term not in sys_hists:
                sys_hists[sys_term] = []
            sys_hists[sys_term].append(hist)
    f.cd()
    for sys_term, hists in sys_hists.items():
        bad_bins = []
        if no_neg_bins:
            # check for negative bins over all systematics and zero them out
            # negative bins cause lots of problem in the limit setting
            # negative bin contents effectively means
            # the same as "no events here..."
            total_hist = sum(hists)
            for bin, content in enumerate(total_hist):
                if content < 0:
                    log.warning("Found negative bin %d (%f) for systematic %s" % (
                            bin, content, sys_term))
                    bad_bins.append(bin)
        for hist in hists:
            for bin in bad_bins:
                # zero out bad bins
                hist[bin] = 0.
            hist.Write()
