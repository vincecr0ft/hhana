import math

BDT_BLIND = {
    2012: {
        'vbf': 3,
        'boosted': 1},
    2011: {
        'vbf': 2,
        'boosted': 2}
}


def get_label(name, units=True, latex=False):
    info = VARIABLES[name]
    if latex:
        label = info['title']
    else:
        label = info['root']
    if units and 'units' in info:
        label += ' [{0}]'.format(info['units'])
    return label


def get_binning(name, category, year):
    binning = VARIABLES[name]['binning']
#    binning = OBSERVABLES[name]['binning']
    if isinstance(binning, dict):
        if year in binning:
            binning = binning[year]
        if isinstance(binning, dict):
            binning = binning.get(category.name.upper(), binning[None])
    return binning


def get_scale(name):
    info = VARIABLES[name]
    return info.get('scale', 1)


def get_units(name):
    info = VARIABLES[name]
    return info.get('units', None)


def blind_hist(name, hist, year=None, category=None):
    if name.upper() == 'BDT':
        unblind = BDT_BLIND[year][category.name]
        hist[-(unblind + 1):] = (0, 0)
        return
    blind = VARIABLES[name].get('blind', None)
    if blind is None:
        return
    left, right = blind
    left_bin = hist.FindBin(left)
    right_bin = hist.FindBin(right)
    for ibin in xrange(left_bin, right_bin + 1):
        hist[ibin] = (0, 0)


WEIGHTS = {
    'pileup_weight': {
        'title': 'Pile-up Weight',
        'root': 'Pile-up Weight',
        'filename': 'weight_pileup',
        'binning': (50, -.2, 3.)
    },
}

YEAR_VARIABLES = {
    2011: {
        'RunNumber' : {
            'title': r'RunNumber',
            'root': 'Run Number',
            'filename': 'runnumber',
            'binning': [
                177531, 177986, 178163, 179710, 180614,
                182013, 182726, 183544, 185353, 186516,
                186873, 188902, 190503],
        }
    },
    2012: {
        'RunNumber' : {
            'title': r'RunNumber',
            'root': 'Run Number',
            'filename': 'runnumber',
            'binning': [
                200804, 202660, 206248, 207447, 209074, 210184, 211522,
                212619, 213431, 213900, 214281, 215414, 216399],
        }
    }
}

OBSERVABLES={
    'o1': {
        'title': r'0$_{1}$',
        'root': '#font[152]{0_1}(#font[52]{Optimal Observable 1})',
        'filename': 'o1',
        'binning': [-15,-3.5,-1.5,0,1.5,3.5,15],
        'cats': ['2J', 'VBF'],
        'legend': 'right',
    },
    'dphi': {
        'title': r'dPhi$_{jj}$',
        'root': '#font[152]{dY}jj(#font[52]{sign delta phi})',
        'filename': 'dphi',
        'binning': [-5,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,5,5.1],
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
}

VARIABLES = {
    'o1': {
        'title': r'0$_{1}$',
        'root': '#font[152]{0_1}(#font[52]{Optimal Observable 1})',
        'filename': 'o1',
        'binning': [-15,-3.5,-1.5,0,1.5,3.5,15],
        'cats': ['2J', 'VBF'],
        'legend': 'right',
    },
    'dphi': {
        'title': r'dPhi$_{jj}$',
        'root': '#font[152]{dY}jj(#font[52]{sign delta phi})',
        'filename': 'dphi',
        'binning': [-5,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,5,5.1],
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    #'ntrack_pv': {
    #    'title': r'Number of Tracks from the Primary Vertex',
    #    'root': '#font[52]{Number of Tracks from the Primary Vertex}',
    #    'filename': 'ntrack_pv',
    #    'bins': 10,
    #    'range': (0.5, 120.5),
    #},
    #'ntrack_nontau_pv': {
    #    'title': r'Number of Non-Tau Tracks from the Primary Vertex',
    #    'root': '#font[52]{Number of Non-Tau Tracks from the Primary Vertex}',
    #    'filename': 'ntrack_nontau_pv',
    #    'bins': 10,
    #    'range': (0.5, 120.5),
    #},
    'averageIntPerXing': {
        'title': r'$\langle\mu\rangle|_{LB,BCID}$',
        'root': '#font[152]{#LT#mu#GT#cbar}_{LB,BCID}',
        'filename': 'averageIntPerXing',
        'binning': (20, 0, 40),
    },
    #'actualIntPerXing': {
    #    'title': r'$\langle\mu\rangle|_{LB}(BCID)$',
    #    'root': '#font[152]{#LT#mu#GT#cbar}_{LB}#font[52]{(BCID)}',
    #    'filename': 'actualIntPerXing',
    #    'bins': 20,
    #    'range': (0, 40),
    #},
    'sum_pt': {
        'title': r'$\sum p_T$ Taus and Two Leading Jets',
        'root': '#font[152]{#sum} #font[52]{p}_{T} #font[52]{Taus and Two Leading Jets}',
        'filename': 'sum_pt',
        'binning': (20, 50, 550),
        'scale': 0.001,
        'units': 'GeV',
    },
    'sum_pt_full': {
        'title': r'$\sum p_T$ Taus and All Selected Jets',
        'root': '#font[152]{#sum} #font[52]{p}_{T} #font[52]{Taus and All Selected Jets}',
        'filename': 'sum_pt_full',
        'binning': (20, 50, 550),
        'scale': 0.001,
        'units': 'GeV',
    },
    'vector_sum_pt': {
        'title': r'$\sum \vec{p}_T$',
        'root': '#font[52]{p}_{T}^{Total}',
        'filename': 'vector_sum_pt',
        'binning': {
            'VBF': (20, 0, 100),
            None: (20, 0, 200)},
        'scale': 0.001,
        'units': 'GeV',
    },
    'numJets': {
        'title': r'Number of Selected Jets',
        'root': '#font[52]{Number of Selected Jets}',
        'filename': 'numjets',
        'binning': (7, -.5, 6.5),
        'integer': True,
    },
    'MET_et': {
        'title': r'$E^{miss}_{T}$',
        'root': '#font[52]{E}^{miss}_{T}',
        'filename': 'MET',
        'binning': {
            'PRESELECTION': (13, 15, 80),
            'REST': (13, 15, 80),
            None: (17, 15, 100)},
        'scale': 0.001,
        'units': 'GeV',
    },
    'MET_etx': {
        'title': r'$E^{miss}_{T_{x}}$',
        'root': '#font[52]{E}^{miss}_{T_{x}}',
        'filename': 'MET_x',
        'binning': (20, -75, 75),
        'scale': 0.001,
        'units': 'GeV',
        'legend': 'left',
    },
    'MET_ety': {
        'title': r'$E^{miss}_{T_{y}}$',
        'root': '#font[52]{E}^{miss}_{T_{y}}',
        'filename': 'MET_y',
        'binning': (20, -75, 75),
        'scale': 0.001,
        'units': 'GeV',
        'legend': 'left',
    },
    'MET_phi': {
        'title': r'$E^{miss}_{T} \phi$',
        'root': '#font[52]{E}^{miss}_{T} #phi',
        'filename': 'MET_phi',
        'binning': (20, -math.pi, math.pi),
    },
    'dPhi_min_tau_MET': {
        'title': r'min[$\Delta\phi$($\tau$,\/$E^{miss}_{T}$)]',
        'root': '#font[52]{min}[#font[152]{#Delta#phi}(#font[152]{#tau},#font[52]{E}^{miss}_{T})]',
        'filename': 'dPhi_min_tau_MET',
        'binning': (20, 0, math.pi),
    },
    'MET_bisecting': {
        'title': r'$E^{miss}_{T}$ bisects',
        'root': '#font[52]{E}^{miss}_{T} bisects',
        'filename': 'MET_bisecting',
        'binning': (2, -0.5, 1.5),
        'legend': 'left',
        'integer': True,
    },
    #'sphericity': {
    #    'title': r'Sphericity',
    #    'root': '#font[52]{Sphericity}',
    #    'filename': 'sphericity',
    #    'bins': 20,
    #    'range': (0, 1),
    #    'cats': ['2J', 'VBF', '1J', '1J_NONBOOSTED'],
    #},
    #'aplanarity': {
    #    'title': r'Aplanarity',
    #    'root': '#font[52]{Aplanarity}',
    #    'filename': 'aplanarity',
    #    'bins': 20,
    #    'range': (0, .15),
    #    'cats': ['2J', 'VBF', '1J', '1J_NONBOOSTED'],
    #},
    'MET_centrality': {
        'title': r'$E^{miss}_{T}$ Centrality',
        'root': '#font[52]{E}^{miss}_{T} #font[152]{#phi} centrality',
        'filename': 'met_centrality',
        'binning': (20, -math.sqrt(2), math.sqrt(2)),
        'legend': 'left',
    },
    'mass_vis_tau1_tau2': {
        'title': r'$m^{vis}_{\tau\tau}$',
        'root': '#font[52]{m}^{vis}_{#font[152]{#tau}#font[152]{#tau}}',
        'filename': 'mass_vis',
        'binning': {
            'PRESELECTION': (20, 30, 150),
            'REST': (20, 30, 150),
            None: (20, 0, 250)},
        'scale': 0.001,
        'units': 'GeV',
        'blind': (70, 110),
    },
    'mass_collinear_tau1_tau2': {
        'title': r'$m^{col}_{\tau\tau}$',
        'root': '#font[52]{m}^{col}_{#font[152]{#tau}#font[152]{#tau}}',
        'filename': 'mass_collinear',
        'binning': (20, 0, 250),
        'units': 'GeV',
        'scale': 0.001,
        'blind': (100, 150),
    },
    'mass_tau1_tau2_jet1': {
        'title': r'$m^{vis}_{j\tau\tau}$',
        'root': '#font[52]{m}_{#font[52]{j}#font[152]{#tau}#font[152]{#tau}}',
        'filename': 'mass_taus_leading_jet',
        'binning': (20, 0, 800),
        'units': 'GeV',
        'scale': 0.001,
    },
    'jet3_centrality': {
        'title': r'j3 Centrality',
        'root': '#font[52]{j}_{3} #font[152]{#eta} centrality',
        'filename': 'jet3_centrality',
        'binning': (20, 0, 1),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'tau_pt_ratio': {
        'title': r'$\tau_{1} p_{T} / \tau_{2} p_{T}$',
        'root': '#font[52]{p}_{T}(#font[152]{#tau}_{1}) / #font[52]{p}_{T}(#font[152]{#tau}_{2})',
        'filename': 'tau_pt_ratio',
        'binning': (16, 1, 5),
    },
    'tau1_pt': {
        'title': r'$\tau_{1} p_{T}$',
        'root': '#font[152]{#tau}_{1} #font[52]{p}_{T}',
        'filename': 'tau1_pt',
        'binning': {
            2011: {
                'PRESELECTION': (10, 35, 90),
                'REST': (10, 35, 90),
                None: (10, 35, 160)},
            2012: {
                'PRESELECTION': (20, 35, 90),
                'REST': (20, 35, 90),
                None: (20, 35, 160)}},
        'scale': 0.001,
        'units': 'GeV',
    },
    'tau2_pt': {
        'title': r'$\tau_{2} p_{T}$',
        'root': '#font[152]{#tau}_{2} #font[52]{p}_{T}',
        'filename': 'tau2_pt',
        'binning': {
            2011: {
                'PRESELECTION': (10, 25, 60),
                'REST': (10, 25, 60),
                None: (10, 25, 100)},
            2012: {
                'PRESELECTION': (20, 25, 60),
                'REST': (20, 25, 60),
                None: (20, 25, 100)}},
        'scale': 0.001,
        'units': 'GeV',
    },
    'tau1_eta': {
        'title': r'$\tau_{1} \eta$',
        'root': '#font[152]{#tau}_{1} #font[152]{#eta}',
        'filename': 'tau1_eta',
        'binning': (20, -3, 3),
        'legend': 'left',
    },
    'tau2_eta': {
        'title': r'$\tau_{2} \eta$',
        'root': '#font[152]{#tau}_{2} #font[152]{#eta}',
        'filename': 'tau2_eta',
        'binning': (20, -3, 3),
        'legend': 'left',
    },
    'tau1_numTrack': {
        'title': r'$\tau_{1}$ Number of Tracks',
        'root': '#font[152]{#tau}_{1} #font[52]{Tracks}',
        'filename': 'tau1_numTrack',
        'binning': (5, -.5, 4.5),
    },
    'tau2_numTrack': {
        'title': r'$\tau_{2}$ Number of Tracks',
        'root': '#font[152]{#tau}_{2} #font[52]{Tracks}',
        'filename': 'tau2_numTrack',
        'binning': (5, -.5, 4.5),
        'integer': True,
    },
    'tau1_numTrack_recounted': {
        'title': r'$\tau_{1}$ Number of Recounted Tracks',
        'root': '#font[152]{#tau}_{1} #font[52]{Recounted Tracks}',
        'filename': 'tau1_numTrack_recounted',
        'binning': (5, 0.5, 5.5),
        'integer': True,
    },
    'tau2_numTrack_recounted': {
        'title': r'$\tau_{2}$ Number of Recounted Tracks',
        'root': '#font[152]{#tau}_{2} #font[52]{Recounted Tracks}',
        'filename': 'tau2_numTrack_recounted',
        'binning': (5, 0.5, 5.5),
        'integer': True,
    },
    #'tau1_nPi0': {
    #    'title': r'$\tau_{1}$ Number of $\pi^0$s',
    #    'root': '#font[152]{#tau}_{1} #font[52]{Number of} #font[152]{#pi}^{0}#font[52]{s}',
    #    'filename': 'tau1_npi0',
    #    'bins': 7,
    #    'range': (-.5, 6.5),
    #    'integer': True,
    #},
    #'tau2_nPi0': {
    #    'title': r'$\tau_{2}$ Number of $\pi^0$s',
    #    'root': '#font[152]{#tau}_{2} #font[52]{Number of} #font[152]{#pi}^{0}#font[52]{s}',
    #    'filename': 'tau2_npi0',
    #    'bins': 7,
    #    'range': (-.5, 6.5),
    #    'integer': True,
    #},
    #'tau_x_product': {
    #    'title': r'Product of $\tau$ Visible Momentum Fractions',
    #    'root': 'Product of #font[152]{#tau} #font[52]{Visible Momentum Fractions}',
    #    'filename': 'tau_x_product',
    #    'bins': 20,
    #    'range': (-9, 16),
    #},
    #'tau_x_sum': {
    #    'title': r'Sum of $\tau$ Visible Momentum Fractions',
    #    'root': 'Sum of #font[152]{#tau} #font[52]{Visible Momentum Fractions}',
    #    'filename': 'tau_x_sum',
    #    'bins': 20,
    #    'range': (-6, 8),
    #},
    #'tau1_collinear_momentum_fraction': {
    #    'title': r'$\tau_{1}$ Visible Momentum Fraction',
    #    'root': '#font[152]{#tau}_{x1}',
    #    'filename': 'tau1_x',
    #    'binning': (20, -3, 4),
    #},
    #'tau2_collinear_momentum_fraction': {
    #    'title': r'$\tau_{2}$ Visible Momentum Fraction',
    #    'root': '#font[152]{#tau}_{x2}',
    #    'filename': 'tau2_x',
    #    'binning': (20, -3, 4),
    #},
    #'tau1_jvtxf': {
    #    'title': r'$\tau_{1}$ JVF',
    #    'root': '#font[152]{#tau}_{1} #font[52]{JVF}',
    #    'filename': 'tau1_jvf',
    #    'bins': 20,
    #    'range': (0, 1),
    #},
    #'tau2_jvtxf': {
    #    'title': r'$\tau_{2}$ JVF',
    #    'root': '#font[152]{#tau}_{2} #font[52]{JVF}',
    #    'filename': 'tau2_jvf',
    #    'bins': 20,
    #    'range': (0, 1),
    #},
    'tau1_BDTJetScore': {
       'title': r'$\tau_{1}$ BDT Score',
       'root': '#font[152]{#tau}_{1} #font[52]{BDT Score}',
       'filename': 'tau1_BDTJetScore',
       'binning': (20, 0.0001, 1.0001),
    },
    'tau2_BDTJetScore': {
       'title': r'$\tau_{2}$ BDT Score',
       'root': '#font[152]{#tau}_{2} #font[52]{BDT Score}',
       'filename': 'tau2_BDTJetScore',
       'binning': (20, 0.0001, 1.0001),
    },
    #'tau1_vertex_prob': {
    #    'title': r'$\tau_{1}$ Primary Vertex Quality',
    #    'root': '#tau_{1} Primary Vertex Quality',
    #    'filename': 'tau1_vertex_quality',
    #    'bins': 20,
    #    'range': (-0.0001, 1.0001),
    #    'cats': ['1J', '2J',]
    #},
    'cos_theta_tau1_tau2': {
        'title': r'$\cos[\alpha(\tau,\tau)]$',
        'root': '#font[52]{cos}(#font[152]{#alpha}_{#font[152]{#tau}#font[152]{#tau}})',
        'filename': 'cos_theta_tau1_tau2',
        'binning': (20, -1, 1),
    },
    'theta_tau1_tau2': {
        'title': r'$\alpha(\tau,\tau)$',
        'root': '#font[152]{#alpha}_{#font[152]{#tau}#font[152]{#tau}}',
        'filename': 'theta_tau1_tau2',
        'binning': (20, 0, math.pi),
    },
    'dR_tau1_tau2': {
        'title': r'$\Delta R(\tau,\tau)$',
        'root': '#font[152]{#Delta}#font[52]{R}(#font[152]{#tau},#font[152]{#tau})',
        'filename': 'dr_tau1_tau2',
        'binning': (8, 0.8, 2.4),
        'ypadding': (0.5, 0),
    },
    'dPhi_tau1_tau2': {
        'title': r'$\Delta \phi(\tau,\tau)$',
        'root': '#font[152]{#Delta#phi}(#font[152]{#tau},#font[152]{#tau})',
        'filename': 'dphi_tau1_tau2',
        'binning': (12, 0., 2.4),
        'legend': 'left',
    },
    'dEta_tau1_tau2': {
        'title': r'$\Delta \eta(\tau,\tau)$',
        'root': '#font[152]{#Delta#eta}(#font[152]{#tau},#font[152]{#tau})',
        'filename': 'deta_tau1_tau2',
        'binning': {
            'BOOSTED': (10, 0, 1.5),
            'VBF': (10, 0, 1.5),
            'REST': (10, 0, 1.5),
            None: (10, 0, 2.5)},
        'ypadding': (0.5, 0),
        'legend': 'right',
    },
    #'tau1_charge': {
    #    'title': r'$\tau_1$ Charge',
    #    'root': '#font[152]{#tau}_{1} #font[52]{Charge}',
    #    'filename': 'tau1_charge',
    #    'bins': 5,
    #    'range': (-2.5, 2.5),
    #    'integer': True,
    #},
    #'tau2_charge': {
    #    'title': r'$\tau_2$ Charge',
    #    'root': '#font[152]{#tau}_{2} #font[52]{Charge}',
    #    'filename': 'tau2_charge',
    #    'bins': 5,
    #    'range': (-2.5, 2.5),
    #    'integer': True,
    #},
    #'tau1_seedCalo_centFrac': {
    #    'title': r'$\tau_1$ Centrality Fraction',
    #    'root': '#font[152]{#tau}_{1} #font[52]{Centrality Fraction}',
    #    'filename': 'tau1_centfrac',
    #    'bins': 20,
    #    'range': (0.5, 1),
    #},
    #'tau_centrality_product': {
    #    'title': r'$\tau$ Centrality Product',
    #    'root': '#font[152]{#tau} centrality product',
    #    'filename': 'tau_centrality_product',
    #    'bins': 20,
    #    'range': (0, 1),
    #    'cats': ['2J', 'VBF']
    #},
    'tau1_centrality': {
        'title': r'$\tau_1$ Centrality',
        'root': '#font[152]{#tau}_{1} #font[152]{#eta} centrality',
        'filename': 'tau1_centrality',
        'binning': (20, 0, 1),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'tau2_centrality': {
        'title': r'$\tau_2$ Centrality',
        'root': '#font[152]{#tau}_{2} #font[152]{#eta} centrality',
        'filename': 'tau2_centrality',
        'binning': (20, 0, 1),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'jet1_eta': {
        'title': r'jet$_{1}$ $\eta$',
        'root': '#font[152]{#eta}(#font[52]{j}1)',
        'filename': 'jet1_eta',
        'binning': (20, -5, 5),
        'cats': ['2J', 'VBF', '1J', '1J_NONBOOSTED'],
        'legend': 'left',
    },
    'jet2_eta': {
        'title': r'jet$_{2}$ $\eta$',
        'root': '#font[152]{#eta}(#font[52]{j}2)',
        'filename': 'jet2_eta',
        'binning': (20, -5, 5),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'jet1_phi': {
        'title': r'jet$_{1}$ $\phi$',
        'root': '#font[152]{#phi}(#font[52]{j}1)',
        'filename': 'jet1_phi',
        'binning': (20, -5, 5),
        'cats': ['2J', 'VBF', '1J', '1J_NONBOOSTED'],
        'legend': 'left',
    },
    'jet2_phi': {
        'title': r'jet$_{2}$ $\phi$',
        'root': '#font[152]{#phi}(#font[52]{j}2)',
        'filename': 'jet2_phi',
        'binning': (20, -5, 5),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'o1': {
        'title': r'0$_{1}$',
        'root': '#font[152]{01}(#font[52]{O}1)',
        'filename': 'o1',
        'binning': (6, -5, 5),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    'jet1_pt': {
        'title': r'jet$_{1}$ $p_{T}$',
        'root': '#font[52]{p}_{T}(#font[52]{j}1)',
        'filename': 'jet1_pt',
        'binning': (20, 20, 200),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['2J', 'VBF', '1J', '1J_NONBOOSTED']
    },
    'jet2_pt': {
        'title': r'jet$_{2}$ $p_{T}$',
        'root': '#font[52]{p}_{T}(#font[52]{j}2)',
        'filename': 'jet2_pt',
        'binning': (20, 20, 200),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['2J', 'VBF']
    },
    'dEta_jets': {
        'title': r'$\Delta\eta(jet_{1},\/jet_{2})$',
        'root': '#font[152]{#Delta#eta}(#font[52]{j}_{1},#font[52]{j}_{2})',
        'filename': 'dEta_jets',
        'cuts': 'dEta_jets > 0', # ignore default value in plot
        'binning': {
            'VBF': (10, 2, 7),
            None: (14, 0, 7)},
        'cats': ['2J', 'VBF', 'PRESELECTION']
    },
    'eta_product_jets': {
        'title': r'jet$_{1}$ $\eta \times \/$ jet$_{2}$ $\eta$',
        'root': '#font[152]{#eta}_{#font[52]{j}_{1}} #times #font[152]{#eta}_{#font[52]{j}_{2}}',
        'filename': 'eta_product_jets',
        'binning': (15, -10, 5),
        'cats': ['2J', 'VBF'],
        'legend': 'left',
    },
    #'eta_product_jets_boosted': {
    #    'title': r'Boosted $\eta_{jet_{1}} \times \/ \eta_{jet_{2}}$',
    #    'root': 'Boosted #eta_{jet_{1}} #times #eta_{jet_{2}}',
    #    'filename': 'eta_product_jets_boosted',
    #    'bins': 20,
    #    'range': (-10, 10),
    #    'cats': ['2J', 'VBF']
    #},
    'mass_jet1_jet2': {
        'title': r'$M(jet_{1},\/jet_{2})$',
        'root': '#font[52]{m}_{#font[52]{j}#font[52]{j}}',
        'filename': 'mass_jet1_jet2',
        'binning': (20, 0, 1000),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['2J', 'VBF']
    },
    'resonance_pt': {
        'title': r'$p_{T}^{H}$',
        'root': '#font[52]{p}_{T}^{H}',
        'filename': 'resonance_pt',
        'binning': {
            'BOOSTED': (18, 70, 250),
            'VBF': (25, 0, 250),
            'REST': (11, 0, 110),
            None: (20, 0, 200)},
        'scale': 0.001,
        'units': 'GeV',
    },
}

from . import MMC_VERSION
mmc = MMC_VERSION

VARIABLES['mmc%d_mass' % mmc] = {
    'title': r'$m^{MMC}_{\tau\tau}$',
    'root': '#font[52]{m}^{MMC}_{#font[152]{#tau}#font[152]{#tau}}',
    'filename': 'mmc%d_mass' % mmc,
    'binning': {
        2011: (25, 0, 250),
        2012: (25, 0, 250)},
    'units': 'GeV',
    'blind': (100, 150),
}

"""
VARIABLES['mmc%d_MET_et' % mmc] = {
    'title': r'$E^{miss}_{T}$ MMC',
    'root': '#font[52]{MMC} #font[52]{E}^{miss}_{T}',
    'filename': 'mmc%d_MET' % mmc,
    'binning': (20, 0, 100),
    'units': 'GeV',
}

VARIABLES['mmc%d_MET_etx' % mmc] = {
    'title': r'MMC $E^{miss}_{T_{x}}$',
    'root': '#font[52]{MMC} #font[52]{E}^{miss}_{T_{x}}',
    'filename': 'mmc%d_MET_x' % mmc,
    'binning': (20, -75, 75),
    'units': 'GeV',
}

VARIABLES['mmc%d_MET_ety' % mmc] = {
    'title': r'MMC $E^{miss}_{T_{y}}$',
    'root': '#font[52]{MMC} #font[52]{E}^{miss}_{T_{y}}',
    'filename': 'mmc%d_MET_y' % mmc,
    'binning': (20, -75, 75),
    'units': 'GeV',
}

VARIABLES['mmc%d_resonance_pt' % mmc] = {
    'title': r'MMC $p_T^H$',
    'root': 'MMC #font[52]{p}_{T}^{H}',
    'filename': 'mmc%d_resonance_pt' % mmc,
    'binning': {
        'BOOSTED': (20, 50, 200),
        None: (20, 0, 200)},
    'units': 'GeV',
}
"""
