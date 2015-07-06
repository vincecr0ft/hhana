from .common import *
from .mva import *
from .cb import *
from .mva_cb_overlap import *

CATEGORIES = {
    # Preselection
    'presel': [
        Category_Preselection,
        ],
    'presel_deta_controls': [
        Category_Preselection_DEta_Control,
        ],
    
    '1j_inclusive': [
        Category_1J_Inclusive,
        ],
    
    # CB Categories
    'cuts' : [
        Category_Cuts_VBF_LowDR,
        Category_Cuts_VBF_HighDR_Tight,
        Category_Cuts_VBF_HighDR_Loose,
        Category_Cuts_Boosted_Tight,
        Category_Cuts_Boosted_Loose,
        ],
    'cuts_vbf' : [
        Category_Cuts_VBF_LowDR,
        Category_Cuts_VBF_HighDR_Tight,
        Category_Cuts_VBF_HighDR_Loose,
        ],
    'cuts_boosted' : [
        Category_Cuts_Boosted_Tight,
        Category_Cuts_Boosted_Loose,
        ],
    'cuts_2011' : [
        Category_Cuts_VBF_LowDR,
        Category_Cuts_VBF_HighDR,
        Category_Cuts_Boosted_Tight,
        Category_Cuts_Boosted_Loose,
        ],
    'cuts_vbf_2011' : [
        Category_Cuts_VBF_LowDR,
        Category_Cuts_VBF_HighDR,
        ],
    'cuts_boosted_2011' : [
        Category_Cuts_Boosted_Tight,
        Category_Cuts_Boosted_Loose,
        ],
    'cuts_merged' : [
        Category_Cuts_VBF,
        Category_Cuts_Boosted,
        ],
    'cuts_cr' : [
        Category_Cuts_VBF_CR,
        Category_Cuts_Boosted,
        ],
    'cuts_studies' : [
        Category_Cuts_Boosted_Tight,
        Category_Cuts_Boosted_Tight_NoDRCut,
        ],    
    # MVA Categories
    'mva': [
        Category_VBF,
        Category_Boosted,
        ],
    'mva_vbf': [
        Category_VBF,
        ],
    'mva_all': [
        Category_VBF,
        Category_Boosted,
        Category_Rest,
        ],
    'mva_deta_controls': [
        Category_VBF_DEta_Control,
        Category_Boosted_DEta_Control,
        ],
    'mva_workspace_controls': [
        Category_Rest,
        ],

    # CB/MVA Overlap Categories
    'overlap': [
        Category_Cut_VBF_MVA_VBF,
        Category_Cut_VBF_MVA_Boosted,
        Category_Cut_Boosted_MVA_VBF,
        Category_Cut_Boosted_MVA_Boosted,
        ],
    'disjonction': [
        Category_Cut_VBF_Not_MVA,
        Category_Cut_Boosted_Not_MVA,
        Category_MVA_VBF_Not_Cut,
        Category_MVA_Boosted_Not_Cut,
        ],
    'overlap_yields': [
        Category_Cut_VBF_MVA_VBF,
        Category_Cut_VBF_MVA_Boosted,
        Category_Cut_Boosted_MVA_VBF,
        Category_Cut_Boosted_MVA_Boosted,
        Category_Cut_VBF_Not_MVA,
        Category_Cut_Boosted_Not_MVA,
        Category_MVA_VBF_Not_Cut,
        Category_MVA_Boosted_Not_Cut,
        ],

    'overlap_details': [
        Category_Cut_VBF_MVA_VBF,
        Category_Cut_VBF_MVA_Boosted,
        Category_Cut_VBF_MVA_Presel,
        Category_Cut_Boosted_MVA_VBF,
        Category_Cut_Boosted_MVA_Boosted,
        Category_Cut_Boosted_MVA_Presel,
        Category_Cut_Presel_MVA_VBF,
        Category_Cut_Presel_MVA_Boosted,
        Category_Cut_Presel_MVA_Presel,
        Category_Cut_VBF_Not_MVA_VBF,
        Category_Cut_VBF_Not_MVA_Boosted,
        #Category_Cut_VBF_Not_MVA_Presel,
        Category_Cut_Boosted_Not_MVA_VBF,
        Category_Cut_Boosted_Not_MVA_Boosted,
        #Category_Cut_Boosted_Not_MVA_Presel,
        Category_Cut_Presel_Not_MVA_VBF,
        Category_Cut_Presel_Not_MVA_Boosted,
        #Category_Cut_Presel_Not_MVA_Presel,
        Category_MVA_Presel_Not_Cut_VBF,
        Category_MVA_Presel_Not_Cut_Boosted,
        Category_MVA_Presel_Not_Cut_Presel,
        Category_MVA_VBF_Not_Cut_VBF,
        Category_MVA_VBF_Not_Cut_Boosted,
        Category_MVA_VBF_Not_Cut_Presel,
        Category_MVA_Boosted_Not_Cut_VBF,
        Category_MVA_Boosted_Not_Cut_Boosted,
        Category_MVA_Boosted_Not_Cut_Presel,
        ]
}
