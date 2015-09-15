import os

# rootpy imports
from rootpy.tree import Cut
from rootpy.io import root_open

# root_numpy imports
from root_numpy import rec2array, evaluate

# local imports
from . import log
from .sample import SystematicsSample, Background, MC
from ..regions import REGIONS
from .. import DAT_DIR


class Ztautau(Background):
    NORM_BY_THEORY = False

    def histfactory(self, sample, category, systematics=False, **kwargs):
        # isolation systematic
        sample.AddOverallSys(
            'ATLAS_ANA_HH_{0:d}_Isolation'.format(self.year),
            1. - 0.06,
            1. + 0.06)
        if self.workspace_norm is False:
            print 'workspacenorm is false'
            return
        elif self.workspace_norm is not None:
            sample.AddNormFactor(
                'ATLAS_norm_HH_{0:d}_Ztt'.format(self.year),
                self.workspace_norm,
                self.workspace_norm,
                self.workspace_norm,
                True) # const
        elif self.constrain_norm:
            # overallsys
            error = self.scale_error / self.scale
            sample.AddOverallSys(
                'ATLAS_norm_HH_{0:d}_Ztt'.format(self.year),
                1. - error,
                1. + error)
        else:
            sample.AddNormFactor(
                'ATLAS_norm_HH_{0:d}_Ztt'.format(self.year),
                1., 0., 50., False) # floating

    def __init__(self, *args, **kwargs):
        """
        Instead of setting the k factor here
        the normalization is determined by a fit to the data
        """
        self.scale_error = 0.
        self.workspace_norm = kwargs.pop('workspace_norm', None)
        self.constrain_norm = kwargs.pop('constrain_norm', False)
        super(Ztautau, self).__init__(*args, **kwargs)


class MC_Ztautau(Ztautau, MC):
    pass


class MC_Ztautau_DY(MC_Ztautau):
    pass


class Pythia_Ztautau(MC_Ztautau):
    def xsec_kfact_effic(self, isample):
        return 1., 1., 1.


class Embedded_Ztautau(Ztautau, SystematicsSample):

    def __init__(self, *args, **kwargs):
        self.mc_weight = kwargs.pop(
            'mc_weight', True)
        self.posterior_trigger_correction = kwargs.pop(
            'posterior_trigger_correction', True)
        self.embedding_spin_weight = kwargs.pop(
            'embedding_spin_weight', True)
        self.embedding_reco_unfold = kwargs.pop(
            'embedding_reco_unfold', True)
        self.embedding_trigger_weight = kwargs.pop(
            'embedding_trigger_weight', True)
        self.tau_trigger_eff = kwargs.pop(
            'tau_trigger_eff', True)
        super(Embedded_Ztautau, self).__init__(*args, **kwargs)
        with root_open(os.path.join(DAT_DIR, 'embedding_corrections.root')) as file:
            self.trigger_correct = file['ebmc_weight_{0}'.format(self.year % 1000)]
            self.trigger_correct.SetDirectory(0)
        if self.systematics:
            # normalize ISOL and MFS variations to same as nominal
            # at preselection
            from ..categories import Category_Preselection
            nps = [
                ('MFS_UP',),
                ('MFS_DOWN',),
                ('ISOL_UP',),
                ('ISOL_DOWN',)]
            nominal_events = self.events(Category_Preselection)[1].value
            for np in nps:
                np_events = self.events(Category_Preselection,
                                        systematic=np)[1].value
                self.norms[np] =  nominal_events / np_events

    def corrections(self, rec):
        # posterior trigger correction
        if not self.posterior_trigger_correction:
            return
        arr = rec2array(rec[['tau1_pt', 'tau2_pt']])
        weights = evaluate(self.trigger_correct, arr)
        return [weights]

    def systematics_components(self):
        # No FAKERATE for embedding since fakes are data
        return super(Embedded_Ztautau, self).systematics_components() + [
            'MFS',
            'ISOL',]

    def weight_fields(self):
        weights = super(Embedded_Ztautau, self).weight_fields()
        if self.mc_weight:
            weights.append('mc_weight')
        if self.embedding_reco_unfold:
            weights.append('embedding_reco_unfold')
        if self.embedding_trigger_weight:
            weights.append('embedding_trigger_weight')
        if self.embedding_spin_weight:
            weights.append('embedding_spin_weight')
        return weights

    def weight_systematics(self):
        systematics = super(Embedded_Ztautau, self).weight_systematics()
        if self.tau_trigger_eff and self.year == 2011:
            systematics.update({
                'TRIGGER': {
                    'UP': [
                        'tau1_trigger_eff_high',
                        'tau2_trigger_eff_high'],
                    'DOWN': [
                        'tau1_trigger_eff_low',
                        'tau2_trigger_eff_low'],
                    'NOMINAL': [
                        'tau1_trigger_eff',
                        'tau2_trigger_eff']}})
        elif self.tau_trigger_eff:
            systematics.update({
                'TRIGGER': {
                    'UP': [
                        'tau1_trigger_eff_sys_high',
                        'tau2_trigger_eff_sys_high'],
                    'DOWN': [
                        'tau1_trigger_eff_sys_low',
                        'tau2_trigger_eff_sys_low'],
                    'NOMINAL': [
                        'tau1_trigger_eff',
                        'tau2_trigger_eff']},
                'TRIGGER_STAT': {
                    'PERIODA_UP': [
                        'tau1_trigger_eff_stat_scale_PeriodA_high',
                        'tau2_trigger_eff_stat_scale_PeriodA_high'],
                    'PERIODA_DOWN': [
                        'tau1_trigger_eff_stat_scale_PeriodA_low',
                        'tau2_trigger_eff_stat_scale_PeriodA_low'],
                    'PERIODBD_BARREL_UP': [
                        'tau1_trigger_eff_stat_scale_PeriodBD_Barrel_high',
                        'tau2_trigger_eff_stat_scale_PeriodBD_Barrel_high'],
                    'PERIODBD_BARREL_DOWN': [
                        'tau1_trigger_eff_stat_scale_PeriodBD_Barrel_low',
                        'tau2_trigger_eff_stat_scale_PeriodBD_Barrel_low'],
                    'PERIODBD_ENDCAP_UP': [
                        'tau1_trigger_eff_stat_scale_PeriodBD_EndCap_high',
                        'tau2_trigger_eff_stat_scale_PeriodBD_EndCap_high'],
                    'PERIODBD_ENDCAP_DOWN': [
                        'tau1_trigger_eff_stat_scale_PeriodBD_EndCap_low',
                        'tau2_trigger_eff_stat_scale_PeriodBD_EndCap_low'],
                    'PERIODEM_BARREL_UP': [
                        'tau1_trigger_eff_stat_scale_PeriodEM_Barrel_high',
                        'tau2_trigger_eff_stat_scale_PeriodEM_Barrel_high'],
                    'PERIODEM_BARREL_DOWN': [
                        'tau1_trigger_eff_stat_scale_PeriodEM_Barrel_low',
                        'tau2_trigger_eff_stat_scale_PeriodEM_Barrel_low'],
                    'PERIODEM_ENDCAP_UP': [
                        'tau1_trigger_eff_stat_scale_PeriodEM_EndCap_high',
                        'tau2_trigger_eff_stat_scale_PeriodEM_EndCap_high'],
                    'PERIODEM_ENDCAP_DOWN': [
                        'tau1_trigger_eff_stat_scale_PeriodEM_EndCap_low',
                        'tau2_trigger_eff_stat_scale_PeriodEM_EndCap_low'],
                    'NOMINAL': []}})
        return systematics

    def cut_systematics(self):
        systematics = super(Embedded_Ztautau, self).cut_systematics()
        if self.year == 2011:
            return systematics
        # isolation treatment in 2012 is different
        systematics.update({
            'ISOL': { # MUON ISOLATION
                'UP': Cut('(embedding_isolation == 2)'),
                'DOWN': Cut(),
                'NOMINAL': Cut('(embedding_isolation >= 1)')}})
        return systematics

    def xsec_kfact_effic(self, isample):
        return 1., 1., 1.


class MC_Embedded_Ztautau(Embedded_Ztautau):
    pass
