from abc import ABC, abstractmethod
from copy import deepcopy
import sys, os, json, shutil
import hashlib
import numpy as np
import joblib
import logging
log = logging.getLogger(__name__)
import ppdg


def scoring_function(score_config):
    """
    Make the scoring function
    """
    if score_config is None:
        return None
    
    score_config = deepcopy(score_config)
    score_dict = {
        "ppdg" : PPDG,
        "singlepoint" : SinglePoint,
        "ppdx" : PPDX
    }
    return score_dict[score_config]

class ScoringFunction(ABC):
    """
    Abstract class for scoring functions
    """
    @abstractmethod
    def __call__(self, ag, ab):
        raise NotImplementedError

class PPDG(ScoringFunction):

    def __init__(self, ag, ab):
        # Working directory where all models will be created
        ppdg.WRKDIR = os.getcwd()

        # Path to the directory where to find the calc_energy executable and the two potential files
        ppdg.RFSPP = '/projects/kasp8735/software/affinity/'

        # Path where to find the charmm executable
        ppdg.CHARMM = '/projects/kasp8735/software/charmm/exec/gnu/'

    def __call__(self, ag, ab):
        cpxseq = '%s/%s' % (ag, ab)
        name = hashlib.md5(cpxseq.encode()).hexdigest()
        print('Directory name is', name)
        dg = ppdg.eval_pkl('pathogen/hiv/regression-hiv-rfha-veryfast-6224.pkl', '%s/%s' % (ag, ab), (1,2), 'pathogen/hiv/tpl_5fyj.pdb', name=None, nmodels=1, ncores=1)
        # dg = ppdg.eval_pkl('pathogen/hiv/regression-hiv-rfha-veryfast-6224.pkl', '%s/%s' % (ag, ab), (1,2), 'pathogen/hiv/tpl_5fyj.pdb', name=None, nmodels=10, ncores=1)
        print('Binding affinity is', dg)
        os.system(". ./RMDIR.sh")
        os.mkdir('temp_dir')
        os.system('mv %s ./temp_dir' %(name))
        if not os.path.exists('/rc_scratch/jofa3453/hiv/master_dir/%s' %(name)):
        # if not os.path.exists('/pl/active/rdi_data/hiv/master_dir/%s' %(name)):
            if os.path.isfile('./temp_dir/%s/descriptors.json' %(name)):
                # shutil.copytree('./temp_dir/%s' %(name), '/pl/active/rdi_data/hiv/master_dir/%s' %(name))
                shutil.copytree('./temp_dir/%s' %(name), '/rc_scratch/jofa3453/hiv/master_dir/%s' %(name))
        shutil.rmtree('temp_dir')
        return dg

class SinglePoint(ScoringFunction):
    def singlepoint_score(self, ab, ag):
        pass

class PPDX(ScoringFunction):
    def ppdx_score(self, ab, ag):
        pass



        