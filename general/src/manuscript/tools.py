import os

import numpy as np
import pandas as pd

from manuscript import inout

import sys

sys.path.append(
    inout.gustav_path()
)


from gustav import nlm



def get_mesh_children(ui_of_interest):
    """
    Will identify the MeSH ui (e..g: D000001)
    that are childen of provided MeSH
    
    Input:
    ui    str e.g.: D000001
    
    Output
    uis   set with childen ui of input ui
    
    """
    
    ui = nlm.mesh('ui2mn')
    branch = ui[ui['ui']==ui_of_interest]['mn'].values

    if len(branch) > 1:
        pattern = r'^'+r'|^'.join(list(branch))
        children = set(ui[ui['mn'].str.contains(pattern, case=False)]['ui'])

    else:
        pattern = branch[0]
        children = set(ui[ui['mn'].str.startswith(pattern)]['ui'])
        
    return children