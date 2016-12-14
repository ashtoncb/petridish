# -*- coding: utf-8 -*-
# 
# @author <ashton@broadinstitute.org>
# ------------------------------------------------

# package imports
# ---------------
import os
from gems import composite

# get Session
# -----------
from petridish.src.utils import Session

# metadata
# --------
__author__ = 'Ashton Berger'
__email__ = 'ashtoncberger@broadinstitute.org'
__version__ = '0.0.1'

# path configuration
# ------------------
__base__ = os.path.dirname(os.path.realpath(__file__))
__resources__ = os.path.join(__base__, 'resources')
__config__ = os.path.join(__base__, 'config.json')
__db__ = os.environ.get('PETRI_DATA', '/data/petri')
__data__ = os.path.join(__db__, 'collections')

if not os.path.exists(__data__):
    os.makedirs(__data__)

# instantiate session
# -------------------
default_dict = {
    'data_home':__data__,
    # 'resource_home':__resources__,
    'host':'localhost', 
    'port':27017, 
    'write':True,
    }
defaults = composite(default_dict)

if os.path.exists(__config__):
    with open(__config__, 'r') as config_reader:
        config = composite(config_reader)
        config_reader.close()
    defaults = defaults + config
    '''
    blake says:
        defaults = composite(defaults)
        config = defaults + config
    '''
global session
session = Session(**defaults)

# get rest of petri imports
# -------------------------
from petridish.src.api import *

