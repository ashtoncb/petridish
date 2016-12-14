# -*- coding: utf-8 -*-
# 
# @author <ashton@broadinstitute.org>
# ------------------------------------------------
'''

only things i want to explicitly allow through api:
    petridish.load_cohort
    petridish.load_dataset
    petridish.find_samples

'''
# package imports
# ---------------
# import subprocess

# petri imports
# -------------
from petridish import session
from petridish.src.models import Cohort, DataSet, Sample, Segment
# from petridish.src.utils import *
# from petridish.src.models import *


# housekeeping functions
#-----------------------
def cohorts():
    return session.available_cohorts

def data_summary():
    available = session.available_datasets
    print '{0:<5}{1:<15}{2}'.format('','cohort', 'datasets')
    for n,pair in enumerate(sorted(available.items()), 1):
        print '{0:<5}{1:<15}{2}'.format(n,pair[0],', '.join(pair[1]))
    return 

# def set_options(**kwargs):
#     pass

# def find_samples():
#     return

# def find_segments():
#     return

def load_cohort(name):
    allowed = name in session.available_cohorts
    if not allowed:
        raise ValueError('This cohort is not available.')
    return Cohort(name)

def load_dataset(cohort, dataset):
    allowed = cohort in session.available_cohorts
    if not allowed:
        raise ValueError('This cohort is not available.')
    cort = Cohort(cohort)
    if dataset not in cort.datasets:
        raise ValueError('This cohort does not have this dataset.')
    return cort[dataset]

# def new_cohort(name):
#     return Cohort(name=name)

# def remove_cohort(name):
#     if name in session.client_cohorts:
#         mongoname = name
#         localname = name.replace('petridish','')
#     elif name in session.local_cohorts:
#         mongoname = 'petridish' + name
#         localname = name
#     else:
#         raise 'This is not a real cohort.'
#     proceed = raw_input("Are you sure? Press 1 to proceed: ")
#     if proceed:
#         session.mongodb_client.drop_database(mongoname)
#         subprocess.call('rm -rf <path>')
#         return
#     else:
#         print '{} was not deleted.'.format(name)
#     return

