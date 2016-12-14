# -*- coding: utf-8 -*-
# 
# @author <ashton@broadinstitute.org>
# ------------------------------------------------

# package imports
# ---------------
import os
import pymongo


# Session
#--------
class Session(object):
    """
    Object for managing connection to internal databases.

    eventually extend to manage Flask server/plotting as well?

    Args:
        host (str): Host with database to connect to.
        port (int): Port to connect to database with.
        write (bool): Whether or not to allow writing to the database.
    """
    def __init__(self, host='localhost', port=27017, write=True, data_home=None):
        self.host = host
        self.port = port
        self.write = write
        self._local_home = data_home if data_home[-1] == '/' else data_home+'/'
        self._client = None
        return

    @property
    def mongodb_client(self):
        if self._client is None:
            try:
                self._client = pymongo.MongoClient(self.host, self.port)
            except pymongo.errors.ServerSelectionTimeoutError:
                raise AssertionError('Could not connect to database! Try using `mongod` to start mongo server.')
        return self._client

    @property
    def client_cohorts(self):
        data_stores = [str(cohort).replace('petridish-', '').rsplit('-',1)[0] for cohort in self.mongodb_client.database_names() if 'petridish' in cohort]
        # assert every collection has a local dir in _local_home
        return sorted(data_stores)

    @property
    def local_cohorts(self):
        data_stores = [cohort for cohort in os.listdir(self._local_home)]
        # assert 'petri-dish' in every collection in data_stores
        return sorted(data_stores)

    @property
    def local_datasets(self):
        data_stores = {x:sorted(os.listdir(self._local_home+x)) for x in self.local_cohorts}
        return data_stores

    @property
    def available_cohorts(self):
        available = set(self.local_cohorts) & set(self.client_cohorts)
        return sorted(list(available))

    @property
    def available_datasets(self):
        return self.local_datasets


# metaclasses
# -----------
class DocRequire(type):
    """
    Metaclass forcing requirement of docstrings on all public
    class methods.
    """

    def __init__(self, name, bases, attrs):
        for key, value in attrs.items():
            if key.startswith("__") or key.startswith("_"):
                continue
            if not hasattr(value, "__call__"):
                continue
            if not getattr(value, '__doc__'):
                raise TypeError("%s must have a docstring" % key)
        type.__init__(self, name, bases, attrs)

