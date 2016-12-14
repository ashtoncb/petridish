# -*- coding: utf-8 -*-
# 
# @author <ashton@broadinstitute.org>
# ------------------------------------------------
'''
#TODO:

# add self._archived feature for collections/dataset/samples/segments to tell user if it requires saving

# delete functions?

# group by histology/other attributes
    - maybe do: DataSet.find_samples({'histology':whatever})

# maybe restructure database organization?
    - people report problems with mongodb once number of databases > 100
    - however, oen database can hold ~12,000+ collections 

# query by transcript, etc.
    - need to get chromosome arm coordinates to do a lot of these

# copy number by segment function
    - cn_by_segment(segments)

# graphs
    - heatmap for each sample
    - heatmap for all samples in dataset

'''
# package imports
# ---------------
#import ensembl.models as models
import json
import os
import pandas as pd
import pymongo
import subprocess
from cached_property import cached_property

# ensembl imports
#----------------
import ensembl.models

# petri imports
#--------------
# from petridish.src.utils import DocRequire
from petridish import session


# base/abstract classes
# ---------------------
class Cohort(object):
    """
    Object for managing individual cohorts.
    """
    def __init__(self, name, write=False, update_upon_start=False):
        self.cohort_name = name
        self._write = write
        self._cohort_path = os.path.join(session._local_home, self.cohort_name)
        self._cohort_mongo_name = '-'.join(['petridish', self.cohort_name])
        self._datadicthidden = {}

        if not os.path.exists(self._cohort_path) and self._write:
            print 'No saved cohort named {} detected - creating from scratch.'.format(self.cohort_name)
            os.makedirs(self._cohort_path)
        else:           
            if update_upon_start:
                print 'Loading {} cohort ...'.format(self.cohort_name)
                self.load_archive()
            print 'Loaded {} cohort'.format(self.cohort_name)
        return

    def __getitem__(self, key):
        if key not in self.datasets:
            raise ValueError('This cohort does not have this dataset!')
        return self._datadict[key]

    @property
    def _cohort_path_contents(self):
        return [f for f in os.listdir(self._cohort_path) if not f.startswith('.')]

    @property
    def _datadict(self):
        if not self._cohort_path_contents:
            return None
        for dtype in self._cohort_path_contents:
            if dtype not in self._datadicthidden.keys():
                ds = DataSet(cohort_name=self.cohort_name,dtype=dtype,write=self._write)
                self._datadicthidden[dtype] = ds
        return self._datadicthidden

    @property
    def datasets(self):
        if not self._datadict:
            return None
        return self._datadict.keys()

        
    def add_dataset(self, dtype, files, force=False, write=False):
        if not self._write:
            raise ValueError('Cohort is set to not allow writing or uploading of new data!')
        if dtype in self.datasets and not force:
            raise Exception('Cohort already has a {} dataset! Use force=True to overwrite.'.format(dtype))
        dtype = str(dtype).lower()
        print 'Creating {} dataset from files...'.format(dtype)
        if dtype == 'absolute':
            assert 'table' in files.keys()
            assert 'segtab' in files.keys()
            assert 'sampleinfo' in files.keys()
            ds = DataSet.from_absolute(cohort_name=self.cohort_name,table=files['table'],segtab=files['segtab'],sampleinfo=files['sampleinfo'],write=True)
        elif dtype == 'gistic':
            # assert 'peaks??' in files.keys()
            # assert 'genes???' idk
            # ds = DataSet.from_gistic()
            pass
        else:
            raise ValueError('Not a valid dataset type!')
        print 'Archiving {} dataset...'.format(dtype)
        ds.update_archive()
        return

    def load_dataset(self, dtype):
        return self._datadict[dtype]

    def update_archive(self):
        """
        Extract archive into data directory.
        """
        if not self._write:
            raise ValueError('Cohort is set to not allow writing!')
        for ds in self._datadict.values():
            print 'Updating {} dataset ...'.format(ds.dataset_name)
            ds.update_archive()
        return

    def load_archive(self):
        """
        Extract archives into data directories.
        """
        for ds in self._datadict.values():
            print 'Loading {} dataset ...'.format(ds._dtype)
            ds.load_archive()
        return

    # def delete_archive(self):
    #     proceed = raw_input("Are you sure? Press 1 to proceed: ")
    #     if proceed:
    #         pass
    #         # 'Deleting cohort {} ...'.format(self.cohort_name)
    #         # subprocess.call('rm -rf ')
    #         # return
    #     else:
    #         print '{} was not deleted.'.format(self.cohort_name)
    #     return


class DataSet(object):
    """

    to do: prohibit manual entry of samplelist and segmentlist? - restrict so that user can only upload from files to prevent errors?e

    Attributes:
        name (str): Data set name
        dtype (str): Type of data
        metadata (str): Dictionary of things about how data set was made. 
            E.g. {'cohort':"CCLE", 'date_generated':"09-22-16", 'ABSOLUTE_arguments':{...}} 
    """
    def __init__(self, cohort_name, dtype, info=None, samplelist=None, segmentlist=None, write=False):
        self.cohort = cohort_name
        self._dtype = dtype
        self.metadata = info
        self._sampledata = samplelist
        self._segmentdata = segmentlist
        self._write = write
        self.dataset_name = '-'.join([self.cohort, self._dtype])
        self._cohort_path = os.path.join(session._local_home, self.cohort)
        self._dataset_path = os.path.join(self._cohort_path, self._dtype)
        self._cohort_mongo_name = '-'.join(['petridish', self.cohort])
        self._dataset_mongo_name = '-'.join([self._cohort_mongo_name, self._dtype])
        self._db = None
        return

    # util-related properties
    @property
    def _needs_archive(self):
        return bool(self._sampledata) + bool(self._segmentdata)

    @property
    def db(self):
        if self._db is None:
            try:
                self._db = session.mongodb_client[self._dataset_mongo_name]
            except pymongo.errors.ServerSelectionTimeoutError:
                self._db = None
                raise AssertionError('Could not connect to database! Try using `mongod` to start mongo server.')
        return self._db
    
    def _dataquery(self, model):
        return self.db[model]

    # data properties
    @cached_property
    def samples(self):
        samps = []
        for sampledict in self._dataquery('samples').find():
            samps.append(Sample(query=False,**sampledict))
        return samps

    @cached_property
    def segments(self):
        '''
        segments
        '''
        segs = []
        for segdict in self._dataquery('segments').find():
            segs.append(Segment(query=False,**segdict))
        return segs

    @property
    def histology_types(self):
        return sorted(list(set([str(sample.info['histology']) for sample in self.samples if sample.info])))

    @property
    def primary_sites(self):
        return sorted(list(set([str(sample.info['primary_site']) for sample in self.samples if sample.info])))

    # methods
    def find_samples(self, by):
        raise AssertionError('Not implemented yet!!')
        return

    def find_segments(self, **kwargs):
        '''
        get segments by category
        '''
        keys = map(lambda x: x.lower(), kwargs.keys())
        if 'gene' in keys:
            symbol = kwargs['gene'].lower()
            gene = ensembl.models.Gene(symbol)
            if not gene.contig or not gene.start or not gene.stop:
                raise ValueError('Unable to provide copy number information for gene "{}".'.format(symbol))
            else:
                segments = [seg for seg in self.segments if seg.chromosome == int(gene.contig.replace('chr', '')) and seg.start <= int(gene.stop) and int(gene.start) <= seg.stop]
                return segments if len(segments) > 0 else None
        else:
            print 'Sorry, other methods have yet to be implemented.'
        return None

    # IO methods 
    @classmethod
    def from_absolute(cls, cohort_name, table, segtab, sampleinfo, info=None, write=False):
        # this whole function is kind of brute force but hey it works lol
        if cohort_name[:4].lower() == 'tcga':
            sampleinfo_cols = ['participant_id', 'sample_name', 'tcga_barcode', 'array_name', 'chemistry_plate', 'tumor_normal', 'tissue_site', 'tumor_type', 'primary_disease'] # make this from your mappings dataframe
            table_cols = ['array','sample','call_status','purity','ploidy','genome_doublings','delta',
                        'coverage_for_80p_power','cancer_DNA_fraction','subclonal_genome_fraction',
                        'tau','total_copy_ratio']
            seg_cols = ['sample','chromosome','start','stop','n_probes','length','seg_sigma','W','total_copy_ratio',
                        'modal_total_cn','expected_total_cn','total_HZ','total_amp','corrected_total_cn','rescaled_total_cn',
                        'biallelic','copy_ratio','hscr_a1','hscr_a2','modal_a1','modal_a2','expected_a1','expected_a2',
                        'subclonal_a1','subclonal_a2','cancer_cell_frac_a1','ccf_ci95_low_a1','ccf_ci95_high_a1','cancer_cell_frac_a2',
                        'ccf_ci95_low_a2','ccf_ci95_high_a2','LOH','HZ','SC_HZ','amp_a1','amp_a2','rescaled_cn_a1','rescaled_cn_a2']

            sampleinfo_df = pd.read_table(sampleinfo, header=0, names=sampleinfo_cols, sep='\t')
            sampleinfodata = sampleinfo_df.set_index('sample_name').transpose().to_dict()
            del sampleinfo_df

            sampledf = pd.read_table(table, header=0, names=table_cols, sep='\t')
            sampledf['cohort'] = cohort_name
            sampledf['dataset'] = 'absolute'
            sampledf['sample_name'] = sampledf['sample']
            sampledata = sampledf.set_index('sample_name').transpose().to_dict()
            del sampledf

            for key in sampledata.keys():
                if key in sampleinfodata.keys(): 
                    sampledata[key]['info'] = sampleinfodata[key]
                    sampledata[key]['name'] = sampleinfodata[key]['tcga_barcode']
            sampledatavals = sampledata.values()

            segmentdf = pd.read_table(segtab, header=0, names=seg_cols, sep='\t')
            segmentdf['cohort'] = cohort_name
            segmentdf['dataset'] = 'absolute'
            segmentdf['interval'] = ['{}:{}-{}'.format(str(x),str(y),str(z)) for x,y,z in zip(segmentdf['chromosome'],segmentdf['start'],segmentdf['stop'])]
            segmentdf['sample_segment'] = ['{}>{}'.format(str(j),str(k)) for j,k in zip(segmentdf['sample'],segmentdf['interval'])]
            segmentdatavals = segmentdf.transpose().to_dict().values()
            del segmentdf

        elif cohort_name[:4].lower() == 'ccle':
            sampleinfo_cols = ['CCLE_name','cell_line_primary_name','cell_line_aliases','gender','primary_site','histology',
                            'hist_subtype1','notes','source','expression_array','snp_array','oncomap','hybrid_capture_sequencing']
            table_cols = ['array','sample','call_status','purity','ploidy','genome_doublings','delta',
                        'coverage_for_80p_power','cancer_DNA_fraction','subclonal_genome_fraction',
                        'tau','total_copy_ratio']
            seg_cols = ['sample','chromosome','start','stop','n_probes','length','seg_sigma','W','total_copy_ratio',
                        'modal_total_cn','expected_total_cn','total_HZ','total_amp','corrected_total_cn','rescaled_total_cn',
                        'biallelic','copy_ratio','hscr_a1','hscr_a2','modal_a1','modal_a2','expected_a1','expected_a2',
                        'subclonal_a1','subclonal_a2','cancer_cell_frac_a1','ccf_ci95_low_a1','ccf_ci95_high_a1','cancer_cell_frac_a2',
                        'ccf_ci95_low_a2','ccf_ci95_high_a2','LOH','HZ','SC_HZ','amp_a1','amp_a2','rescaled_cn_a1','rescaled_cn_a2']
            
            sampleinfo_df = pd.read_table(sampleinfo, header=0, names=sampleinfo_cols, sep='\t')
            sampleinfodata = sampleinfo_df.set_index('CCLE_name').transpose().to_dict()
            del sampleinfo_df
            
            sampledf = pd.read_table(table, header=0, names=table_cols, sep='\t')
            sampledf['cohort'] = cohort_name
            sampledf['dataset'] = 'absolute'
            sampledf['CCLE_name'] = sampledf['sample']
            sampledata = sampledf.set_index('CCLE_name').transpose().to_dict()
            del sampledf

            for key in sampleinfodata.keys():
                if key in sampledata.keys(): 
                    sampledata[key]['info'] = sampleinfodata[key]
                    sampledata[key]['name'] = sampleinfodata[key]['cell_line_primary_name']
            sampledatavals = sampledata.values()

            segmentdf = pd.read_table(segtab, header=0, names=seg_cols, sep='\t')
            segmentdf['cohort'] = cohort_name
            segmentdf['dataset'] = 'absolute'
            segmentdf['interval'] = ['{}:{}-{}'.format(str(x),str(y),str(z)) for x,y,z in zip(segmentdf['chromosome'],segmentdf['start'],segmentdf['stop'])]
            segmentdf['sample_segment'] = ['{}>{}'.format(str(j),str(k)) for j,k in zip(segmentdf['sample'],segmentdf['interval'])]
            segmentdatavals = segmentdf.transpose().to_dict().values()
            del segmentdf

        else:
            raise ValueError('No from_absolute() function implemented for cohort type {} yet!'.format(cohort_name))
        return cls(cohort_name, dtype='absolute', samplelist=sampledatavals, segmentlist=segmentdatavals, info=info, write=write)

    # @classmethod
    # def from_gistic(cls, all_lesions, amp_genes, del_genes, broad, focal):
    #     pass

    def update_archive(self):
        if not self._write:
            raise ValueError('Dataset is set to not allow writing of data!')
        if not self._sampledata or not self._segmentdata:
            raise AttributeError('Cannot create or update archive! Did you provide data for this dataset?')
        Sample._archive(dataset_path=self._dataset_path, data=self._sampledata)
        setattr(self, '_sampledata', None)
        Segment._archive(dataset_path=self._dataset_path, data=self._segmentdata)
        setattr(self, '_segmentdata', None)
        self.load_archive()
        return

    def load_archive(self):
        Sample._update(dataquery=self._dataquery, dataset_path=self._dataset_path, dtype=self._dtype)
        Segment._update(dataquery=self._dataquery, dataset_path=self._dataset_path, dtype=self._dtype)
        return

    # def delete_archive(self, filename):
    #     return


class ModelBase(object):
    """
    Abstract class for performing database operations on models.

    right now things default to __defaults__ always, but maybe allow flexilibity in future?

    """
    # __metaclass__ = DocRequire
    __model__ = 'base'
    __indexes__ = []
    __defaults__ = {}

    def __init__(self, query=True, dataquery=None, *args, **kwargs):
        if self.__model__ == 'base':
            raise NotImplementedError('__model__ property must be specified by subclasses!')
        indexers = [t[0] for t in self.__indexes__]
        kwargs = self._consolidate_args(*args, **kwargs)
        self._data = {}

        # if we're querying data
        if query or dataquery:
            # set up query
            if dataquery:
                querydb = dataquery
            elif 'cohort' in kwargs.keys() and 'dataset' in kwargs.keys():
                mongoname = '-'.join(['petridish', kwargs.get('cohort'), kwargs.get('dataset')])
                try:
                    querydb = session.mongodb_client[mongoname]
                except pymongo.errors.ServerSelectionTimeoutError:
                    raise ValueError('Could not connect to database! Try using `mongod` to start mongo server.')
            else:
                raise ValueError('Need to provide a cohort and dataset or a database client to query for an object!')

            # try to select items to query by
            if '_id' in kwargs.keys():
                self._data = self._query(dataquery=querydb, **kwargs)
            elif indexers and any([kwargs.has_key(idxr) for idxr in indexers]):
                to_use = {idxr:kwargs[idxr] for idxr in indexers if idxr in kwargs.keys()}
                self._data = self._query(dataquery=querydb, **to_use)
            else:
                self._data = self._query(dataquery=querydb, **kwargs)

            # if query wasn't successful
            if not self._data:
                print 'Warning: could not find object with provided attributes.'
                self._data = {}
       
        # if we're not querying
        else:
            for key in self.__defaults__:
                if kwargs.has_key(key):
                    self._data[key] = kwargs[key]

        # extend data by defaults
        for key in self.__defaults__:
            if self._data.get(key) is None:
                self._data[key] = self.__defaults__[key]
    
        # save data as attributes
        for key in self._data.keys():
            setattr(self, key, self._data[key])
        return

    def _query(self, dataquery, **kwargs):
        return dataquery[self.__model__].find_one(kwargs)

    def _consolidate_args(self, *args, **kwargs):
        """
        Used to resolve unspecified args as belonging to kwarg set.
        This is particularly useful for resolving identifiers or names
        on the models below, where a user can simply specify either and
        their input will be properly resolved to instantiate the right
        object.
        """
        if len(args) == 1:
            kwargs['id'] = args[0]
        return kwargs

    @staticmethod
    def _archive(dataset_path, data):
        """
        Archive data in database-ready format.
        
        Args:
            datafile (str): Name of archive file.
        """
        raise NotImplementedError('Download method must be implemented by subclasses!')

    @classmethod
    def _update(self, dataquery, dataset_path, dtype):
        """
        Update database with archived data.
        """
        jsonfile = os.path.join(dataset_path, self.__model__ + '.json')
        if not os.path.exists(jsonfile) or os.stat(jsonfile).st_size == 0:
            print 'No {} dataset archive for {} present.'.format(dtype, self.__model__)
        else:
            print 'Loading {} archive...'.format(self.__model__)
            data = json.load(open(jsonfile, 'r'))
            dataquery(self.__model__).drop()
            dataquery(self.__model__).insert_many(data)
            del data
            if len(self.__indexes__) > 0:
                dataquery(self.__model__).create_index(self.__indexes__)
        return None


class Sample(ModelBase):
    '''
    potential re-work of finding CN values?

        estimate_CN(by=, type=modal)
            by can be gene, chromosome, segment, interval, chromosome arm, snp, etc

            if CN value not in genes or snps or whatever, store it?

        genes = {gene:CN}
        snps = {snp:CN}

    '''
    __model__ = 'samples'
    __indexes__ = [
        ('sample', pymongo.ASCENDING),
        ('name', pymongo.ASCENDING)
    ]
    __defaults__ = {
        'name': None,
        'sample': None,
        'array': None,
        'dataset': None,
        'cohort': None,
        'info': {},
        'call_status': None,
        'purity': None,
        'ploidy': None,
        'genome_doublings': None,
        'coverage_for_80p_power': None,
        'cancer_DNA_fraction': None,
        'subclonal_genome_fraction': None,
        'total_copy_ratio': None,
    }

    @cached_property
    def segments(self):
        '''
        all segments belonging to sample in seg file

        look at getattr for gene and transcript in ensembl?
        '''
        mongoname = '-'.join(['petridish', self.cohort, self.dataset])
        try:
            querydb = session.mongodb_client[mongoname]
        except pymongo.errors.ServerSelectionTimeoutError:
            raise AssertionError('Could not connect to database! Try using `mongod` to start mongo server.')
        segobs = []
        for segdict in querydb['segments'].find({'sample':self.sample}):
            segobs.append(Segment(query=True,**segdict))
        return sorted(segobs, key=lambda x: (x.chromosome))

    def find_segments(self, **kwargs):
        '''
        get segments by category

        right now it only does one arguments at a time but i want to extend return_segs and return that at the end instead of returning in middle of if statement
        for key in keys...
        '''
        keys = map(lambda x: x.lower(), kwargs.keys())

        return_segs = []
        if 'contig' in keys or 'chromosome' in keys:
            pass
        if 'amplified' in keys:
            pass
        if 'deleted' in keys:
            pass
        if 'gene' in keys:
            symbol = kwargs['gene'].lower()
            gene = ensembl.models.Gene(symbol)
            if not gene.contig or not gene.start or not gene.stop:
                raise ValueError('Unable to provide copy number information for gene "{}".'.format(symbol))
            else:
                segments = [seg for seg in self.segments if seg.chromosome == int(gene.contig.replace('chr', '')) and seg.start <= int(gene.stop) and int(gene.start) <= seg.stop]
                return segments if len(segments) > 0 else None
        else:
            print 'Sorry, other methods have yet to be implemented.'
        return None

    def modal_ploidy(self):
        cns = [seg.modal_total_cn for seg in self.segments]
        return max(set(cns), key=cns.count)

    @staticmethod
    def CN_by_segment(segments):
        if segments:
            return {segment.sample_segment:segment.modal_total_cn for segment in segments}
        return None

    def CN_by_position(self, contig, coord, quiet=False):
        segments = [seg for seg in self.segments if seg.chromosome == int(contig) and seg.start <= int(coord) <= seg.stop]
        if segments:
            return self.CN_by_segment(segments)
        if not quiet:
            print 'No segment found for sample on chromosome {} at position {}.'.format(str(contig),str(coord))
        return None

    def CN_by_interval(self, ival, quiet=False):
        contig, coords = ival.split(':')
        contig = contig.lower()
        if contig[:3] == 'chr':
            contig = contig.strip('chr')
        if contig[:6] == '_hschr':
            contig = contig.strip('_hschr').split('_')[0]
        try:
            contig = int(contig)
        except ValueError:
            print 'Unable to locate a segment for contig {}'.format(contig)
            return None
        start, stop = coords.split('-')
        segments = [seg for seg in self.segments if seg.chromosome == int(contig) and seg.start <= int(stop) and int(start) <= seg.stop]
        if segments:
            return self.CN_by_segment(segments)
        if not quiet:
            print 'No segment found for sample that intersects this interval.'
        return None

    # def CN_by_band(self, band):
    #     pass

    # def CN_by_arm(self, arm):
    #     pass

    def CN_by_chromosome(self, contig, quiet=False):
        segments = [seg for seg in self.segments if seg.chromosome == int(contig)]
        if segments:
            return self.CN_by_segment(segments)
        if not quiet:
            print 'No segment found for sample on chromosome {}.'.format(str(contig))
        return None

    def CN_by_gene(self, symbol, use_ensembl=True):
        if use_ensembl:
            gene = ensembl.models.Gene(symbol)
            if not gene.contig or not gene.start or not gene.stop:
                raise ValueError('Unable to provide copy number information for gene "{}".'.format(symbol))
            else:
                ival = '{}:{}-{}'.format(str(gene.contig).strip('chr'),str(gene.start),str(gene.stop))
                return self.CN_by_interval(ival, quiet=True)
        print 'Sorry, other gene querying methods have yet to be implemented.'
        return None

    # def CN_by_variant(self):
        # idk??
        # pass

    # def CN_by_transcript(self, use_ensembl=True):
    #     if use_ensembl:
    #         # get gene interval
    #         return self.CN_by_interval(ival)
    #     else:
    #         print 'Sorry, other transcript querying methods have yet to be implemented.'
    #     return

    @staticmethod
    def _archive(dataset_path, data):
        """
        Update table data and relationships.
        """
        datafile = os.path.join(dataset_path, 'samples.json')
        datadir = os.path.dirname(datafile)

        if not os.path.exists(datadir):
            os.makedirs(datadir)

        if not data:
            raise AttributeError('No samples to archive!')

        # save final data structure
        with open(datafile, 'w') as outfile:
            json.dump(data, outfile, indent=2)
        return


class Segment(ModelBase):
    __model__ = 'segments'
    __indexes__ = [
        ('sample', pymongo.ASCENDING),
        ('chromosome', pymongo.ASCENDING)
    ]
    __defaults__ = {
        'sample': None,
        'dataset': None,
        'cohort': None,
        'chromosome': None,
        'start': None,
        'stop': None,
        'interval': None,
        'sample_interval': None,
        'expected_total_cn': None,
        'rescaled_total_cn': None,
        'amp_a1': None,
        'amp_a2': None,
        'HZ': None,
        'SC_HZ': None,
        'LOH': None,
        'total_amp': None,
        'total_HZ': None,
        'total_copy_ratio': None,
        'modal_total_cn': None
    }

    @property
    def info(self):
        return {'cohort': self.cohort, 'dataset': self.dataset, 'sample': self.sample, 'interval': self.interval}

    def print_copy_number(self, rounded=True):
        print 'Sample: {}'.format(self.sample)
        print 'Segment interval: {}'.format(self.interval)
        if not rounded:
            print 'Expected total CN: {}'.format(self.expected_total_cn)
            return self.expected_total_cn
        print 'Modal total CN: {}'.format(self.modal_total_cn)
        return self.modal_total_cn

    @staticmethod
    def _archive(dataset_path, data):
        """
        Update table data and relationships.
        """
        datafile = os.path.join(dataset_path, 'segments.json')
        datadir = os.path.dirname(datafile)
        
        if not os.path.exists(datadir):
            os.makedirs(datadir)

        if not data:
            raise AttributeError('No segments to archive!')

        # save final data structure
        with open(datafile, 'w') as outfile:
            json.dump(data, outfile, indent=2)
        return

