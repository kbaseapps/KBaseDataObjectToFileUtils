# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
import numpy as np
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from biokbase.workspace.client import Workspace as workspaceService
from requests_toolbelt import MultipartEncoder
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

# silence whining
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class KBaseDataObjectToFileUtils:
    '''
    Module Name:
    KBaseDataObjectToFileUtils

    Module Description:
    ** A KBase module: kb_blast
**
** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
**
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/KBaseDataObjectToFileUtils.git"
    GIT_COMMIT_HASH = "3c859e1531e002b5fc9f4483bea68b8019efc4fc"
    
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.scratch = os.path.abspath(config['scratch'])
        # HACK!! temporary hack for issue where megahit fails on mac because of silent named pipe error
        #self.host_scratch = self.scratch
        self.scratch = os.path.join('/kb','module','local_scratch')
        # end hack
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass
    

    def TranslateNucToProtSeq(self, ctx, params):
        """
        Methods for converting KBase Data Objects to common bioinformatics format files
        **
        :param params: instance of type "TranslateNucToProtSeq_Params"
           (TranslateNucToProtSeq() Params) -> structure: parameter "nuc_seq"
           of String, parameter "genetic_code" of String
        :returns: instance of type "TranslateNucToProtSeq_Output"
           (TranslateNucToProtSeq() Output) -> structure: parameter
           "prot_seq" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN TranslateNucToProtSeq
        #END TranslateNucToProtSeq

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method TranslateNucToProtSeq return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeAnnotationToFASTA_Params"
           (GenomeAnnotationToFASTA() Params) -> structure: parameter
           "genome_ref" of type "data_obj_ref", parameter "file" of type
           "path_type", parameter "dir" of type "path_type", parameter
           "console" of list of type "log_msg", parameter "invalid_msgs" of
           list of type "log_msg", parameter "residue_type" of String,
           parameter "feature_type" of String, parameter "record_id_pattern"
           of type "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type "GenomeAnnotationToFASTA_Output"
           (GenomeAnnotationToFASTA() Output) -> structure: parameter
           "output_file" of type "path_type", parameter "feature_ids" of list
           of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeToFASTA
        #END GenomeToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeAnnotationToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeAnnotationToFASTA_Params"
           (GenomeAnnotationToFASTA() Params) -> structure: parameter
           "genome_ref" of type "data_obj_ref", parameter "file" of type
           "path_type", parameter "dir" of type "path_type", parameter
           "console" of list of type "log_msg", parameter "invalid_msgs" of
           list of type "log_msg", parameter "residue_type" of String,
           parameter "feature_type" of String, parameter "record_id_pattern"
           of type "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type "GenomeAnnotationToFASTA_Output"
           (GenomeAnnotationToFASTA() Output) -> structure: parameter
           "output_file" of type "path_type", parameter "feature_ids" of list
           of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeAnnotationToFASTA
        #END GenomeAnnotationToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeAnnotationToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
