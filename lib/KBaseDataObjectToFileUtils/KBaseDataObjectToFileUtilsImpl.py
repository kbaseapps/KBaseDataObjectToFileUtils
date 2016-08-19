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
    GIT_COMMIT_HASH = "db0fb60fa160af3cced34b7a0e6fe1db33cb2781"
    
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
        if 'nuc_seq' != params or params['nuc_seq'] == None:
            raise ValueError('Method TranslateNucToProtSeq() requires nuc_seq parameter')
        if 'genetic_code' != params or params['genetic_code'] == None:
            params['genetic_code'] = '11'

        if params['genetic_cde'] != '11':
            raise ValueError('Method TranslateNucToProtSeq() only knows genetic code 11')
        
        nuc_seq = nuc_seq.upper()
        prot_seq = ''

        genetic_code = dict()
        genetic_code['11'] = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
            }
        prot_seq = ''.join([genetic_code[table].get(nuc_seq[3*i:3*i+3],'X') for i in range(len(nuc_seq)//3)])

        returnVal = dict()
        returnVal['prot_seq'] = prot_seq
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
           "fasta_file_path" of type "path_type", parameter "feature_ids" of
           list of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeToFASTA
        # init and clean up params
        genome_ref = params['genome_ref']
        file = params['file']
        dir = params['dir']
        console = params['console']
        invalid_msgs = params['invalid_msgs']
        residue_type = params['residue_type']
        feature_type = params['feature_type']
        record_id_pattern = params['record_id_pattern']
        record_desc_pattern = params['record_desc_pattern']
        case = params['case']
        linewrap = params['linewrap']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'nuc'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = '%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_id%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0

        # init and simplify
        feature_ids = []
        feature_sequence_found = False
        residue_type = residue_type[0:3].lower()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch
        fasta_file_path = os.path.join(dir, file)
        self.log(console, 'KB SDK data2file Genome2Fasta: writing fasta file: '+fasta_file_path)

        # get genome object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            genome_object = ws.get_objects([{'ref':genome_ref}])[0]['data']
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()
            

        # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
        #
        #records = []
        self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG

        with open(fasta_file_path, 'w', 0) as fasta_file_handle:
                        
            for feature in genome_object['features']:

                if feature_type == 'ALL' or feature_type == feature['type']:

                    # protein recs
                    if residue_type == 'pro' or residue_type == 'pep':
                        if feature['type'] != 'CDS':
                            continue
                        elif 'protein_translation' not in feature or feature['protein_translation'] == None:
                            self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_id_pattern
                            rec_desc = record_desc_pattern
                            rec_id = record_header_sub(rec_id, feature['id'], genome_object['id'])
                            rec_desc = record_header_sub(rec_desc, feature['id'], genome_object['id'])
                            seq = feature['protein_translation']
                            seq = seq.upper() if case == 'U' else seq.lower()

                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            #record = SeqRecord(Seq(seq), id=rec_id, description=rec_desc)
                            #records.append(record)
                            feature_ids.append(feature['id'])
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if 'dna_sequence' not in feature or feature['dna_sequence'] == None:
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_id_pattern
                            rec_desc = record_desc_pattern
                            rec_id = record_header_sub(rec_id, feature['id'], genome_object['id'])
                            rec_desc = record_header_sub(rec_desc, feature['id'], genome_object['id'])
                            seq = feature['dna_sequence']
                            seq = seq.upper() if case == 'U' else seq.lower()

                            rec_rows = []
                            rec_rows.append('>'+rec_id+' '+rec_desc)
                            if linewrap == None or linewrap == 0:
                                rec_rows.append(seq)
                            else:
                                seq_len = len(seq)
                                base_rows_cnt = seq_len//linewrap
                                for i in range(base_rows_cnt):
                                    rec_rows.append(seq[i*linewrap:(i+1)*linewrap])
                                rec_rows.append(seq[base_rows_cnt*linewrap:])
                            rec = "\n".join(rec_rows)+"\n"

                            #record = SeqRecord(Seq(seq), id=rec_id, description=rec_desc)
                            #records.append(record)
                            feature_ids.append(feature['id'])
                            fasta_file_handle.write(rec)

        # report if no features found
        if not feature_sequence_found:
            self.log(invalid_msgs, "No sequence records found in Genome "+genome_object['id']+" of residue_type: "+residue_type+", feature_type: "+feature_type)
        #else:
        #    SeqIO.write(records, fasta_file_path, "fasta")


        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path'] = fasta_file_path
        returnVal['feature_ids'] = feature_ids
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
           "fasta_file_path" of type "path_type", parameter "feature_ids" of
           list of String
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
