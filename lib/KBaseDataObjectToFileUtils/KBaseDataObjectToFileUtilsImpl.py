# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import requests
import json
import re
import traceback
import uuid
from datetime import datetime
from pprint import pprint, pformat
#import numpy as np
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

#from biokbase.workspace.client import Workspace as workspaceService
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.DataFileUtilClient import DataFileUtil

from requests_toolbelt import MultipartEncoder
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
#from doekbase.data_api.annotation.genome_annotation.api import GenomeAnnotationAPI as GenomeAnnotationAPI

# silence whining
import requests
requests.packages.urllib3.disable_warnings()

#END_HEADER


class KBaseDataObjectToFileUtils:
    '''
    Module Name:
    KBaseDataObjectToFileUtils

    Module Description:
    ** A KBase module: KBaseDataObjectToFileUtils
**
** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
**
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.1.0"
    GIT_URL = "https://github.com/kbaseapps/KBaseDataObjectToFileUtils"
    GIT_COMMIT_HASH = "def8d74204a95d1fe9a6c974aa9c28bee4c6b231"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None
    callbackURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    #def _get_ama_features_as_json (self, features_handle_ref, gff_handle_ref, protein_handle_ref):
    def _get_ama_features_as_json (self, features_handle_ref):
        this_id = str(uuid.uuid4())
        this_scratch_dir = os.path.join (self.scratch, this_id)
        json_features_file_path = os.path.join (this_scratch_dir, 'features.json')
        #gff_file_path = os.path.join (this_scratch_dir, 'genes.gff')
        #protein_file_path = os.path.join (this_scratch_dir, 'protein.fasta')

        try:
            dfu = DataFileUtil (self.callbackURL)
        except Exception as e:
            raise ValueError('Unable to connect to DFU: ' + str(e))

        try:
            dfu.shock_to_file({'handle_id': features_handle_ref,
                               'file_path': json_features_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly features from SHOCK: ' + str(e))
        """
        try:
            dfu.shock_to_file({'handle_id': gff_handle_ref,
                               'file_path': gff_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly gffs from SHOCK: ' + str(e))

        try:
            dfu.shock_to_file({'handle_id': protein_handle_ref,
                               'file_path': protein_file_path+'.gz',
                               'unpack': 'uncompress'
                           })
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly protein FASTA from SHOCK: ' + str(e))
        """

        # DEBUG
        """
        print ("SCRATCH CONTENTS")
        sys.stdout.flush()
        for this_file in os.listdir (this_scratch_dir):
            print ("\t"+this_file)
            sys.stdout.flush()

        buf = []
        #with open(json_features_file_path, 'r') as f:
        with open(protein_file_path, 'r') as f:
            for line in f.readlines():
                buf.append (line)
            #features_json = json.load(f)

        print ("FEATURES_JSON:\n"+"\n".join(buf))
        sys.stdout.flush()
        """
        
        with open(json_features_file_path, 'r') as f:
            features_json = json.load(f)


        os.remove(json_features_file_path+'.gz')
        os.remove(json_features_file_path)
        #os.remove(gff_file_path+'.gz')
        #os.remove(gff_file_path)
        #os.remove(protein_file_path+'.gz')
        #os.remove(protein_file_path)

        return features_json


    def _get_features_from_AnnotatedMetagenomeAssembly(self, ctx, ama_ref):

        # get ama object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            ama_object = ws.get_objects2({'objects':[{'ref':ama_ref}]})['data'][0]
            ama_object_data = ama_object['data']
            ama_object_info = ama_object['info']
        except Exception as e:
            raise ValueError('Unable to fetch AnnotatedMetagenomeAssembly object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        # get features from json
        features_handle_ref = ama_object_data['features_handle_ref']
        #gff_handle_ref = ama_object_data['gff_handle_ref']
        #protein_handle_ref = ama_object_data['protein_handle_ref']
        #features_json = self._get_ama_features_as_json (features_handle_ref, gff_handle_ref, protein_handle_ref)
        features_json = self._get_ama_features_as_json (features_handle_ref)

        return features_json


    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.scratch = os.path.abspath(config['scratch'])
        #self.scratch = os.path.join('/kb','module','local_scratch')
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
        if 'nuc_seq' not in params or params['nuc_seq'] == None or params['nuc_seq'] == '':
            raise ValueError('Method TranslateNucToProtSeq() requires nuc_seq parameter')
        if 'genetic_code' not in params or params['genetic_code'] == None or params['genetic_code'] == '':
            params['genetic_code'] = '11'

        if params['genetic_code'] != '11':
            raise ValueError('Method TranslateNucToProtSeq() only knows genetic code 11')
        
        nuc_seq = params['nuc_seq'].upper()
        prot_seq = ''

        genetic_code = params['genetic_code']
        genetic_code_table = dict()
        genetic_code_table['11'] = {
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
        if genetic_code not in genetic_code_table:
            raise ValueError ("genetic code '"+str(genetic_code)+"' not configured in genetic_code_table")

        prot_seq = ''.join([genetic_code_table[genetic_code].get(nuc_seq[3*i:3*i+3],'X') for i in range(len(nuc_seq)//3)])
        if prot_seq.endswith('_'):
            prot_seq = prot_seq.rstrip('_')

        returnVal = dict()
        returnVal['prot_seq'] = prot_seq
        #END TranslateNucToProtSeq

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method TranslateNucToProtSeq return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def ParseFastaStr(self, ctx, params):
        """
        :param params: instance of type "ParseFastaStr_Params"
           (ParseFastaStr() Params) -> structure: parameter "fasta_str" of
           String, parameter "residue_type" of String, parameter "case" of
           String, parameter "console" of type "log_msg", parameter
           "invalid_msgs" of type "log_msg"
        :returns: instance of type "ParseFastaStr_Output" (ParseFastaStr()
           Output) -> structure: parameter "id" of String, parameter "desc"
           of String, parameter "seq" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN ParseFastaStr

        # init
        if 'fasta_str' not in params or params['fasta_str'] == None or params['fasta_str'] == '':
            raise ValueError('Method ParseFastaStr() requires fasta_str parameter')
        input_sequence_buf = params['fasta_str']
        residue_type       = params['residue_type']
        case               = params['case']
        console            = params['console']
        invalid_msgs       = params['invalid_msgs']

        now = int(100*(datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds())
        header_id = 'id.'+str(now)
        header_desc = 'desc.'+str(now)
        if residue_type == None:
            residue_type = 'NUC'
        if case == None:
            case = 'UPPER'

        residue_type = residue_type[0:1].upper()
        case = case[0:1].upper()

        PROT_pattern  = re.compile("^[acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWYxX ]+$")
        DNA_pattern   = re.compile("^[acgtuACGTUnryNRY ]+$")
        space_pattern = re.compile("^[ \t]*$")

        
        # money rock on
        #
        sequence_str_buf = ''
        fastq_format = False
        input_sequence_buf = input_sequence_buf.strip()
        if input_sequence_buf.startswith('@'):
            fastq_format = True
        input_sequence_buf = re.sub ('&apos;', "'", input_sequence_buf)
        input_sequence_buf = re.sub ('&quot;', '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#39;',  "'", input_sequence_buf)
#        input_sequence_buf = re.sub ('&#34;',  '"', input_sequence_buf)
#        input_sequence_buf = re.sub ('&lt;;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#60;',  '<', input_sequence_buf)
#        input_sequence_buf = re.sub ('&gt;',   '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#62;',  '>', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#36;',  '$', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#37;',  '%', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#47;',  '/', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#63;',  '?', input_sequence_buf)
##        input_sequence_buf = re.sub ('&#92;',  chr(92), input_sequence_buf)  # FIX LATER
#        input_sequence_buf = re.sub ('&#96;',  '`', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#124;', '|', input_sequence_buf)
#        input_sequence_buf = re.sub ('&amp;', '&', input_sequence_buf)
#        input_sequence_buf = re.sub ('&#38;', '&', input_sequence_buf)
#        self.log(console,"INPUT_SEQ AFTER: '''\n"+input_sequence_buf+"\n'''")  # DEBUG

        split_input_sequence_buf = input_sequence_buf.split("\n")

        # no header rows, just sequence
        if not input_sequence_buf.startswith('>') and not input_sequence_buf.startswith('@'):
            for line in split_input_sequence_buf:
                if space_pattern.match(line):
                    continue
                line = re.sub (" ","",line)
                line = re.sub ("\t","",line)
                line = re.sub ("\r","",line)
                if (residue_type == 'N' and not DNA_pattern.match(line)) \
                        or (residue_type == 'P' and not PROT_pattern.match(line)):
                    self.log(invalid_msgs,"BAD record:\n"+line+"\n")
                    break
                sequence_str_buf += line
        else:
            # format checks
            for i,line in enumerate(split_input_sequence_buf):
                if fastq_format and line.startswith('@'):
                    format_ok = True
                    seq_len = len(split_input_sequence_buf[i+1])
                    if not DNA_pattern.match(split_input_sequence_buf[i+1]) \
                            or not seq_len > 0 \
                            or not split_input_sequence_buf[i+2].startswith('+') \
                            or not seq_len == len(split_input_sequence_buf[i+3]):
                        format_ok = False
                    if not format_ok:
                        bad_record = "\n".join([split_input_sequence_buf[i],
                                                split_input_sequence_buf[i+1],
                                                split_input_sequence_buf[i+2],
                                                split_input_sequence_buf[i+3]])
                        self.log(invalid_msgs,"BAD record:\n"+bad_record+"\n")
                        break
                elif line.startswith('>'):
                    continue
                else:
                    if space_pattern.match(line):
                        continue
                    if (residue_type == 'N' and not DNA_pattern.match(line)) \
                            or (residue_type == 'P' and not PROT_pattern.match(line)):
                        self.log(invalid_msgs,"BAD record:\n"+("\n".join(split_input_sequence_buf))+"\n")
                        break

            # store that sucker, removing spaces
            for i,line in enumerate(split_input_sequence_buf):
                if line.startswith('>'):
                    if line.find(" ") < 0 and line.find("\t") < 0:
                        header_id = line
                    elif line.find(" ") < 0:
                        (header_id, header_desc) = line.split("\t", 1)
                    elif line.find("\t") < 0:
                        (header_id, header_desc) = line.split(" ", 1)
                    elif line.find(" ") < line.find("\t"):
                        (header_id, header_desc) = line.split(" ", 1)
                    else:
                        (header_id, header_desc) = line.split("\t", 1)

                    for j in range(i+1,len(split_input_sequence_buf)):
                        if space_pattern.match(split_input_sequence_buf[j]):
                            continue
                        if split_input_sequence_buf[j].startswith('>'):
                            break
                        seq_line = re.sub (" ","",split_input_sequence_buf[j])
                        seq_line = re.sub ("\t","",seq_line)
                        seq_line = re.sub ("\r","",seq_line)
                        seq_line = seq_line.lower()
                        sequence_str_buf += seq_line

                    break  # only want first record
                elif line.startswith('@'):
                    if line.find(" ") < 0 and line.find("\t") < 0:
                        header_id = line
                    elif line.find(" ") < 0:
                        (header_id, header_desc) = line.split("\t", 1)
                    elif line.find("\t") < 0:
                        (header_id, header_desc) = line.split(" ", 1)
                    elif line.find(" ") < line.find("\t"):
                        (header_id, header_desc) = line.split(" ", 1)
                    else:
                        (header_id, header_desc) = line.split("\t", 1)

                    seq_line = re.sub (" ","",split_input_sequence_buf[i+1])
                    seq_line = re.sub ("\t","",seq_line)
                    seq_line = re.sub ("\r","",seq_line)
                    seq_line = seq_line.lower()
                    #qual_line = re.sub (" ","",split_input_sequence_buf[i+3])
                    #qual_line = re.sub ("\t","",qual_line)
                    sequence_str_buf += seq_line
                    break  # only want first record

        # tidy up
        #
        if case == 'U':
            sequence_str_buf = sequence_str_buf.upper()
        else:
            sequence_str_buf = sequence_str_buf.lower()
        if header_id[0:1] == '@':
            header_id = re.sub('^@','',header_id)
        elif header_id[0:1] == '>':
            header_id = re.sub('^>','',header_id)

        # Done
        #
        if sequence_str_buf == '':
            raise ValueError ("No sequence found in fasta_str: '"+params['fasta_str']+"'")
        returnVal = dict()
        returnVal['id']   = header_id
        returnVal['desc'] = header_desc
        returnVal['seq']  = sequence_str_buf
        #END ParseFastaStr

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method ParseFastaStr return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeToFASTA_Params"
           (GenomeToFASTA() Params) -> structure: parameter "genome_ref" of
           type "data_obj_ref", parameter "file" of type "path_type",
           parameter "dir" of type "path_type", parameter "console" of list
           of type "log_msg", parameter "invalid_msgs" of list of type
           "log_msg", parameter "residue_type" of String, parameter
           "feature_type" of String, parameter "record_id_pattern" of type
           "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type "GenomeToFASTA_Output" (GenomeToFASTA()
           Output) -> structure: parameter "fasta_file_path" of type
           "path_type", parameter "feature_ids" of list of type "feature_id",
           parameter "feature_id_to_function" of mapping from type
           "feature_id" to String, parameter "genome_ref_to_sci_name" of
           mapping from type "data_obj_ref" to String, parameter
           "genome_ref_to_obj_name" of mapping from type "data_obj_ref" to
           String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeToFASTA

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11) 

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
        feature_id_to_function = dict()
        genome_ref_to_sci_name = dict()
        genome_ref_to_obj_name = dict()
        feature_sequence_found = False
        residue_type = residue_type[0:1].upper()
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
        self.log(console, 'KB SDK data2file Genome2FASTA(): writing fasta file: '+fasta_file_path)

        # get genome object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #genome_object = ws.get_objects([{'ref':genome_ref}])[0]['data']
            genome_obj_base = ws.get_objects2({'objects':[{'ref':genome_ref}]})['data'][0]
            genome_object = genome_obj_base['data']
            genome_info = genome_obj_base['info']
            genome_ref_to_obj_name[genome_ref] = genome_info[NAME_I]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        # set sci name
        genome_ref_to_sci_name[genome_ref] = genome_object['scientific_name']
        feature_id_to_function[genome_ref] = dict()


        # set features list ('features' only contains CDS features')
        target_features = []
        if feature_type == 'CDS' or residue_type == 'P' or 'non_coding_features' not in genome_object:
            target_features = genome_object['features']
        else:
            target_features = genome_object['features'] + genome_object['non_coding_features']

            
        # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
        #
        #records = []
        self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG

        with open(fasta_file_path, 'w', 0) as fasta_file_handle:
                        
            for feature in target_features:
                
                # set function
                if 'function' in feature:  # Feature <= 2.3 (Genome <= 8.2)
                    feature_id_to_function[genome_ref][feature['id']] = feature['function']
                elif 'functions' in feature:  # Feature >= 3.0 (Genome >= 9.0)
                    feature_id_to_function[genome_ref][feature['id']] = "; ".join(feature['functions'])
                else:
                    feature_id_to_function[genome_ref][feature['id']] = 'N/A'
                
                #if feature_type == 'ALL' or feature_type == feature['type']:
                if True:  # don't want to deal with changing indentation

                    # protein recs
                    if residue_type == 'P':
                        #if feature['type'] != 'CDS':
                        #    continue
                        if 'protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == '':
                            #self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, feature['id'], genome_object['id'])
                            rec_desc = record_header_sub(record_desc_pattern, feature['id'], genome_object['id'])
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
                        if feature_type == 'CDS' and ('protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == ''):
                            continue
                        elif 'dna_sequence' not in feature or feature['dna_sequence'] == None or feature['dna_sequence'] == '':
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, feature['id'], genome_object['id'])
                            rec_desc = record_header_sub(record_desc_pattern, feature['id'], genome_object['id'])
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
        returnVal['feature_id_to_function'] = feature_id_to_function
        returnVal['genome_ref_to_sci_name'] = genome_ref_to_sci_name
        returnVal['genome_ref_to_obj_name'] = genome_ref_to_obj_name
        #END GenomeToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def GenomeSetToFASTA(self, ctx, params):
        """
        :param params: instance of type "GenomeSetToFASTA_Params"
           (GenomeSetToFASTA() Params) -> structure: parameter
           "genomeSet_ref" of type "data_obj_ref", parameter "file" of type
           "path_type", parameter "dir" of type "path_type", parameter
           "console" of list of type "log_msg", parameter "invalid_msgs" of
           list of type "log_msg", parameter "residue_type" of String,
           parameter "feature_type" of String, parameter "record_id_pattern"
           of type "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long, parameter "merge_fasta_files" of type "true_false"
        :returns: instance of type "GenomeSetToFASTA_Output"
           (GenomeSetToFASTA() Output) -> structure: parameter
           "fasta_file_path_list" of list of type "path_type", parameter
           "feature_ids_by_genome_id" of mapping from type "genome_id" to
           list of type "feature_id", parameter "feature_id_to_function" of
           mapping from type "feature_id" to String, parameter
           "genome_ref_to_sci_name" of mapping from type "data_obj_ref" to
           String, parameter "genome_ref_to_obj_name" of mapping from type
           "data_obj_ref" to String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN GenomeSetToFASTA

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11) 

        # init and clean up params
        genomeSet_ref = params['genomeSet_ref']
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
        merge_fasta_files = params['merge_fasta_files']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'NUC'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = 'g:%%genome_ref%%.f:%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0
        if merge_fasta_files == None or merge_fasta_files[0:1].upper() == 'T':
            merge_fasta_files = True
        else:
            merge_fasta_files = False

        # init and simplify
        fasta_file_path_list = []
        feature_ids_by_genome_id = dict()
        feature_id_to_function = dict()
        genome_ref_to_sci_name = dict()
        genome_ref_to_obj_name = dict()
        residue_type = residue_type[0:1].upper()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch

        # get genomeSet object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #genomeSet_object = ws.get_objects([{'ref':genomeSet_ref}])[0]['data']
            genomeSet_object = ws.get_objects2({'objects':[{'ref':genomeSet_ref}]})['data'][0]['data']
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        # iterate through genomeSet members
        genome_ids = genomeSet_object['elements'].keys()
        for genome_i in range(len(genome_ids)):
            feature_sequence_found = False
            genome_id = genome_ids[genome_i]
            feature_ids_by_genome_id[genome_id] = []

            if 'ref' not in genomeSet_object['elements'][genome_id] or \
                    genomeSet_object['elements'][genome_id]['ref'] == None or \
                    genomeSet_object['elements'][genome_id]['ref'] == '':
                raise ValueError('GenomeSetToFASTA() cannot handle GenomeSet objects with embedded genome.  Must be a set of genome references')
                #to get the full stack trace: traceback.format_exc()       
            else:
                genome_ref = genomeSet_object['elements'][genome_id]['ref']


            # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
            #
            #records = []
            this_file = file
            if len(genome_ids) > 1 and not merge_fasta_files:
                this_file = this_file+'.'+genome_id
                this_file.replace('/', '-')
            fasta_file_path = os.path.join(dir, this_file)
            #self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG
            
            if genome_i == 0 or not merge_fasta_files:
                fasta_file_handle = open(fasta_file_path, 'w', 0)
                self.log(console, 'KB SDK data2file GenomeSet2FASTA(): writing fasta file: '+fasta_file_path)
            # DEBUG
            self.log(console, "ADDING GENOME: "+str(genome_i+1)+" of "+str(len(genome_ids))+" "+genome_ids[genome_i])


            # get genome object
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                #genome_object = ws.get_objects([{'ref':genome_ref}])[0]['data']
                genome_obj_base = ws.get_objects2({'objects':[{'ref':genome_ref}]})['data'][0]
                genome_object = genome_obj_base['data']
                genome_info = genome_obj_base['info']
                genome_ref_to_obj_name[genome_ref] = genome_info[NAME_I]
            except Exception as e:
                raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            # set sci name
            genome_ref_to_sci_name[genome_ref] = genome_object['scientific_name']
            feature_id_to_function[genome_ref] = dict()

            # set features list ('features' only contains CDS features')
            target_features = []
            if feature_type == 'CDS' or residue_type == 'P' or 'non_coding_features' not in genome_object:
                target_features = genome_object['features']
            else:
                target_features = genome_object['features'] + genome_object['non_coding_features']


            # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
            #
            #records = []
                        
            for feature in target_features:
                fid = feature['id']

                # set function
                if 'function' in feature:  # Feature <= 2.3 (Genome <= 8.2)
                    feature_id_to_function[genome_ref][fid] = feature['function']
                elif 'functions' in feature:  # Feature >= 3.0 (Genome >= 9.0)
                    feature_id_to_function[genome_ref][fid] = "; ".join(feature['functions'])
                else:
                    feature_id_to_function[genome_ref][fid] = 'N/A'

                #self.log (console, "FID:'"+str(fid)+"' FXN: '"+str(feature_id_to_function[genome_ref][fid]))  # DEBUG

                #if feature_type == 'ALL' or feature_type == feature['type']:
                if True:  # don't want to deal with changing indentation

                    # protein recs
                    if residue_type == 'P':
                        #if feature['type'] != 'CDS':
                        #    continue
                        if 'protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == '':
                            #self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
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
                            feature_ids_by_genome_id[genome_id].append(feature['id'])
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if feature_type == 'CDS' and ('protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == ''):
                            continue
                        elif 'dna_sequence' not in feature or feature['dna_sequence'] == None or feature['dna_sequence'] == '':
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
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
                            feature_ids_by_genome_id[genome_id].append(feature['id'])
                            fasta_file_handle.write(rec)


            if genome_i == (len(genome_ids)-1) or not merge_fasta_files:
                self.log(console,"CLOSING FILE: '"+fasta_file_path+"'")  # DEBUG
                fasta_file_handle.close()
                fasta_file_path_list.append(fasta_file_path)


            # report if no features found
            if not feature_sequence_found:
                self.log(invalid_msgs, "No sequence records found in Genome "+genome_id+" of residue_type: "+residue_type+", feature_type: "+feature_type)


        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path_list'] = fasta_file_path_list
        returnVal['feature_ids_by_genome_id'] = feature_ids_by_genome_id
        returnVal['feature_id_to_function'] = feature_id_to_function
        returnVal['genome_ref_to_sci_name'] = genome_ref_to_sci_name
        returnVal['genome_ref_to_obj_name'] = genome_ref_to_obj_name
        #END GenomeSetToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method GenomeSetToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def SpeciesTreeToFASTA(self, ctx, params):
        """
        :param params: instance of type "SpeciesTreeToFASTA_Params"
           (SpeciesTreeToFASTA() Params) -> structure: parameter "tree_ref"
           of type "data_obj_ref", parameter "file" of type "path_type",
           parameter "dir" of type "path_type", parameter "console" of list
           of type "log_msg", parameter "invalid_msgs" of list of type
           "log_msg", parameter "residue_type" of String, parameter
           "feature_type" of String, parameter "record_id_pattern" of type
           "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long, parameter "merge_fasta_files" of type "true_false"
        :returns: instance of type "SpeciesTreeToFASTA_Output"
           (SpeciesTreeToFASTA() Output) -> structure: parameter
           "fasta_file_path_list" of list of type "path_type", parameter
           "feature_ids_by_genome_id" of mapping from type "genome_id" to
           list of type "feature_id", parameter "feature_id_to_function" of
           mapping from type "feature_id" to String, parameter
           "genome_ref_to_sci_name" of mapping from type "data_obj_ref" to
           String, parameter "genome_ref_to_obj_name" of mapping from type
           "data_obj_ref" to String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN SpeciesTreeToFASTA

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11) 

        # init and clean up params
        tree_ref = params['tree_ref']
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
        merge_fasta_files = params['merge_fasta_files']

        # Defaults
        if console == None:
            console = []
        if invalid_msgs == None:
            invalid_msgs = []
        if residue_type == None:
            residue_type = 'NUC'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = 'g:%%genome_ref%%.f:%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0
        if merge_fasta_files == None or merge_fasta_files[0:1].upper() == 'T':
            merge_fasta_files = True
        else:
            merge_fasta_files = False

        # init and simplify
        fasta_file_path_list = []
        feature_ids_by_genome_id = dict()
        feature_id_to_function = dict()
        genome_ref_to_sci_name = dict()
        genome_ref_to_obj_name = dict()
        residue_type = residue_type[0:1].upper()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch

        # get SpeciesTree object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            tree_object = ws.get_objects2({'objects':[{'ref':tree_ref}]})['data'][0]['data']
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        if 'type' not in tree_object or tree_object['type'] != 'SpeciesTree':
            raise ValueError ("must have SpeciesTree type for tree object "+tree_ref)
        if 'ws_refs' not in tree_object:
            raise ValueError ("must have ws_refs for tree object "+tree_ref)
            
            
        # iterate through tree genome members
        genome_ids = tree_object['ws_refs'].keys()
        for genome_i in range(len(genome_ids)):
            feature_sequence_found = False
            genome_id = genome_ids[genome_i]
            feature_ids_by_genome_id[genome_id] = []

            if not tree_object['ws_refs'][genome_id].get('g'):
                raise ValueError('SpeciesTreeToFASTA() malformed ws_refs')
                #to get the full stack trace: traceback.format_exc()       
            else:
                genome_ref = tree_object['ws_refs'][genome_id]['g'][0]


            # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
            #
            #records = []
            this_file = file
            if len(genome_ids) > 1 and not merge_fasta_files:
                this_file = this_file+'.'+genome_id
                this_file.replace('/', '-')
            fasta_file_path = os.path.join(dir, this_file)
            #self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG
            
            if genome_i == 0 or not merge_fasta_files:
                fasta_file_handle = open(fasta_file_path, 'w', 0)
                self.log(console, 'KB SDK data2file GenomeSet2FASTA(): writing fasta file: '+fasta_file_path)
            # DEBUG
            self.log(console, "ADDING GENOME: "+str(genome_i+1)+" of "+str(len(genome_ids))+" "+genome_ids[genome_i])


            # get genome object
            try:
                ws = workspaceService(self.workspaceURL, token=ctx['token'])
                #genome_object = ws.get_objects([{'ref':genome_ref}])[0]['data']
                genome_obj_base = ws.get_objects2({'objects':[{'ref':genome_ref}]})['data'][0]
                genome_object = genome_obj_base['data']
                genome_info = genome_obj_base['info']
                genome_ref_to_obj_name[genome_ref] = genome_info[NAME_I]
            except Exception as e:
                raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                #to get the full stack trace: traceback.format_exc()

            # set sci name
            genome_ref_to_sci_name[genome_ref] = genome_object['scientific_name']
            feature_id_to_function[genome_ref] = dict()

            # set features list ('features' only contains CDS features')
            target_features = []
            if feature_type == 'CDS' or residue_type == 'P' or 'non_coding_features' not in genome_object:
                target_features = genome_object['features']
            else:
                target_features = genome_object['features'] + genome_object['non_coding_features']


            # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
            #
            #records = []
                        
            for feature in target_features:
                fid = feature['id']

                # set function
                if 'function' in feature:  # Feature <= 2.3 (Genome <= 8.2)
                    feature_id_to_function[genome_ref][fid] = feature['function']
                elif 'functions' in feature:  # Feature >= 3.0 (Genome >= 9.0)
                    feature_id_to_function[genome_ref][fid] = "; ".join(feature['functions'])
                else:
                    feature_id_to_function[genome_ref][fid] = 'N/A'

                #self.log (console, "FID:'"+str(fid)+"' FXN: '"+str(feature_id_to_function[genome_ref][fid]))  # DEBUG

                #if feature_type == 'ALL' or feature_type == feature['type']:
                if True:  # don't want to deal with changing indentation

                    # protein recs
                    if residue_type == 'P':
                        #if feature['type'] != 'CDS':
                        #    continue
                        if 'protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == '':
                            #self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
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
                            feature_ids_by_genome_id[genome_id].append(feature['id'])
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if feature_type == 'CDS' and ('protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == ''):
                            continue
                        elif 'dna_sequence' not in feature or feature['dna_sequence'] == None or feature['dna_sequence'] == '':
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_id, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_id, genome_ref)
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
                            feature_ids_by_genome_id[genome_id].append(feature['id'])
                            fasta_file_handle.write(rec)


            if genome_i == (len(genome_ids)-1) or not merge_fasta_files:
                self.log(console,"CLOSING FILE: '"+fasta_file_path+"'")  # DEBUG
                fasta_file_handle.close()
                fasta_file_path_list.append(fasta_file_path)


            # report if no features found
            if not feature_sequence_found:
                self.log(invalid_msgs, "No sequence records found in Genome "+genome_id+" of residue_type: "+residue_type+", feature_type: "+feature_type)


        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path_list'] = fasta_file_path_list
        returnVal['feature_ids_by_genome_id'] = feature_ids_by_genome_id
        returnVal['feature_id_to_function'] = feature_id_to_function
        returnVal['genome_ref_to_sci_name'] = genome_ref_to_sci_name
        returnVal['genome_ref_to_obj_name'] = genome_ref_to_obj_name
        #END SpeciesTreeToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method SpeciesTreeToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def FeatureSetToFASTA(self, ctx, params):
        """
        :param params: instance of type "FeatureSetToFASTA_Params"
           (FeatureSetToFASTA() Params) -> structure: parameter
           "featureSet_ref" of type "data_obj_ref", parameter "file" of type
           "path_type", parameter "dir" of type "path_type", parameter
           "console" of list of type "log_msg", parameter "invalid_msgs" of
           list of type "log_msg", parameter "residue_type" of String,
           parameter "feature_type" of String, parameter "record_id_pattern"
           of type "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type "FeatureSetToFASTA_Output"
           (FeatureSetToFASTA() Output) -> structure: parameter
           "fasta_file_path" of type "path_type", parameter
           "feature_ids_by_genome_ref" of mapping from type "data_obj_ref" to
           list of type "feature_id", parameter "feature_id_to_function" of
           mapping from type "feature_id" to String, parameter
           "genome_ref_to_sci_name" of mapping from type "data_obj_ref" to
           String, parameter "genome_ref_to_obj_name" of mapping from type
           "data_obj_ref" to String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN FeatureSetToFASTA

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11) 

        # init and clean up params
        featureSet_ref = params['featureSet_ref']
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
            residue_type = 'NUC'
        if feature_type == None:
            feature_type = 'ALL';
        if record_id_pattern == None:
            record_id_pattern = 'g:%%genome_ref%%.f:%%feature_id%%'
        if record_desc_pattern == None:
            record_desc_pattern = '[%%genome_ref%%]'
        if case == None:
            case = 'UPPER'
        if linewrap == None:
            linewrap = 0

        # init and simplify
        feature_ids_by_genome_ref = dict()
        feature_id_to_function = dict()
        genome_ref_to_sci_name = dict()
        genome_ref_to_obj_name = dict()
        residue_type = residue_type[0:1].upper()
        feature_type = feature_type.upper()
        case = case[0:1].upper()
        feature_sequence_found = False
        
        def record_header_sub(str, feature_id, genome_id, genome_ref):
            str = str.replace('%%feature_id%%', feature_id)
            str = str.replace('%%genome_id%%', genome_id)
            str = str.replace('%%genome_ref%%', genome_ref)
            return str

        if file == None:
            file = 'runfile.fasta'
        if dir == None:
            dir = self.scratch
        fasta_file_path = os.path.join(dir, file)
        self.log(console, 'KB SDK data2file FeatureSetToFASTA(): writing fasta file: '+fasta_file_path)

        # get featureSet object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #featureSet_object = ws.get_objects([{'ref':featureSet_ref}])[0]['data']
            featureSet_object = ws.get_objects2({'objects':[{'ref':featureSet_ref}]})['data'][0]['data']
        except Exception as e:
            raise ValueError('Unable to fetch featureSet object from workspace: ' + str(e))
            #to get the full stack trace: traceback.format_exc()

        featureSet_features = featureSet_object['elements']
        genome2features = {}
        featureSetLookup = {}
        for fId in featureSet_features.keys():
            for genome_ref in featureSet_features[fId]:
                #genome_ref = featureSet_features[fId][0]  # get more than first one
                if genome_ref not in genome2features.keys():
                    genome2features[genome_ref] = []
                    featureSetLookup[genome_ref] = dict()
                genome2features[genome_ref].append(fId)
                featureSetLookup[genome_ref][fId] = True

        # write file
        with open (fasta_file_path, 'w', 0)  as fasta_file_handle:

            for genome_ref in genome2features.keys():
                feature_ids_by_genome_ref[genome_ref] = []

                # get genome object
                try:
                    ws = workspaceService(self.workspaceURL, token=ctx['token'])
                    #genome_object = ws.get_objects([{'ref':genome_ref}])[0]['data']
                    genome_obj_base = ws.get_objects2({'objects':[{'ref':genome_ref}]})['data'][0]
                    genome_object = genome_obj_base['data']
                    genome_info = genome_obj_base['info']
                    genome_ref_to_obj_name[genome_ref] = genome_info[NAME_I]

                    obj_type = genome_info[TYPE_I].split('-')[0]
                except Exception as e:
                    raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
                    #to get the full stack trace: traceback.format_exc()
            
                # different behavior for genome and ama
                these_features = []
                if obj_type == 'KBaseGenomes.Genome':
                    genome_ref_to_sci_name[genome_ref] = genome_object['scientific_name']
                    # set features list ('features' only contains CDS features')
                    if feature_type == 'CDS' or residue_type == 'P' or 'non_coding_features' not in genome_object:
                        these_features = genome_object['features']
                    else:
                        these_features = genome_object['features'] + genome_object['non_coding_features']
                else:
                    genome_ref_to_sci_name[genome_ref] = genome_info[NAME_I]
                    these_features = self._get_features_from_AnnotatedMetagenomeAssembly (ctx, genome_ref)
                    
                feature_id_to_function[genome_ref] = dict()

                    
                # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
                #
                #records = []
                for feature in these_features:
                    fid = feature['id']
                
                    try:
                        in_set = featureSetLookup[genome_ref][fid]
                    except:
                        continue

                    # set function
                    if 'function' in feature:  # Feature <= 2.3 (Genome <= 8.2)
                        feature_id_to_function[genome_ref][fid] = feature['function']
                    elif 'functions' in feature:  # Feature >= 3.0 (Genome >= 9.0)
                        feature_id_to_function[genome_ref][fid] = "; ".join(feature['functions'])
                    else:
                        feature_id_to_function[genome_ref][fid] = 'N/A'

                    # protein recs
                    if residue_type == 'P':
                        #if feature['type'] != 'CDS':
                        #    continue
                        if 'protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == '':
                            #self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
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
                            feature_ids_by_genome_ref[genome_ref].append(feature['id'])
                            fasta_file_handle.write(rec)

                    # nuc recs
                    else:
                        if feature_type == 'CDS' and ('protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == ''):
                            continue
                        elif 'dna_sequence' not in feature or feature['dna_sequence'] == None or feature['dna_sequence'] == '':
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        elif 'parent_gene' in feature:
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, fid, genome_ref, genome_ref)
                            rec_desc = record_header_sub(record_desc_pattern, fid, genome_ref, genome_ref)
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
                            feature_ids_by_genome_ref[genome_ref].append(feature['id'])
                            fasta_file_handle.write(rec)



# HERE
#            if genome_i == (len(genome_ids)-1) or not merge_fasta_files:
#                self.log(console,"CLOSING FILE: '"+fasta_file_path+"'")  # DEBUG
#                fasta_file_handle.close()
#                fasta_file_path_list.append(fasta_file_path)


        # report if no features found
        if not feature_sequence_found:
            self.log(invalid_msgs, "No sequence records found in Genome "+genome_ref+" of residue_type: "+residue_type+", feature_type: "+feature_type)

        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path'] = fasta_file_path
        returnVal['feature_ids_by_genome_ref'] = feature_ids_by_genome_ref
        returnVal['feature_id_to_function'] = feature_id_to_function
        returnVal['genome_ref_to_sci_name'] = genome_ref_to_sci_name
        returnVal['genome_ref_to_obj_name'] = genome_ref_to_obj_name
        #END FeatureSetToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method FeatureSetToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def AnnotatedMetagenomeAssemblyToFASTA(self, ctx, params):
        """
        :param params: instance of type
           "AnnotatedMetagenomeAssemblyToFASTA_Params"
           (AnnotatedMetagenomeAssemblyToFASTA() Params) -> structure:
           parameter "ama_ref" of type "data_obj_ref", parameter "file" of
           type "path_type", parameter "dir" of type "path_type", parameter
           "console" of list of type "log_msg", parameter "invalid_msgs" of
           list of type "log_msg", parameter "residue_type" of String,
           parameter "feature_type" of String, parameter "record_id_pattern"
           of type "pattern_type", parameter "record_desc_pattern" of type
           "pattern_type", parameter "case" of String, parameter "linewrap"
           of Long
        :returns: instance of type
           "AnnotatedMetagenomeAssemblyToFASTA_Output"
           (AnnotatedMetagenomeAssemblyToFASTA() Output) -> structure:
           parameter "fasta_file_path" of type "path_type", parameter
           "feature_ids" of list of type "feature_id", parameter
           "feature_id_to_function" of mapping from type "feature_id" to
           String, parameter "ama_ref_to_obj_name" of mapping from type
           "data_obj_ref" to String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN AnnotatedMetagenomeAssemblyToFASTA

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11) 

        # init and clean up params
        ama_ref = params['ama_ref']
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
        feature_id_to_function = dict()
        ama_ref_to_obj_name = dict()
        feature_sequence_found = False
        residue_type = residue_type[0:1].upper()
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
        self.log(console, 'KB SDK data2file AnnotatedMetagenomeAssemby2FASTA(): writing fasta file: '+fasta_file_path)

        # get genome object
        try:
            ws = workspaceService(self.workspaceURL, token=ctx['token'])
            #ama_object = ws.get_objects2({'objects':[{'ref':ama_ref}]})['data'][0]['data']
            #ama_object_data = ama_object['data']
            #ama_object_info = ama_object['info']
            
            ama_object_info = ws.get_object_info3 ({'objects':[{'ref':ama_ref}]})['infos'][0]
            ama_object_name = ama_object_info[NAME_I]
        except Exception as e:
            raise ValueError('Unable to fetch input_one_name object from workspace: ' + str(e))
        #to get the full stack trace: traceback.format_exc()

        # set obj name
        ama_ref_to_obj_name[ama_ref] = ama_object_name
        feature_id_to_function[ama_ref] = dict()


        # Get feature recs from AMA (this is very different than for Genome)
        #
        ama_features = self._get_features_from_AnnotatedMetagenomeAssembly (ctx, ama_ref)


        # FIX: should I write recs as we go to reduce memory footprint, or is a single buffer write much faster?  Check later.
        #
        #records = []
        self.log(console,"FASTA_FILE_PATH'"+fasta_file_path+"'\n")  # DEBUG

        with open(fasta_file_path, 'w', 0) as fasta_file_handle:
                        
            for feature in ama_features:
                # set function
                if 'function' in feature:  # Feature <= 2.3
                    feature_id_to_function[ama_ref][feature['id']] = feature['function']
                elif 'functions' in feature:  # Feature >= 3.0
                    feature_id_to_function[ama_ref][feature['id']] = "; ".join(feature['functions'])
                else:
                    feature_id_to_function[ama_ref][feature['id']] = 'N/A'
                
                #if feature_type == 'ALL' or feature_type == feature['type']:
                if True:  # don't want to deal with changing indentation

                    # protein recs
                    if residue_type == 'P':
                        if feature['type'] != 'CDS':
                            continue
                        #if 'protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == '':
                        #    #self.log(invalid_msgs, "bad CDS feature "+feature['id']+": No protein_translation field.")
                        #    continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, feature['id'], ama_object_name)
                            rec_desc = record_header_sub(record_desc_pattern, feature['id'], ama_object_name)

                            if feature.get('protein_translation'):
                                seq = feature['protein_translation']
                            elif feature.get('dna_sequence'):
                                seq = self.TranslateNucToProtSeq(ctx,
                                                                 {'nuc_seq': feature['dna_sequence'],
                                                                  'genetic_code': '11'}
                                                                 )[0]['prot_seq']
                            else:
                                raise ValueError ("no sequence information for gene "+feature['id'])
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
                        if feature_type == 'CDS' and ('protein_translation' not in feature or feature['protein_translation'] == None or feature['protein_translation'] == ''):
                            continue
                        elif 'dna_sequence' not in feature or feature['dna_sequence'] == None or feature['dna_sequence'] == '':
                            self.log(invalid_msgs, "bad feature "+feature['id']+": No dna_sequence field.")
                        elif 'parent_gene' in feature:
                            continue
                        else:
                            feature_sequence_found = True
                            rec_id = record_header_sub(record_id_pattern, feature['id'], ama_object_name)
                            rec_desc = record_header_sub(record_desc_pattern, feature['id'], ama_object_name)
                            if feature.get('dna_sequence'):
                                seq = feature['dna_sequence']
                            else:
                                raise ValueError ("no sequence information for gene "+feature['id'])
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
            self.log(invalid_msgs, "No sequence records found in AnnotatedMetagenomeAssembly "+ama_object_name+" of residue_type: "+residue_type+", feature_type: "+feature_type)
        #else:
        #    SeqIO.write(records, fasta_file_path, "fasta")


        # build returnVal
        #
        returnVal = dict()
        returnVal['fasta_file_path'] = fasta_file_path
        returnVal['feature_ids'] = feature_ids
        returnVal['feature_id_to_function'] = feature_id_to_function
        returnVal['ama_ref_to_obj_name'] = ama_ref_to_obj_name
        #END AnnotatedMetagenomeAssemblyToFASTA

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method AnnotatedMetagenomeAssemblyToFASTA return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK", 'message': "", 'version': self.VERSION, 
                     'git_url': self.GIT_URL, 'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
