# -*- coding: utf-8 -*-
import unittest
import os
import json
import time
import requests
import uuid
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint

from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from KBaseDataObjectToFileUtils.KBaseDataObjectToFileUtilsImpl import KBaseDataObjectToFileUtils
from KBaseDataObjectToFileUtils.KBaseDataObjectToFileUtilsServer import MethodContext
from installed_clients.authclient import KBaseAuth as _KBaseAuth


class KBaseDataObjectToFileUtilsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('KBaseDataObjectToFileUtils'):
            cls.cfg[nameval[0]] = nameval[1]
        authServiceUrl = cls.cfg.get('auth-service-url',
                "https://kbase.us/services/authorization/Sessions/Login")
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'KBaseDataObjectToFileUtils',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.scratch = cls.cfg['scratch']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = KBaseDataObjectToFileUtils(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_KBaseDataObjectToFileUtils_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # call this method to get the WS object info of a Genome
    #   (will upload the example data if this is the first time the method is called during tests)
    def getGenomeInfo(self, genome_basename, item_i=0):
        if hasattr(self.__class__, 'genomeInfo_list'):
            try:
                info = self.__class__.genomeInfo_list[item_i]
                name = self.__class__.genomeName_list[item_i]
                if info != None:
                    if name != genome_basename:
                        self.__class__.genomeInfo_list[item_i] = None
                        self.__class__.genomeName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform genbank to kbase genome object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        genome_data_file = 'data/genomes/'+genome_basename+'.gbff.gz'
        genome_file = os.path.join(shared_dir, os.path.basename(genome_data_file))
        shutil.copy(genome_data_file, genome_file)

        SERVICE_VER = 'release'
        #SERVICE_VER = 'dev'
        GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                             token=self.getContext()['token'],
                             service_ver=SERVICE_VER
                         )
        print ("UPLOADING genome: "+genome_basename+" to WORKSPACE "+self.getWsName()+" ...")
        genome_upload_result = GFU.genbank_to_genome({'file': {'path': genome_file },
                                                      'workspace_name': self.getWsName(),
                                                      'genome_name': genome_basename
                                                  })
#                                                  })[0]
        pprint(genome_upload_result)
        genome_ref = genome_upload_result['genome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': genome_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'genomeInfo_list'):
            self.__class__.genomeInfo_list = []
            self.__class__.genomeName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.genomeInfo_list[i]
            except:
                self.__class__.genomeInfo_list.append(None)
                self.__class__.genomeName_list.append(None)

        self.__class__.genomeInfo_list[item_i] = new_obj_info
        self.__class__.genomeName_list[item_i] = genome_basename
        return new_obj_info

    # call this method to get the WS object info of an AnnotatedMetagenomeAssembly
    #   (will upload the example data if this is the first time the method is called during tests)
    def getAMAInfo(self, ama_basename, item_i=0):
        if hasattr(self.__class__, 'amaInfo_list'):
            try:
                info = self.__class__.amaInfo_list[item_i]
                name = self.__class__.amaName_list[item_i]
                if info != None:
                    if name != ama_basename:
                        self.__class__.amaInfo_list[item_i] = None
                        self.__class__.amaName_list[item_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform GFF+FNA to kbase AMA object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        ama_gff_srcfile = 'data/amas/'+ama_basename+'.gff'
        ama_fna_srcfile = 'data/amas/'+ama_basename+'.fa'
        ama_gff_dstfile = os.path.join(shared_dir, os.path.basename(ama_gff_srcfile))
        ama_fna_dstfile = os.path.join(shared_dir, os.path.basename(ama_fna_srcfile))
        shutil.copy(ama_gff_srcfile, ama_gff_dstfile)
        shutil.copy(ama_fna_srcfile, ama_fna_dstfile)

        try:
            SERVICE_VER = 'release'
            #SERVICE_VER = 'dev'
            GFU = GenomeFileUtil(os.environ['SDK_CALLBACK_URL'],
                                 token=self.getContext()['token'],
                                 service_ver=SERVICE_VER
            )
        except:
            raise ValueError ("unable to obtain GenomeFileUtil client")
        print ("UPLOADING AMA: "+ama_basename+" to WORKSPACE "+self.getWsName()+" ...")
        ama_upload_params = {
            "workspace_name": self.getWsName(),
            "genome_name": ama_basename,
            "fasta_file": {"path": ama_fna_dstfile},
            "gff_file": {"path": ama_gff_dstfile},
            "source": "GFF",
            "scientific_name": "TEST AMA",
            "generate_missing_genes": "True"
        }        
        try:
            ama_upload_result = GFU.fasta_gff_to_metagenome(ama_upload_params)
        except:
            raise ValueError("unable to upload test AMA data object")
        print ("AMA UPLOADED")
        pprint(ama_upload_result)

        ama_ref = ama_upload_result['metagenome_ref']
        new_obj_info = self.getWsClient().get_object_info_new({'objects': [{'ref': ama_ref}]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'amaInfo_list'):
            self.__class__.amaInfo_list = []
            self.__class__.amaName_list = []
        for i in range(item_i+1):
            try:
                assigned = self.__class__.amaInfo_list[i]
            except:
                self.__class__.amaInfo_list.append(None)
                self.__class__.amaName_list.append(None)

        self.__class__.amaInfo_list[item_i] = new_obj_info
        self.__class__.amaName_list[item_i] = ama_basename
        return new_obj_info


    # call this method to get the WS object info of a Tree
    #   (will upload the example data if this is the first time the method is called during tests)
    def getTreeInfo(self, tree_basename, lib_i=0, genome_ref_map=None):
        if hasattr(self.__class__, 'treeInfo_list'):
            try:
                info = self.__class__.treeInfo_list[lib_i]
                name = self.__class__.treeName_list[lib_i]
                if info != None:
                    if name != tree_basename:
                        self.__class__.treeInfo_list[lib_i] = None
                        self.__class__.treeName_list[lib_i] = None
                    else:
                        return info
            except:
                pass

        # 1) transform json to kbase Tree object and upload to ws
        shared_dir = "/kb/module/work/tmp"
        tree_data_file = 'data/trees/'+tree_basename+'.json'
        tree_file = os.path.join(shared_dir, os.path.basename(tree_data_file))
        shutil.copy(tree_data_file, tree_file)

        # create object
        with open (tree_file, 'r', 0) as tree_fh:
            tree_obj = json.load(tree_fh)

        # update genome_refs
        if genome_ref_map != None:
            for label_id in tree_obj['default_node_labels']:
                for old_genome_ref in genome_ref_map.keys():
                    tree_obj['default_node_labels'][label_id] = tree_obj['default_node_labels'][label_id].replace(old_genome_ref, genome_ref_map[old_genome_ref])
            for label_id in tree_obj['ws_refs'].keys():
                new_genome_refs = []
                for old_genome_ref in tree_obj['ws_refs'][label_id]['g']:
                    new_genome_refs.append(genome_ref_map[old_genome_ref])
                tree_obj['ws_refs'][label_id]['g'] = new_genome_refs

        provenance = [{}]
        new_obj_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseTrees.Tree',
                    'data': tree_obj,
                    'name': tree_basename+'.test_TREE',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        # 2) store it
        if not hasattr(self.__class__, 'treeInfo_list'):
            self.__class__.treeInfo_list = []
            self.__class__.treeName_list = []
        for i in range(lib_i+1):
            try:
                assigned = self.__class__.treeInfo_list[i]
            except:
                self.__class__.treeInfo_list.append(None)
                self.__class__.treeName_list.append(None)

        self.__class__.treeInfo_list[lib_i] = new_obj_info
        self.__class__.treeName_list[lib_i] = tree_basename
        return new_obj_info


    #
    # NOTE: According to Python unittest naming rules test method names should start from 'test'.
    #

    #### TranslateNucToProtSeq_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_TranslateNucToProtSeq_01()')
    def test_KBaseDataObjectToFileUtils_TranslateNucToProtSeq_01(self):
        # Prepare test objects in workspace if needed using 
        # self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        
        # E. coli K-12 MG1655 (start modified from GTG to ATG)
        nuc_seq = 'ATGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAACACGCTGGCCCTGTACGCGCCAAACCGTTTTGTCCTCGATTGGGTACGGGACAAGTACCTTAATAATATCAATGGACTGCTAACCAGTTTCTGCGGAGCGGATGCCCCACAGCTGCGTTTTGAAGTCGGCACCAAACCGGTGACGCAAACGCCACAAGCGGCAGTGACGAGCAACGTCGCGGCCCCTGCACAGGTGGCGCAAACGCAGCCGCAACGTGCTGCGCCTTCTACGCGCTCAGGTTGGGATAACGTCCCGGCCCCGGCAGAACCGACCTATCGTTCTAACGTAAACGTCAAACACACGTTTGATAACTTCGTTGAAGGTAAATCTAACCAACTGGCGCGCGCGGCGGCTCGCCAGGTGGCGGATAACCCTGGCGGTGCCTATAACCCGTTGTTCCTTTATGGCGGCACGGGTCTGGGTAAAACTCACCTGCTGCATGCGGTGGGTAACGGCATTATGGCGCGCAAGCCGAATGCCAAAGTGGTTTATATGCACTCCGAGCGCTTTGTTCAGGACATGGTTAAAGCCCTGCAAAACAACGCGATCGAAGAGTTTAAACGCTACTACCGTTCCGTAGATGCACTGCTGATCGACGATATTCAGTTTTTTGCTAATAAAGAACGATCTCAGGAAGAGTTTTTCCACACCTTCAACGCCCTGCTGGAAGGTAATCAACAGATCATTCTCACCTCGGATCGCTATCCGAAAGAGATCAACGGCGTTGAGGATCGTTTGAAATCCCGCTTCGGTTGGGGACTGACTGTGGCGATCGAACCGCCAGAGCTGGAAACCCGTGTGGCGATCCTGATGAAAAAGGCCGACGAAAACGACATTCGTTTGCCGGGCGAAGTGGCGTTCTTTATCGCCAAGCGTCTACGATCTAACGTACGTGAGCTGGAAGGGGCGCTGAACCGCGTCATTGCCAATGCCAACTTTACCGGACGGGCGATCACCATCGACTTCGTGCGTGAGGCGCTGCGCGACTTGCTGGCATTGCAGGAAAAACTGGTCACCATCGACAATATTCAGAAGACGGTGGCGGAGTACTACAAGATCAAAGTCGCGGATCTCCTTTCCAAGCGTCGATCCCGCTCGGTGGCGCGTCCGCGCCAGATGGCGATGGCGCTGGCGAAAGAGCTGACTAACCACAGTCTGCCGGAGATTGGCGATGCGTTTGGTGGCCGTGACCACACGACGGTGCTTCATGCCTGCCGTAAGATCGAGCAGTTGCGTGAAGAGAGCCACGATATCAAAGAAGATTTTTCAAATTTAATCAGAACATTGTCATCGTAA'
        prot_seq = 'MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLARAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVYMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEFFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELETRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFTGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRSRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREESHDIKEDFSNLIRTLSS'

        parameters = { 'nuc_seq': nuc_seq,
                       'genetic_code': '11'
                     }
        ret = self.getImpl().TranslateNucToProtSeq(self.getContext(), parameters)[0]
        self.assertEqual(ret['prot_seq'], prot_seq)
        pass


    #### ParseFastaStr_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_ParseFastaStr_01()')
    def test_KBaseDataObjectToFileUtils_ParseFastaStr_01(self):
        
        # E. coli K-12 MG1655 (start modified from GTG to ATG)
        seq_id = 'test1'
        seq_desc = '[E. coli K-12 MG1655] DnaA'
        prot_seq_nl = "MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDW\nVRDKYLNNINGLLTSFCGADAPQLRFEVGTKPVTQTPQAAVTSNVAAPAQ\nVAQTQPQRAAPSTRSGWDNVPAPAEPTYRSNVNVKHTFDNFVEGKSNQLA\nRAAARQVADNPGGAYNPLFLYGGTGLGKTHLLHAVGNGIMARKPNAKVVY\nMHSERFVQDMVKALQNNAIEEFKRYYRSVDALLIDDIQFFANKERSQEEF\nFHTFNALLEGNQQIILTSDRYPKEINGVEDRLKSRFGWGLTVAIEPPELE\nTRVAILMKKADENDIRLPGEVAFFIAKRLRSNVRELEGALNRVIANANFT\nGRAITIDFVREALRDLLALQEKLVTIDNIQKTVAEYYKIKVADLLSKRRS\nRSVARPRQMAMALAKELTNHSLPEIGDAFGGRDHTTVLHACRKIEQLREE\nSHDIKEDFSNLIRTLSS\n"
        prot_seq = prot_seq_nl.replace("\n",'')
        fasta_str = '>'+seq_id+' '+seq_desc+"\n"+prot_seq_nl

        parameters = { 'fasta_str': fasta_str,
                       'residue_type': 'prot',
                       'case': 'U',
                       'console': [],
                       'invalid_msgs': []
                     }
        ret = self.getImpl().ParseFastaStr(self.getContext(), parameters)[0]
        self.assertEqual(ret['id'], seq_id)
        self.assertEqual(ret['desc'], seq_desc)
        self.assertEqual(ret['seq'], prot_seq)
        pass
        

    #### GenomeToFASTA_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_GenomeToFASTA_01()')
    def test_KBaseDataObjectToFileUtils_GenomeToFASTA_01(self):
        test_name = 'GenomeToFASTA_01'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        """
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])        

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'genome_ref':          genome_ref_0,
                'file':                test_name+'-test_genome.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }
        ret = self.getImpl().GenomeToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'])
        self.assertIsNotNone(ret['feature_ids'])
        self.assertNotEqual(len(ret['feature_ids']), 0)
        pass
        

    #### GenomeToFASTA_02()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_GenomeToFASTA_02()')
    def test_KBaseDataObjectToFileUtils_GenomeToFASTA_02(self):
        test_name = 'GenomeToFASTA_02'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        """
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])        

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'genome_ref':          genome_ref_0,
                'file':                test_name+'-test_genome.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'nuc',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }
        ret = self.getImpl().GenomeToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'])
        self.assertIsNotNone(ret['feature_ids'])
        self.assertNotEqual(len(ret['feature_ids']), 0)
        pass
        

    #### GenomeSetToFASTA_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_GenomeSetToFASTA_01()')
    def test_KBaseDataObjectToFileUtils_GenomeSetToFASTA_01(self):
        test_name = 'GenomeSetToFASTA_01'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genome_ref_2 = reference_prok_genomes_WS+'/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'
        genome_ref_3 = reference_prok_genomes_WS+'/GCF_900129775.1/1'  # Halobaculum gomorrense (16 contigs)
        """
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])        
        
        genome_id_feature_id_delim = '.f:'

        genomeSet_obj = { 'description': 'test genomeSet',
                          'elements': { 'genome_0': { 'ref': genome_ref_0 },
                                        'genome_1': { 'ref': genome_ref_1 },
                                        'genome_2': { 'ref': genome_ref_2 }
                                    }
                        }
        provenance = [{}]
        genomeSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseSearch.GenomeSet',
                    'data': genomeSet_obj,
                    'name': 'test_genomeSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        genomeSet_ref = str(genomeSet_info[WSID_I])+'/'+str(genomeSet_info[OBJID_I])+'/'+str(genomeSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'genomeSet_ref':       genomeSet_ref,
                'file':                test_name+'-test_genomeSet.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }
        ret = self.getImpl().GenomeSetToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path_list'][0])
        self.assertIsNotNone(ret['feature_ids_by_genome_id'])
        self.assertNotEqual(len(ret['feature_ids_by_genome_id'].keys()), 0)
        pass
        

    #### GenomeSetToFASTA_02()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_GenomeSetToFASTA_02()')
    def test_KBaseDataObjectToFileUtils_GenomeSetToFASTA_02(self):
        test_name = 'GenomeSetToFASTA_02'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genome_ref_2 = reference_prok_genomes_WS+'/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'
        genome_ref_3 = reference_prok_genomes_WS+'/GCF_900129775.1/1'  # Halobaculum gomorrense (16 contigs)
        """
        
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])        
        
        genome_id_feature_id_delim = '.f:'

        genomeSet_obj = { 'description': 'test genomeSet',
                          'elements': { 'genome_0': { 'ref': genome_ref_0 },
                                        'genome_1': { 'ref': genome_ref_1 },
                                        'genome_2': { 'ref': genome_ref_2 }
                                    }
                        }
        provenance = [{}]
        genomeSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseSearch.GenomeSet',
                    'data': genomeSet_obj,
                    'name': 'test_genomeSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        genomeSet_ref = str(genomeSet_info[WSID_I])+'/'+str(genomeSet_info[OBJID_I])+'/'+str(genomeSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'genomeSet_ref':       genomeSet_ref,
                'file':                test_name+'-test_genomeSet.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'nuc',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }
        ret = self.getImpl().GenomeSetToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path_list'][0])
        self.assertIsNotNone(ret['feature_ids_by_genome_id'])
        self.assertNotEqual(len(ret['feature_ids_by_genome_id'].keys()), 0)
        pass
        

    #### SpeciesTreeToFASTA_01()
    ##
    # HIDE @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_SpeciesTreeToFASTA_01()')
    def test_KBaseDataObjectToFileUtils_SpeciesTreeToFASTA_01(self):
        test_name = 'SpeciesTreeToFASTA_01'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple


        # input_data
        genomeInfo_3 = self.getGenomeInfo('GCF_000287295.1_ASM28729v1_genomic', 3)  # Candidatus Carsonella ruddii HT isolate Thao2000
        genomeInfo_4 = self.getGenomeInfo('GCF_000306885.1_ASM30688v1_genomic', 4)  # Wolbachia endosymbiont of Onchocerca ochengi
        genomeInfo_5 = self.getGenomeInfo('GCF_001439985.1_wTPRE_1.0_genomic',  5)  # Wolbachia endosymbiont of Trichogramma pretiosum
        genomeInfo_6 = self.getGenomeInfo('GCF_000022285.1_ASM2228v1_genomic',  6)  # Wolbachia sp. wRi

        genome_ref_3 = '/'.join([str(genomeInfo_3[WSID_I]),
                                 str(genomeInfo_3[OBJID_I]),
                                 str(genomeInfo_3[VERSION_I])])
        genome_ref_4 = '/'.join([str(genomeInfo_4[WSID_I]),
                                 str(genomeInfo_4[OBJID_I]),
                                 str(genomeInfo_4[VERSION_I])])
        genome_ref_5 = '/'.join([str(genomeInfo_5[WSID_I]),
                                 str(genomeInfo_5[OBJID_I]),
                                 str(genomeInfo_5[VERSION_I])])
        genome_ref_6 = '/'.join([str(genomeInfo_6[WSID_I]),
                                 str(genomeInfo_6[OBJID_I]),
                                 str(genomeInfo_6[VERSION_I])])

        # upload Tree
        genome_refs_map = { '23880/3/1': genome_ref_3,
                            '23880/4/1': genome_ref_4,
                            '23880/5/1': genome_ref_5,
                            '23880/6/1': genome_ref_6
                        }
        obj_info = self.getTreeInfo('Tiny_things.SpeciesTree', 0, genome_refs_map)
        tree_ref = str(obj_info[WSID_I])+'/'+str(obj_info[OBJID_I])+'/'+str(obj_info[VERSION_I])

        #feature_id_0 = 'A355_RS00030'   # F0F1 ATP Synthase subunit B
        #feature_id_1 = 'WOO_RS00195'    # F0 ATP Synthase subunit B
        #feature_id_2 = 'AOR14_RS04755'  # F0 ATP Synthase subunit B
        #feature_id_3 = 'WRI_RS01560'    # F0 ATP Synthase subunit B


        genome_id_feature_id_delim = '.f:'

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'tree_ref':            tree_ref,
                'file':                test_name+'-test_genomeSet.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+genome_id_feature_id_delim+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }
        ret = self.getImpl().SpeciesTreeToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path_list'][0])
        self.assertIsNotNone(ret['feature_ids_by_genome_id'])
        self.assertNotEqual(len(ret['feature_ids_by_genome_id'].keys()), 0)
        pass
        

    #### FeatureSetToFASTA_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_FeatureSetToFASTA_01()')
    def test_KBaseDataObjectToFileUtils_FeatureSetToFASTA_01(self):
        test_name = 'FeatureSetToFASTA_01'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genome_ref_2 = reference_prok_genomes_WS+'/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'
        genome_ref_3 = reference_prok_genomes_WS+'/GCF_900129775.1/1'  # Halobaculum gomorrense (16 contigs)
        """
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])        
        

        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_0 = '/'.join([str(amaInfo_0[WSID_I]),
                              str(amaInfo_0[OBJID_I]),
                              str(amaInfo_0[VERSION_I])])

        genome_id_feature_id_delim = '.f:'
        feature_id_1 = 'AWN69_RS07145' # DnaA
        feature_id_2 = 'DVMF_RS00005'  # DnaA
        feature_id_3 = 'A6701_RS00005' # DnaA
        ama_feature_id_1 = '4_27'
        
        featureSet_obj = { 'description': 'test genomeSet',
                           'element_ordering': [
                               feature_id_1,
                               feature_id_2,
                               feature_id_3,
                               ama_feature_id_1
                           ],
                           'elements': { 
                               feature_id_1: [genome_ref_0],
                               feature_id_2: [genome_ref_1],
                               feature_id_3: [genome_ref_2],
                               ama_feature_id_1: [ama_ref_0]
                           }
                        }
        provenance = [{}]
        featureSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseCollections.FeatureSet',
                    'data': featureSet_obj,
                    'name': 'test_featureSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        featureSet_ref = str(featureSet_info[WSID_I])+'/'+str(featureSet_info[OBJID_I])+'/'+str(featureSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'featureSet_ref':      featureSet_ref,
                'file':                test_name+'-test_featureSet.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%genome_ref%%'+'.f:'+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }
        ret = self.getImpl().FeatureSetToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'][0])
        self.assertIsNotNone(ret['feature_ids_by_genome_ref'])
        self.assertNotEqual(len(ret['feature_ids_by_genome_ref'].keys()), 0)
        self.assertIsNotNone(ret['genome_ref_to_sci_name'])
        self.assertNotEqual(len(ret['genome_ref_to_sci_name'].keys()), 0)
        pass
        

    #### FeatureSetToFASTA_02()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_FeatureSetToFASTA_02()')
    def test_KBaseDataObjectToFileUtils_FeatureSetToFASTA_02(self):
        test_name = 'FeatureSetToFASTA_02'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        """
        reference_prok_genomes_WS = 'ReferenceDataManager'  # PROD and CI
        genome_ref_1 = reference_prok_genomes_WS+'/GCF_001566335.1/1'  # E. coli K-12 MG1655
        genome_ref_2 = reference_prok_genomes_WS+'/GCF_000021385.1/1'  # D. vulgaris str. 'Miyazaki F'
        genome_ref_3 = reference_prok_genomes_WS+'/GCF_900129775.1/1'  # Halobaculum gomorrense (16 contigs)
        """
        genomeInfo_0 = self.getGenomeInfo('GCF_001566335.1_ASM156633v1_genomic', 0)  # E. coli K-12 MG1655
        genomeInfo_1 = self.getGenomeInfo('GCF_000021385.1_ASM2138v1_genomic', 1)    # D. vulgaris str. 'Miyazaki F
        genomeInfo_2 = self.getGenomeInfo('GCF_001721825.1_ASM172182v1_genomic', 2)  # Pseudomonas aeruginosa
        genome_ref_0 = '/'.join([str(genomeInfo_0[WSID_I]),
                                 str(genomeInfo_0[OBJID_I]),
                                 str(genomeInfo_0[VERSION_I])])
        genome_ref_1 = '/'.join([str(genomeInfo_1[WSID_I]),
                                 str(genomeInfo_1[OBJID_I]),
                                 str(genomeInfo_1[VERSION_I])])
        genome_ref_2 = '/'.join([str(genomeInfo_2[WSID_I]),
                                 str(genomeInfo_2[OBJID_I]),
                                 str(genomeInfo_2[VERSION_I])])        
        

        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_0 = '/'.join([str(amaInfo_0[WSID_I]),
                              str(amaInfo_0[OBJID_I]),
                              str(amaInfo_0[VERSION_I])])

        genome_id_feature_id_delim = '.f:'
        feature_id_1 = 'AWN69_RS07145' # DnaA
        feature_id_2 = 'DVMF_RS00005'  # DnaA
        feature_id_3 = 'A6701_RS00005' # DnaA
        ama_feature_id_1 = '4_27'
        
        featureSet_obj = { 'description': 'test genomeSet',
                           'element_ordering': [
                               feature_id_1,
                               feature_id_2,
                               feature_id_3,
                               ama_feature_id_1
                           ],
                           'elements': { 
                               feature_id_1: [genome_ref_0],
                               feature_id_2: [genome_ref_1],
                               feature_id_3: [genome_ref_2],
                               ama_feature_id_1: [ama_ref_0]
                           }
                        }
        provenance = [{}]
        featureSet_info = self.getWsClient().save_objects({
            'workspace': self.getWsName(), 
            'objects': [
                {
                    'type': 'KBaseCollections.FeatureSet',
                    'data': featureSet_obj,
                    'name': 'test_featureSet',
                    'meta': {},
                    'provenance': provenance
                }
            ]})[0]

        featureSet_ref = str(featureSet_info[WSID_I])+'/'+str(featureSet_info[OBJID_I])+'/'+str(featureSet_info[VERSION_I])

        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'featureSet_ref':      featureSet_ref,
                'file':                test_name+'-test_featureSet.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'nuc',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%genome_ref%%'+'.f:'+'%%feature_id%%',
                'record_desc_pattern': '[%%genome_ref%%]',
                'case':                'upper',
                'linewrap':            50,
                'merge_fasta_files':   'TRUE'
                }
        ret = self.getImpl().FeatureSetToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'][0])
        self.assertIsNotNone(ret['feature_ids_by_genome_ref'])
        self.assertNotEqual(len(ret['feature_ids_by_genome_ref'].keys()), 0)
        self.assertIsNotNone(ret['genome_ref_to_sci_name'])
        self.assertNotEqual(len(ret['genome_ref_to_sci_name'].keys()), 0)
        pass
        

    #### AnnotatedMetagenomeAssemblyToFASTA_01()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_AnnotatedMetagenomeAssemblyToFASTA_01()')
    def test_KBaseDataObjectToFileUtils_AnnotatedMetagenomeAssemblyToFASTA_01(self):
        test_name = 'AnnotatedMetagenomeAssemblyToFASTA_01'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # upload test data
        ama_feature_cnt = 888

        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_0 = '/'.join([str(amaInfo_0[WSID_I]),
                              str(amaInfo_0[OBJID_I]),
                              str(amaInfo_0[VERSION_I])])
        
        # get protein fastas
        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'ama_ref':             ama_ref_0,
                'file':                test_name+'-test_ama.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'protein',
                'feature_type':        'CDS',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }
        ret = self.getImpl().AnnotatedMetagenomeAssemblyToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'])
        self.assertIsNotNone(ret['feature_ids'])
        self.assertNotEqual(len(ret['feature_ids']), 0)
        self.assertEqual(len(ret['feature_ids']), ama_feature_cnt)
        pass


    #### AnnotatedMetagenomeAssemblyToFASTA_02()
    ##
    @unittest.skip ('skipping test_KBaseDataObjectToFileUtils_AnnotatedMetagenomeAssemblyToFASTA_02()')
    def test_KBaseDataObjectToFileUtils_AnnotatedMetagenomeAssemblyToFASTA_02(self):
        test_name = 'AnnotatedMetagenomeAssemblyToFASTA_02'
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # upload test data
        ama_feature_cnt = 888

        amaInfo_0 = self.getAMAInfo("test_ama", 0)
        ama_ref_0 = '/'.join([str(amaInfo_0[WSID_I]),
                              str(amaInfo_0[OBJID_I]),
                              str(amaInfo_0[VERSION_I])])
        
        # get protein fastas
        output_dir = os.path.join(self.scratch,'fasta_out.'+str(uuid.uuid4()))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        parameters = {
                'ama_ref':             ama_ref_0,
                'file':                test_name+'-test_ama.fasta',
                'dir':                 output_dir,
                'console':             [],
                'invalid_msgs':        [],
                'residue_type':        'nuc',
                'feature_type':        'ALL',
                'record_id_pattern':   '%%feature_id%%',
                'record_desc_pattern': '[%%genome_id%%]',
                'case':                'upper',
                'linewrap':            50
                }
        ret = self.getImpl().AnnotatedMetagenomeAssemblyToFASTA(self.getContext(), parameters)[0]
        self.assertIsNotNone(ret['fasta_file_path'])
        self.assertIsNotNone(ret['feature_ids'])
        self.assertNotEqual(len(ret['feature_ids']), 0)
        self.assertEqual(len(ret['feature_ids']), ama_feature_cnt)
        pass
