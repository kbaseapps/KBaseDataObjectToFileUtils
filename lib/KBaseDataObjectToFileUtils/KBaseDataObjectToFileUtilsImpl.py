# -*- coding: utf-8 -*-
#BEGIN_HEADER
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
    GIT_COMMIT_HASH = "700a94610e9910fcac6758c91999eb5713901032"
    
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass
    

    def GenomeAnnotationToFASTA(self, ctx, params):
        """
        Methods for converting KBase Data Objects to common bioinformatics format files
        **
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
