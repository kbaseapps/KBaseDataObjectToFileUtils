/*
** A KBase module: kb_blast
**
** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
** 
*/

module KBaseDataObjectToFileUtils {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;
    typedef string feature_id;
    typedef string genome_id;
    typedef string path_type;
    typedef string pattern_type;
    typedef string log_msg;
    typedef string true_false;


    /* TranslateNucToProtSeq() Params
    */
    typedef structure {
	string  nuc_seq;
	string  genetic_code;
    } TranslateNucToProtSeq_Params;


    /* TranslateNucToProtSeq() Output
    */
    typedef structure {
	string  prot_seq;
    } TranslateNucToProtSeq_Output;


    /* GenomeAnnotationToFASTA() Params
    */
    typedef structure {
	data_obj_ref   genome_ref;
	path_type      file;
	path_type      dir;
	list<log_msg>  console;
	list<log_msg>  invalid_msgs;
	string         residue_type;
	string         feature_type;
	pattern_type   record_id_pattern;
	pattern_type   record_desc_pattern;
	string         case;
	int            linewrap;
    } GenomeAnnotationToFASTA_Params;


    /* GenomeAnnotationToFASTA() Output
    */
    typedef structure {
	path_type         fasta_file_path;
	list<feature_id>  feature_ids;
    } GenomeAnnotationToFASTA_Output;


    /* GenomeSetToFASTA() Params
    */
    typedef structure {
	data_obj_ref   genomeSet_ref;
	path_type      file;
	path_type      dir;
	list<log_msg>  console;
	list<log_msg>  invalid_msgs;
	string         residue_type;
	string         feature_type;
	pattern_type   record_id_pattern;
	pattern_type   record_desc_pattern;
	string         case;
	int            linewrap;
	true_false     merge_fasta_files;
    } GenomeSetToFASTA_Params;


    /* GenomeSetToFASTA() Output
    */
    typedef structure {
	list<path_type>                       fasta_file_path_list;
	mapping<genome_id, list<feature_id>>  feature_ids_by_genome_id;
    } GenomeSetToFASTA_Output;


    /* FeatureSetToFASTA() Params
    */
    typedef structure {
	data_obj_ref   featureSet_ref;
	path_type      file;
	path_type      dir;
	list<log_msg>  console;
	list<log_msg>  invalid_msgs;
	string         residue_type;
	string         feature_type;
	pattern_type   record_id_pattern;
	pattern_type   record_desc_pattern;
	string         case;
	int            linewrap;
    } FeatureSetToFASTA_Params;


    /* FeatureSetToFASTA() Output
    */
    typedef structure {
	path_type                             fasta_file_path;
	mapping<genome_id, list<feature_id>>  feature_ids_by_genome_id;
    } FeatureSetToFASTA_Output;


    /*  Methods for converting KBase Data Objects to common bioinformatics format files
    **
    */
    funcdef TranslateNucToProtSeq (TranslateNucToProtSeq_Params params)  returns (TranslateNucToProtSeq_Output) authentication required;

    /* this should not be used, but is temporarily being retained to compare speed */
    funcdef GenomeToFASTA (GenomeAnnotationToFASTA_Params params)  returns (GenomeAnnotationToFASTA_Output) authentication required;

    funcdef GenomeAnnotationToFASTA (GenomeAnnotationToFASTA_Params params)  returns (GenomeAnnotationToFASTA_Output) authentication required;

    funcdef GenomeSetToFASTA (GenomeSetToFASTA_Params params)  returns (GenomeSetToFASTA_Output) authentication required;

    funcdef FeatureSetToFASTA (FeatureSetToFASTA_Params params)  returns (FeatureSetToFASTA_Output) authentication required;
};
