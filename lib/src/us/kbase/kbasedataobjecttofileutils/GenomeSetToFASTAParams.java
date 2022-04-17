
package us.kbase.kbasedataobjecttofileutils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: GenomeSetToFASTA_Params</p>
 * <pre>
 * GenomeSetToFASTA() Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "genomeSet_ref",
    "file",
    "dir",
    "console",
    "invalid_msgs",
    "residue_type",
    "feature_type",
    "record_id_pattern",
    "record_desc_pattern",
    "case",
    "linewrap",
    "id_len_limit",
    "merge_fasta_files"
})
public class GenomeSetToFASTAParams {

    @JsonProperty("genomeSet_ref")
    private java.lang.String genomeSetRef;
    @JsonProperty("file")
    private java.lang.String file;
    @JsonProperty("dir")
    private java.lang.String dir;
    @JsonProperty("console")
    private List<String> console;
    @JsonProperty("invalid_msgs")
    private List<String> invalidMsgs;
    @JsonProperty("residue_type")
    private java.lang.String residueType;
    @JsonProperty("feature_type")
    private java.lang.String featureType;
    @JsonProperty("record_id_pattern")
    private java.lang.String recordIdPattern;
    @JsonProperty("record_desc_pattern")
    private java.lang.String recordDescPattern;
    @JsonProperty("case")
    private java.lang.String _case;
    @JsonProperty("linewrap")
    private Long linewrap;
    @JsonProperty("id_len_limit")
    private Long idLenLimit;
    @JsonProperty("merge_fasta_files")
    private java.lang.String mergeFastaFiles;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("genomeSet_ref")
    public java.lang.String getGenomeSetRef() {
        return genomeSetRef;
    }

    @JsonProperty("genomeSet_ref")
    public void setGenomeSetRef(java.lang.String genomeSetRef) {
        this.genomeSetRef = genomeSetRef;
    }

    public GenomeSetToFASTAParams withGenomeSetRef(java.lang.String genomeSetRef) {
        this.genomeSetRef = genomeSetRef;
        return this;
    }

    @JsonProperty("file")
    public java.lang.String getFile() {
        return file;
    }

    @JsonProperty("file")
    public void setFile(java.lang.String file) {
        this.file = file;
    }

    public GenomeSetToFASTAParams withFile(java.lang.String file) {
        this.file = file;
        return this;
    }

    @JsonProperty("dir")
    public java.lang.String getDir() {
        return dir;
    }

    @JsonProperty("dir")
    public void setDir(java.lang.String dir) {
        this.dir = dir;
    }

    public GenomeSetToFASTAParams withDir(java.lang.String dir) {
        this.dir = dir;
        return this;
    }

    @JsonProperty("console")
    public List<String> getConsole() {
        return console;
    }

    @JsonProperty("console")
    public void setConsole(List<String> console) {
        this.console = console;
    }

    public GenomeSetToFASTAParams withConsole(List<String> console) {
        this.console = console;
        return this;
    }

    @JsonProperty("invalid_msgs")
    public List<String> getInvalidMsgs() {
        return invalidMsgs;
    }

    @JsonProperty("invalid_msgs")
    public void setInvalidMsgs(List<String> invalidMsgs) {
        this.invalidMsgs = invalidMsgs;
    }

    public GenomeSetToFASTAParams withInvalidMsgs(List<String> invalidMsgs) {
        this.invalidMsgs = invalidMsgs;
        return this;
    }

    @JsonProperty("residue_type")
    public java.lang.String getResidueType() {
        return residueType;
    }

    @JsonProperty("residue_type")
    public void setResidueType(java.lang.String residueType) {
        this.residueType = residueType;
    }

    public GenomeSetToFASTAParams withResidueType(java.lang.String residueType) {
        this.residueType = residueType;
        return this;
    }

    @JsonProperty("feature_type")
    public java.lang.String getFeatureType() {
        return featureType;
    }

    @JsonProperty("feature_type")
    public void setFeatureType(java.lang.String featureType) {
        this.featureType = featureType;
    }

    public GenomeSetToFASTAParams withFeatureType(java.lang.String featureType) {
        this.featureType = featureType;
        return this;
    }

    @JsonProperty("record_id_pattern")
    public java.lang.String getRecordIdPattern() {
        return recordIdPattern;
    }

    @JsonProperty("record_id_pattern")
    public void setRecordIdPattern(java.lang.String recordIdPattern) {
        this.recordIdPattern = recordIdPattern;
    }

    public GenomeSetToFASTAParams withRecordIdPattern(java.lang.String recordIdPattern) {
        this.recordIdPattern = recordIdPattern;
        return this;
    }

    @JsonProperty("record_desc_pattern")
    public java.lang.String getRecordDescPattern() {
        return recordDescPattern;
    }

    @JsonProperty("record_desc_pattern")
    public void setRecordDescPattern(java.lang.String recordDescPattern) {
        this.recordDescPattern = recordDescPattern;
    }

    public GenomeSetToFASTAParams withRecordDescPattern(java.lang.String recordDescPattern) {
        this.recordDescPattern = recordDescPattern;
        return this;
    }

    @JsonProperty("case")
    public java.lang.String getCase() {
        return _case;
    }

    @JsonProperty("case")
    public void setCase(java.lang.String _case) {
        this._case = _case;
    }

    public GenomeSetToFASTAParams withCase(java.lang.String _case) {
        this._case = _case;
        return this;
    }

    @JsonProperty("linewrap")
    public Long getLinewrap() {
        return linewrap;
    }

    @JsonProperty("linewrap")
    public void setLinewrap(Long linewrap) {
        this.linewrap = linewrap;
    }

    public GenomeSetToFASTAParams withLinewrap(Long linewrap) {
        this.linewrap = linewrap;
        return this;
    }

    @JsonProperty("id_len_limit")
    public Long getIdLenLimit() {
        return idLenLimit;
    }

    @JsonProperty("id_len_limit")
    public void setIdLenLimit(Long idLenLimit) {
        this.idLenLimit = idLenLimit;
    }

    public GenomeSetToFASTAParams withIdLenLimit(Long idLenLimit) {
        this.idLenLimit = idLenLimit;
        return this;
    }

    @JsonProperty("merge_fasta_files")
    public java.lang.String getMergeFastaFiles() {
        return mergeFastaFiles;
    }

    @JsonProperty("merge_fasta_files")
    public void setMergeFastaFiles(java.lang.String mergeFastaFiles) {
        this.mergeFastaFiles = mergeFastaFiles;
    }

    public GenomeSetToFASTAParams withMergeFastaFiles(java.lang.String mergeFastaFiles) {
        this.mergeFastaFiles = mergeFastaFiles;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((((((("GenomeSetToFASTAParams"+" [genomeSetRef=")+ genomeSetRef)+", file=")+ file)+", dir=")+ dir)+", console=")+ console)+", invalidMsgs=")+ invalidMsgs)+", residueType=")+ residueType)+", featureType=")+ featureType)+", recordIdPattern=")+ recordIdPattern)+", recordDescPattern=")+ recordDescPattern)+", _case=")+ _case)+", linewrap=")+ linewrap)+", idLenLimit=")+ idLenLimit)+", mergeFastaFiles=")+ mergeFastaFiles)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
