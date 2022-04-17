
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
 * <p>Original spec-file type: AnnotatedMetagenomeAssemblyToFASTA_Params</p>
 * <pre>
 * AnnotatedMetagenomeAssemblyToFASTA() Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "ama_ref",
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
    "id_len_limit"
})
public class AnnotatedMetagenomeAssemblyToFASTAParams {

    @JsonProperty("ama_ref")
    private java.lang.String amaRef;
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
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("ama_ref")
    public java.lang.String getAmaRef() {
        return amaRef;
    }

    @JsonProperty("ama_ref")
    public void setAmaRef(java.lang.String amaRef) {
        this.amaRef = amaRef;
    }

    public AnnotatedMetagenomeAssemblyToFASTAParams withAmaRef(java.lang.String amaRef) {
        this.amaRef = amaRef;
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withFile(java.lang.String file) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withDir(java.lang.String dir) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withConsole(List<String> console) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withInvalidMsgs(List<String> invalidMsgs) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withResidueType(java.lang.String residueType) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withFeatureType(java.lang.String featureType) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withRecordIdPattern(java.lang.String recordIdPattern) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withRecordDescPattern(java.lang.String recordDescPattern) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withCase(java.lang.String _case) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withLinewrap(Long linewrap) {
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

    public AnnotatedMetagenomeAssemblyToFASTAParams withIdLenLimit(Long idLenLimit) {
        this.idLenLimit = idLenLimit;
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
        return ((((((((((((((((((((((((((("AnnotatedMetagenomeAssemblyToFASTAParams"+" [amaRef=")+ amaRef)+", file=")+ file)+", dir=")+ dir)+", console=")+ console)+", invalidMsgs=")+ invalidMsgs)+", residueType=")+ residueType)+", featureType=")+ featureType)+", recordIdPattern=")+ recordIdPattern)+", recordDescPattern=")+ recordDescPattern)+", _case=")+ _case)+", linewrap=")+ linewrap)+", idLenLimit=")+ idLenLimit)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
