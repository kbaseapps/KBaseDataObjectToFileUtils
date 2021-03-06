
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
 * <p>Original spec-file type: FeatureSetToFASTA_Params</p>
 * <pre>
 * FeatureSetToFASTA() Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "featureSet_ref",
    "file",
    "dir",
    "console",
    "invalid_msgs",
    "residue_type",
    "feature_type",
    "record_id_pattern",
    "record_desc_pattern",
    "case",
    "linewrap"
})
public class FeatureSetToFASTAParams {

    @JsonProperty("featureSet_ref")
    private java.lang.String featureSetRef;
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
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("featureSet_ref")
    public java.lang.String getFeatureSetRef() {
        return featureSetRef;
    }

    @JsonProperty("featureSet_ref")
    public void setFeatureSetRef(java.lang.String featureSetRef) {
        this.featureSetRef = featureSetRef;
    }

    public FeatureSetToFASTAParams withFeatureSetRef(java.lang.String featureSetRef) {
        this.featureSetRef = featureSetRef;
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

    public FeatureSetToFASTAParams withFile(java.lang.String file) {
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

    public FeatureSetToFASTAParams withDir(java.lang.String dir) {
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

    public FeatureSetToFASTAParams withConsole(List<String> console) {
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

    public FeatureSetToFASTAParams withInvalidMsgs(List<String> invalidMsgs) {
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

    public FeatureSetToFASTAParams withResidueType(java.lang.String residueType) {
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

    public FeatureSetToFASTAParams withFeatureType(java.lang.String featureType) {
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

    public FeatureSetToFASTAParams withRecordIdPattern(java.lang.String recordIdPattern) {
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

    public FeatureSetToFASTAParams withRecordDescPattern(java.lang.String recordDescPattern) {
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

    public FeatureSetToFASTAParams withCase(java.lang.String _case) {
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

    public FeatureSetToFASTAParams withLinewrap(Long linewrap) {
        this.linewrap = linewrap;
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
        return ((((((((((((((((((((((((("FeatureSetToFASTAParams"+" [featureSetRef=")+ featureSetRef)+", file=")+ file)+", dir=")+ dir)+", console=")+ console)+", invalidMsgs=")+ invalidMsgs)+", residueType=")+ residueType)+", featureType=")+ featureType)+", recordIdPattern=")+ recordIdPattern)+", recordDescPattern=")+ recordDescPattern)+", _case=")+ _case)+", linewrap=")+ linewrap)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
