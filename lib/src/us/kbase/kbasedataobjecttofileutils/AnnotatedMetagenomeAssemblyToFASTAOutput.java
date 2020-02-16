
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
 * <p>Original spec-file type: AnnotatedMetagenomeAssemblyToFASTA_Output</p>
 * <pre>
 * AnnotatedMetagenomeAssemblyToFASTA() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file_path",
    "feature_ids"
})
public class AnnotatedMetagenomeAssemblyToFASTAOutput {

    @JsonProperty("fasta_file_path")
    private java.lang.String fastaFilePath;
    @JsonProperty("feature_ids")
    private List<String> featureIds;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fasta_file_path")
    public java.lang.String getFastaFilePath() {
        return fastaFilePath;
    }

    @JsonProperty("fasta_file_path")
    public void setFastaFilePath(java.lang.String fastaFilePath) {
        this.fastaFilePath = fastaFilePath;
    }

    public AnnotatedMetagenomeAssemblyToFASTAOutput withFastaFilePath(java.lang.String fastaFilePath) {
        this.fastaFilePath = fastaFilePath;
        return this;
    }

    @JsonProperty("feature_ids")
    public List<String> getFeatureIds() {
        return featureIds;
    }

    @JsonProperty("feature_ids")
    public void setFeatureIds(List<String> featureIds) {
        this.featureIds = featureIds;
    }

    public AnnotatedMetagenomeAssemblyToFASTAOutput withFeatureIds(List<String> featureIds) {
        this.featureIds = featureIds;
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
        return ((((((("AnnotatedMetagenomeAssemblyToFASTAOutput"+" [fastaFilePath=")+ fastaFilePath)+", featureIds=")+ featureIds)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
