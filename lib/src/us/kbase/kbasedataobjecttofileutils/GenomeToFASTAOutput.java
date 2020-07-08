
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
 * <p>Original spec-file type: GenomeToFASTA_Output</p>
 * <pre>
 * GenomeToFASTA() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file_path",
    "feature_ids",
    "feature_id_to_function",
    "genome_ref_to_sci_name",
    "genome_ref_to_obj_name"
})
public class GenomeToFASTAOutput {

    @JsonProperty("fasta_file_path")
    private java.lang.String fastaFilePath;
    @JsonProperty("feature_ids")
    private List<String> featureIds;
    @JsonProperty("feature_id_to_function")
    private Map<String, String> featureIdToFunction;
    @JsonProperty("genome_ref_to_sci_name")
    private Map<String, String> genomeRefToSciName;
    @JsonProperty("genome_ref_to_obj_name")
    private Map<String, String> genomeRefToObjName;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fasta_file_path")
    public java.lang.String getFastaFilePath() {
        return fastaFilePath;
    }

    @JsonProperty("fasta_file_path")
    public void setFastaFilePath(java.lang.String fastaFilePath) {
        this.fastaFilePath = fastaFilePath;
    }

    public GenomeToFASTAOutput withFastaFilePath(java.lang.String fastaFilePath) {
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

    public GenomeToFASTAOutput withFeatureIds(List<String> featureIds) {
        this.featureIds = featureIds;
        return this;
    }

    @JsonProperty("feature_id_to_function")
    public Map<String, String> getFeatureIdToFunction() {
        return featureIdToFunction;
    }

    @JsonProperty("feature_id_to_function")
    public void setFeatureIdToFunction(Map<String, String> featureIdToFunction) {
        this.featureIdToFunction = featureIdToFunction;
    }

    public GenomeToFASTAOutput withFeatureIdToFunction(Map<String, String> featureIdToFunction) {
        this.featureIdToFunction = featureIdToFunction;
        return this;
    }

    @JsonProperty("genome_ref_to_sci_name")
    public Map<String, String> getGenomeRefToSciName() {
        return genomeRefToSciName;
    }

    @JsonProperty("genome_ref_to_sci_name")
    public void setGenomeRefToSciName(Map<String, String> genomeRefToSciName) {
        this.genomeRefToSciName = genomeRefToSciName;
    }

    public GenomeToFASTAOutput withGenomeRefToSciName(Map<String, String> genomeRefToSciName) {
        this.genomeRefToSciName = genomeRefToSciName;
        return this;
    }

    @JsonProperty("genome_ref_to_obj_name")
    public Map<String, String> getGenomeRefToObjName() {
        return genomeRefToObjName;
    }

    @JsonProperty("genome_ref_to_obj_name")
    public void setGenomeRefToObjName(Map<String, String> genomeRefToObjName) {
        this.genomeRefToObjName = genomeRefToObjName;
    }

    public GenomeToFASTAOutput withGenomeRefToObjName(Map<String, String> genomeRefToObjName) {
        this.genomeRefToObjName = genomeRefToObjName;
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
        return ((((((((((((("GenomeToFASTAOutput"+" [fastaFilePath=")+ fastaFilePath)+", featureIds=")+ featureIds)+", featureIdToFunction=")+ featureIdToFunction)+", genomeRefToSciName=")+ genomeRefToSciName)+", genomeRefToObjName=")+ genomeRefToObjName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
