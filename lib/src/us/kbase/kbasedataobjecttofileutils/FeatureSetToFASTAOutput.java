
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
 * <p>Original spec-file type: FeatureSetToFASTA_Output</p>
 * <pre>
 * FeatureSetToFASTA() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file_path",
    "feature_ids_by_genome_ref",
    "feature_id_to_function",
    "genome_ref_to_sci_name"
})
public class FeatureSetToFASTAOutput {

    @JsonProperty("fasta_file_path")
    private java.lang.String fastaFilePath;
    @JsonProperty("feature_ids_by_genome_ref")
    private Map<String, List<String>> featureIdsByGenomeRef;
    @JsonProperty("feature_id_to_function")
    private Map<String, String> featureIdToFunction;
    @JsonProperty("genome_ref_to_sci_name")
    private Map<String, String> genomeRefToSciName;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fasta_file_path")
    public java.lang.String getFastaFilePath() {
        return fastaFilePath;
    }

    @JsonProperty("fasta_file_path")
    public void setFastaFilePath(java.lang.String fastaFilePath) {
        this.fastaFilePath = fastaFilePath;
    }

    public FeatureSetToFASTAOutput withFastaFilePath(java.lang.String fastaFilePath) {
        this.fastaFilePath = fastaFilePath;
        return this;
    }

    @JsonProperty("feature_ids_by_genome_ref")
    public Map<String, List<String>> getFeatureIdsByGenomeRef() {
        return featureIdsByGenomeRef;
    }

    @JsonProperty("feature_ids_by_genome_ref")
    public void setFeatureIdsByGenomeRef(Map<String, List<String>> featureIdsByGenomeRef) {
        this.featureIdsByGenomeRef = featureIdsByGenomeRef;
    }

    public FeatureSetToFASTAOutput withFeatureIdsByGenomeRef(Map<String, List<String>> featureIdsByGenomeRef) {
        this.featureIdsByGenomeRef = featureIdsByGenomeRef;
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

    public FeatureSetToFASTAOutput withFeatureIdToFunction(Map<String, String> featureIdToFunction) {
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

    public FeatureSetToFASTAOutput withGenomeRefToSciName(Map<String, String> genomeRefToSciName) {
        this.genomeRefToSciName = genomeRefToSciName;
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
        return ((((((((((("FeatureSetToFASTAOutput"+" [fastaFilePath=")+ fastaFilePath)+", featureIdsByGenomeRef=")+ featureIdsByGenomeRef)+", featureIdToFunction=")+ featureIdToFunction)+", genomeRefToSciName=")+ genomeRefToSciName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
