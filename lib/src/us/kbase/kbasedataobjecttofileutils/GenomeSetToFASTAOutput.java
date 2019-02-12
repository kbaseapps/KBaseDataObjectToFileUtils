
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
 * <p>Original spec-file type: GenomeSetToFASTA_Output</p>
 * <pre>
 * GenomeSetToFASTA() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_file_path_list",
    "feature_ids_by_genome_id"
})
public class GenomeSetToFASTAOutput {

    @JsonProperty("fasta_file_path_list")
    private List<String> fastaFilePathList;
    @JsonProperty("feature_ids_by_genome_id")
    private Map<String, List<String>> featureIdsByGenomeId;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fasta_file_path_list")
    public List<String> getFastaFilePathList() {
        return fastaFilePathList;
    }

    @JsonProperty("fasta_file_path_list")
    public void setFastaFilePathList(List<String> fastaFilePathList) {
        this.fastaFilePathList = fastaFilePathList;
    }

    public GenomeSetToFASTAOutput withFastaFilePathList(List<String> fastaFilePathList) {
        this.fastaFilePathList = fastaFilePathList;
        return this;
    }

    @JsonProperty("feature_ids_by_genome_id")
    public Map<String, List<String>> getFeatureIdsByGenomeId() {
        return featureIdsByGenomeId;
    }

    @JsonProperty("feature_ids_by_genome_id")
    public void setFeatureIdsByGenomeId(Map<String, List<String>> featureIdsByGenomeId) {
        this.featureIdsByGenomeId = featureIdsByGenomeId;
    }

    public GenomeSetToFASTAOutput withFeatureIdsByGenomeId(Map<String, List<String>> featureIdsByGenomeId) {
        this.featureIdsByGenomeId = featureIdsByGenomeId;
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
        return ((((((("GenomeSetToFASTAOutput"+" [fastaFilePathList=")+ fastaFilePathList)+", featureIdsByGenomeId=")+ featureIdsByGenomeId)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
