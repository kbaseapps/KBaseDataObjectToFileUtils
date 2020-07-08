
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
    "feature_ids_by_genome_id",
    "feature_id_to_function",
    "genome_ref_to_sci_name",
    "genome_ref_to_obj_name"
})
public class GenomeSetToFASTAOutput {

    @JsonProperty("fasta_file_path_list")
    private List<String> fastaFilePathList;
    @JsonProperty("feature_ids_by_genome_id")
    private Map<String, List<String>> featureIdsByGenomeId;
    @JsonProperty("feature_id_to_function")
    private Map<String, String> featureIdToFunction;
    @JsonProperty("genome_ref_to_sci_name")
    private Map<String, String> genomeRefToSciName;
    @JsonProperty("genome_ref_to_obj_name")
    private Map<String, String> genomeRefToObjName;
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

    @JsonProperty("feature_id_to_function")
    public Map<String, String> getFeatureIdToFunction() {
        return featureIdToFunction;
    }

    @JsonProperty("feature_id_to_function")
    public void setFeatureIdToFunction(Map<String, String> featureIdToFunction) {
        this.featureIdToFunction = featureIdToFunction;
    }

    public GenomeSetToFASTAOutput withFeatureIdToFunction(Map<String, String> featureIdToFunction) {
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

    public GenomeSetToFASTAOutput withGenomeRefToSciName(Map<String, String> genomeRefToSciName) {
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

    public GenomeSetToFASTAOutput withGenomeRefToObjName(Map<String, String> genomeRefToObjName) {
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
        return ((((((((((((("GenomeSetToFASTAOutput"+" [fastaFilePathList=")+ fastaFilePathList)+", featureIdsByGenomeId=")+ featureIdsByGenomeId)+", featureIdToFunction=")+ featureIdToFunction)+", genomeRefToSciName=")+ genomeRefToSciName)+", genomeRefToObjName=")+ genomeRefToObjName)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
