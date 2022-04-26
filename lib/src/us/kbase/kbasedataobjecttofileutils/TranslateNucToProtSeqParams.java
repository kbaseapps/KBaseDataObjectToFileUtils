
package us.kbase.kbasedataobjecttofileutils;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: TranslateNucToProtSeq_Params</p>
 * <pre>
 * TranslateNucToProtSeq() Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "nuc_seq",
    "genetic_code",
    "write_off_code_prot_seq"
})
public class TranslateNucToProtSeqParams {

    @JsonProperty("nuc_seq")
    private String nucSeq;
    @JsonProperty("genetic_code")
    private String geneticCode;
    @JsonProperty("write_off_code_prot_seq")
    private Long writeOffCodeProtSeq;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("nuc_seq")
    public String getNucSeq() {
        return nucSeq;
    }

    @JsonProperty("nuc_seq")
    public void setNucSeq(String nucSeq) {
        this.nucSeq = nucSeq;
    }

    public TranslateNucToProtSeqParams withNucSeq(String nucSeq) {
        this.nucSeq = nucSeq;
        return this;
    }

    @JsonProperty("genetic_code")
    public String getGeneticCode() {
        return geneticCode;
    }

    @JsonProperty("genetic_code")
    public void setGeneticCode(String geneticCode) {
        this.geneticCode = geneticCode;
    }

    public TranslateNucToProtSeqParams withGeneticCode(String geneticCode) {
        this.geneticCode = geneticCode;
        return this;
    }

    @JsonProperty("write_off_code_prot_seq")
    public Long getWriteOffCodeProtSeq() {
        return writeOffCodeProtSeq;
    }

    @JsonProperty("write_off_code_prot_seq")
    public void setWriteOffCodeProtSeq(Long writeOffCodeProtSeq) {
        this.writeOffCodeProtSeq = writeOffCodeProtSeq;
    }

    public TranslateNucToProtSeqParams withWriteOffCodeProtSeq(Long writeOffCodeProtSeq) {
        this.writeOffCodeProtSeq = writeOffCodeProtSeq;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((("TranslateNucToProtSeqParams"+" [nucSeq=")+ nucSeq)+", geneticCode=")+ geneticCode)+", writeOffCodeProtSeq=")+ writeOffCodeProtSeq)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
