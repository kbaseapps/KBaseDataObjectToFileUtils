
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
 * <p>Original spec-file type: TranslateNucToProtSeq_Output</p>
 * <pre>
 * TranslateNucToProtSeq() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "prot_seq"
})
public class TranslateNucToProtSeqOutput {

    @JsonProperty("prot_seq")
    private String protSeq;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("prot_seq")
    public String getProtSeq() {
        return protSeq;
    }

    @JsonProperty("prot_seq")
    public void setProtSeq(String protSeq) {
        this.protSeq = protSeq;
    }

    public TranslateNucToProtSeqOutput withProtSeq(String protSeq) {
        this.protSeq = protSeq;
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
        return ((((("TranslateNucToProtSeqOutput"+" [protSeq=")+ protSeq)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
