
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
 * <p>Original spec-file type: ParseFastaStr_Output</p>
 * <pre>
 * ParseFastaStr() Output
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "id",
    "desc",
    "seq"
})
public class ParseFastaStrOutput {

    @JsonProperty("id")
    private String id;
    @JsonProperty("desc")
    private String desc;
    @JsonProperty("seq")
    private String seq;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("id")
    public String getId() {
        return id;
    }

    @JsonProperty("id")
    public void setId(String id) {
        this.id = id;
    }

    public ParseFastaStrOutput withId(String id) {
        this.id = id;
        return this;
    }

    @JsonProperty("desc")
    public String getDesc() {
        return desc;
    }

    @JsonProperty("desc")
    public void setDesc(String desc) {
        this.desc = desc;
    }

    public ParseFastaStrOutput withDesc(String desc) {
        this.desc = desc;
        return this;
    }

    @JsonProperty("seq")
    public String getSeq() {
        return seq;
    }

    @JsonProperty("seq")
    public void setSeq(String seq) {
        this.seq = seq;
    }

    public ParseFastaStrOutput withSeq(String seq) {
        this.seq = seq;
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
        return ((((((((("ParseFastaStrOutput"+" [id=")+ id)+", desc=")+ desc)+", seq=")+ seq)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
