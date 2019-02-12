
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
 * <p>Original spec-file type: ParseFastaStr_Params</p>
 * <pre>
 * ParseFastaStr() Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fasta_str",
    "residue_type",
    "case",
    "console",
    "invalid_msgs"
})
public class ParseFastaStrParams {

    @JsonProperty("fasta_str")
    private String fastaStr;
    @JsonProperty("residue_type")
    private String residueType;
    @JsonProperty("case")
    private String _case;
    @JsonProperty("console")
    private String console;
    @JsonProperty("invalid_msgs")
    private String invalidMsgs;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("fasta_str")
    public String getFastaStr() {
        return fastaStr;
    }

    @JsonProperty("fasta_str")
    public void setFastaStr(String fastaStr) {
        this.fastaStr = fastaStr;
    }

    public ParseFastaStrParams withFastaStr(String fastaStr) {
        this.fastaStr = fastaStr;
        return this;
    }

    @JsonProperty("residue_type")
    public String getResidueType() {
        return residueType;
    }

    @JsonProperty("residue_type")
    public void setResidueType(String residueType) {
        this.residueType = residueType;
    }

    public ParseFastaStrParams withResidueType(String residueType) {
        this.residueType = residueType;
        return this;
    }

    @JsonProperty("case")
    public String getCase() {
        return _case;
    }

    @JsonProperty("case")
    public void setCase(String _case) {
        this._case = _case;
    }

    public ParseFastaStrParams withCase(String _case) {
        this._case = _case;
        return this;
    }

    @JsonProperty("console")
    public String getConsole() {
        return console;
    }

    @JsonProperty("console")
    public void setConsole(String console) {
        this.console = console;
    }

    public ParseFastaStrParams withConsole(String console) {
        this.console = console;
        return this;
    }

    @JsonProperty("invalid_msgs")
    public String getInvalidMsgs() {
        return invalidMsgs;
    }

    @JsonProperty("invalid_msgs")
    public void setInvalidMsgs(String invalidMsgs) {
        this.invalidMsgs = invalidMsgs;
    }

    public ParseFastaStrParams withInvalidMsgs(String invalidMsgs) {
        this.invalidMsgs = invalidMsgs;
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
        return ((((((((((((("ParseFastaStrParams"+" [fastaStr=")+ fastaStr)+", residueType=")+ residueType)+", _case=")+ _case)+", console=")+ console)+", invalidMsgs=")+ invalidMsgs)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
