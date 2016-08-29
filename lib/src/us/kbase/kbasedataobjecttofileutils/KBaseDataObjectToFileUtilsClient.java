package us.kbase.kbasedataobjecttofileutils;

import com.fasterxml.jackson.core.type.TypeReference;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import us.kbase.auth.AuthToken;
import us.kbase.common.service.JsonClientCaller;
import us.kbase.common.service.JsonClientException;
import us.kbase.common.service.RpcContext;
import us.kbase.common.service.UnauthorizedException;

/**
 * <p>Original spec-file module name: KBaseDataObjectToFileUtils</p>
 * <pre>
 * ** A KBase module: kb_blast
 * **
 * ** This module contains methods for converting KBase Data Objects to common bioinformatics file formats
 * **
 * </pre>
 */
public class KBaseDataObjectToFileUtilsClient {
    private JsonClientCaller caller;
    private String serviceVersion = null;


    /** Constructs a client with a custom URL and no user credentials.
     * @param url the URL of the service.
     */
    public KBaseDataObjectToFileUtilsClient(URL url) {
        caller = new JsonClientCaller(url);
    }
    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param token the user's authorization token.
     * @throws UnauthorizedException if the token is not valid.
     * @throws IOException if an IOException occurs when checking the token's
     * validity.
     */
    public KBaseDataObjectToFileUtilsClient(URL url, AuthToken token) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, token);
    }

    /** Constructs a client with a custom URL
     * and a custom authorization service URL.
     * @param url the URL of the service.
     * @param token the user's authorization token.
     * @param auth the URL of the authorization server.
     * @throws UnauthorizedException if the token is not valid.
     * @throws IOException if an IOException occurs when checking the token's
     * validity.
     */
    public KBaseDataObjectToFileUtilsClient(URL url, AuthToken token, URL auth) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, token, auth);
    }

    /** Constructs a client with a custom URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KBaseDataObjectToFileUtilsClient(URL url, String user, String password) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password);
    }

    /** Constructs a client with a custom URL
     * and a custom authorization service URL.
     * @param url the URL of the service.
     * @param user the user name.
     * @param password the password for the user name.
     * @param auth the URL of the authorization server.
     * @throws UnauthorizedException if the credentials are not valid.
     * @throws IOException if an IOException occurs when checking the user's
     * credentials.
     */
    public KBaseDataObjectToFileUtilsClient(URL url, String user, String password, URL auth) throws UnauthorizedException, IOException {
        caller = new JsonClientCaller(url, user, password, auth);
    }

    /** Get the token this client uses to communicate with the server.
     * @return the authorization token.
     */
    public AuthToken getToken() {
        return caller.getToken();
    }

    /** Get the URL of the service with which this client communicates.
     * @return the service URL.
     */
    public URL getURL() {
        return caller.getURL();
    }

    /** Set the timeout between establishing a connection to a server and
     * receiving a response. A value of zero or null implies no timeout.
     * @param milliseconds the milliseconds to wait before timing out when
     * attempting to read from a server.
     */
    public void setConnectionReadTimeOut(Integer milliseconds) {
        this.caller.setConnectionReadTimeOut(milliseconds);
    }

    /** Check if this client allows insecure http (vs https) connections.
     * @return true if insecure connections are allowed.
     */
    public boolean isInsecureHttpConnectionAllowed() {
        return caller.isInsecureHttpConnectionAllowed();
    }

    /** Deprecated. Use isInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public boolean isAuthAllowedForHttp() {
        return caller.isAuthAllowedForHttp();
    }

    /** Set whether insecure http (vs https) connections should be allowed by
     * this client.
     * @param allowed true to allow insecure connections. Default false
     */
    public void setIsInsecureHttpConnectionAllowed(boolean allowed) {
        caller.setInsecureHttpConnectionAllowed(allowed);
    }

    /** Deprecated. Use setIsInsecureHttpConnectionAllowed().
     * @deprecated
     */
    public void setAuthAllowedForHttp(boolean isAuthAllowedForHttp) {
        caller.setAuthAllowedForHttp(isAuthAllowedForHttp);
    }

    /** Set whether all SSL certificates, including self-signed certificates,
     * should be trusted.
     * @param trustAll true to trust all certificates. Default false.
     */
    public void setAllSSLCertificatesTrusted(final boolean trustAll) {
        caller.setAllSSLCertificatesTrusted(trustAll);
    }
    
    /** Check if this client trusts all SSL certificates, including
     * self-signed certificates.
     * @return true if all certificates are trusted.
     */
    public boolean isAllSSLCertificatesTrusted() {
        return caller.isAllSSLCertificatesTrusted();
    }
    /** Sets streaming mode on. In this case, the data will be streamed to
     * the server in chunks as it is read from disk rather than buffered in
     * memory. Many servers are not compatible with this feature.
     * @param streamRequest true to set streaming mode on, false otherwise.
     */
    public void setStreamingModeOn(boolean streamRequest) {
        caller.setStreamingModeOn(streamRequest);
    }

    /** Returns true if streaming mode is on.
     * @return true if streaming mode is on.
     */
    public boolean isStreamingModeOn() {
        return caller.isStreamingModeOn();
    }

    public void _setFileForNextRpcResponse(File f) {
        caller.setFileForNextRpcResponse(f);
    }

    public String getServiceVersion() {
        return this.serviceVersion;
    }

    public void setServiceVersion(String newValue) {
        this.serviceVersion = newValue;
    }

    /**
     * <p>Original spec-file function name: TranslateNucToProtSeq</p>
     * <pre>
     * Methods for converting KBase Data Objects to common bioinformatics format files
     * **
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbasedataobjecttofileutils.TranslateNucToProtSeqParams TranslateNucToProtSeqParams} (original type "TranslateNucToProtSeq_Params")
     * @return   instance of type {@link us.kbase.kbasedataobjecttofileutils.TranslateNucToProtSeqOutput TranslateNucToProtSeqOutput} (original type "TranslateNucToProtSeq_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public TranslateNucToProtSeqOutput translateNucToProtSeq(TranslateNucToProtSeqParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<TranslateNucToProtSeqOutput>> retType = new TypeReference<List<TranslateNucToProtSeqOutput>>() {};
        List<TranslateNucToProtSeqOutput> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.TranslateNucToProtSeq", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: GenomeToFASTA</p>
     * <pre>
     * this should not be used, but is temporarily being retained to compare speed
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeAnnotationToFASTAParams GenomeAnnotationToFASTAParams} (original type "GenomeAnnotationToFASTA_Params")
     * @return   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeAnnotationToFASTAOutput GenomeAnnotationToFASTAOutput} (original type "GenomeAnnotationToFASTA_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public GenomeAnnotationToFASTAOutput genomeToFASTA(GenomeAnnotationToFASTAParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<GenomeAnnotationToFASTAOutput>> retType = new TypeReference<List<GenomeAnnotationToFASTAOutput>>() {};
        List<GenomeAnnotationToFASTAOutput> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.GenomeToFASTA", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: GenomeAnnotationToFASTA</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeAnnotationToFASTAParams GenomeAnnotationToFASTAParams} (original type "GenomeAnnotationToFASTA_Params")
     * @return   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeAnnotationToFASTAOutput GenomeAnnotationToFASTAOutput} (original type "GenomeAnnotationToFASTA_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public GenomeAnnotationToFASTAOutput genomeAnnotationToFASTA(GenomeAnnotationToFASTAParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<GenomeAnnotationToFASTAOutput>> retType = new TypeReference<List<GenomeAnnotationToFASTAOutput>>() {};
        List<GenomeAnnotationToFASTAOutput> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.GenomeAnnotationToFASTA", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: GenomeSetToFASTA</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeSetToFASTAParams GenomeSetToFASTAParams} (original type "GenomeSetToFASTA_Params")
     * @return   instance of type {@link us.kbase.kbasedataobjecttofileutils.GenomeSetToFASTAOutput GenomeSetToFASTAOutput} (original type "GenomeSetToFASTA_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public GenomeSetToFASTAOutput genomeSetToFASTA(GenomeSetToFASTAParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<GenomeSetToFASTAOutput>> retType = new TypeReference<List<GenomeSetToFASTAOutput>>() {};
        List<GenomeSetToFASTAOutput> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.GenomeSetToFASTA", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    /**
     * <p>Original spec-file function name: FeatureSetToFASTA</p>
     * <pre>
     * </pre>
     * @param   params   instance of type {@link us.kbase.kbasedataobjecttofileutils.FeatureSetToFASTAParams FeatureSetToFASTAParams} (original type "FeatureSetToFASTA_Params")
     * @return   instance of type {@link us.kbase.kbasedataobjecttofileutils.FeatureSetToFASTAOutput FeatureSetToFASTAOutput} (original type "FeatureSetToFASTA_Output")
     * @throws IOException if an IO exception occurs
     * @throws JsonClientException if a JSON RPC exception occurs
     */
    public FeatureSetToFASTAOutput featureSetToFASTA(FeatureSetToFASTAParams params, RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        args.add(params);
        TypeReference<List<FeatureSetToFASTAOutput>> retType = new TypeReference<List<FeatureSetToFASTAOutput>>() {};
        List<FeatureSetToFASTAOutput> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.FeatureSetToFASTA", args, retType, true, true, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }

    public Map<String, Object> status(RpcContext... jsonRpcContext) throws IOException, JsonClientException {
        List<Object> args = new ArrayList<Object>();
        TypeReference<List<Map<String, Object>>> retType = new TypeReference<List<Map<String, Object>>>() {};
        List<Map<String, Object>> res = caller.jsonrpcCall("KBaseDataObjectToFileUtils.status", args, retType, true, false, jsonRpcContext, this.serviceVersion);
        return res.get(0);
    }
}
