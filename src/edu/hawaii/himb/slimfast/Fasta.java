package edu.hawaii.himb.slimfast;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.List;
import java.util.ArrayList;


public class Fasta {
	private String inFileName;
	private int nbSequences;
	private BufferedReader buffReader;

	
    private List<String> headers;
    private List<String> sequences;
    private String nextline;
    private List<Sequence> seqs;
    private Utils utils;
    private int nbSubsets;
    public Fasta(String filename, int nbSubsets) {
    	
        try {
            this.buffReader = new BufferedReader(new InputStreamReader((InputStream) new FileInputStream(filename), Charset.defaultCharset()));
        }	
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        //this.headers = new ArrayList<String>();
        //this.sequences = new ArrayList<String>();
        this.nextline = "";
        this.seqs = new ArrayList<Sequence>();
        // use kmer size of 4
        this.utils = new Utils((byte) 4);
        this.nbSubsets = nbSubsets;
    }
   
    public double returnSignature(byte[] kmers, int [] band, int subsetNumber){
    	return utils.returnSignature(kmers, band, subsetNumber);
    	
    }
    
    public int [][] getRandomSamples(int nbSubsets, int subsetSize){
    	return utils.getRandomSamples(nbSubsets, subsetSize);
    }
    
       
    public void parse(int records) throws IOException {

    	Sequence seq = null; 

    	String line;
        StringBuilder seqStr = null;
    
        
        if (!this.nextline.isEmpty()) {
            this.nextline = "";
            seq = new Sequence(this.nbSubsets);
            seq.setId(Integer.parseInt(this.nextline.substring(1)));
            seqStr = new StringBuilder();
        }
        int processed = 0;
        while ((line = this.buffReader.readLine()) != null) {
            if ((line = line.trim()).isEmpty()) continue;
            if (line.startsWith(">")) {
                if (seqStr != null) {
                    seq.setSeqString(seqStr.toString());
                    seq.setKmersVec(this.utils.getKmers(seq.seqString));
                    
//                    System.out.println("-----");
                    seq.setKmersVec_k6(this.utils.getKmers_6k(seq.seqString));
                    seqs.add(seq);
                    ++processed;
                }
                if (processed == records) {
                    this.nextline = line;
                    return;
                }
                seqStr = new StringBuilder();
                seq = new Sequence(nbSubsets);
                seq.setId(Integer.parseInt(line.substring(1)));
                continue;
            }
            seqStr.append(line);
        }
        if (seqStr != null) {
            this.nextline = "";
            seq.setSeqString(seqStr.toString());
            seq.setKmersVec(this.utils.getKmers(seq.seqString));
            seq.setKmersVec_k6(this.utils.getKmers_6k(seq.seqString));
            seqs.add(seq);
        }
    }
    
    
    
    
    public int getNumSequences(){
    	return seqs.size();
    }
    	
    public List<Sequence> getSequences(){
    	return this.seqs;
    }
	
}
