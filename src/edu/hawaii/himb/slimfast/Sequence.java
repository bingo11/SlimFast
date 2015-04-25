package edu.hawaii.himb.slimfast;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;




public class Sequence {

	
	int seqId;
	String seqString;
	byte[] kmersVec;
	byte[] kmersVec_k6;
	
	Signature[] signatures;
	int[] sigSizes;
	
	Multiset<Integer> seqsCounts = HashMultiset.create();

	int pSize=0; // size of the percentile signature that this seq. belongs to
	int nbKmers=0; // number of non null kmers in the sequence	
	
	int nbSubsets;
	
	public Sequence(int nbSubsets){
		signatures = new Signature[nbSubsets];
		sigSizes = new int[nbSubsets];
	}
	

	public void setId(int seqId) {
		this.seqId= seqId;
	}
	
	public void setSeqString(String seqString) {
		this.seqString = seqString;
	}
	
	public void setKmersVec(byte[] kmersVec){
		this.kmersVec = kmersVec;
	}

	public void setKmersVec_k6(byte[] kmersVec){
		this.kmersVec_k6 = kmersVec;
	}
	
	public void addSignature(Signature signature, int subsetNumber){
		signatures[subsetNumber] = signature;
	}

	

	
}
