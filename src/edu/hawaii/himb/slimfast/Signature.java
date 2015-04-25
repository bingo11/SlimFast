package edu.hawaii.himb.slimfast;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;

public class Signature{
	public boolean visited= false;
	public double signature;
	public List<Integer> sequences = new ArrayList<Integer>();
	
	public int bandId; // the band number that was used to generate this signature
	
	Comparator<Integer> c = new Comparator<Integer>() {
	      public int compare(Integer s1, Integer s2) {
	    	return (s1).compareTo(s2);
	      }
	};
	
	
	public Signature(double signature){
		this.signature = signature;
	}
	
	public void addSequenceSign(int seq){
		sequences.add(seq);
	}
	
	/**
	 * Finds the poisiton of the sequence id seq in the list of sequences 
	 * belonging to the current signature
	 * @param  url  an absolute URL giving the base location of the image
	 * @param  name the location of the image, relative to the url argument
	 * @return      the image at the specified URL
	 * @see         Image
	 */	public int getPosition(int seq){
		return Collections.binarySearch(this.sequences, seq, c);
	}
	

}
