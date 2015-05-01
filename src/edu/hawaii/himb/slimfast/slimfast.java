package edu.hawaii.himb.slimfast;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.IntStream;

import org.la4j.iterator.RowMajorMatrixIterator;
import org.la4j.matrix.sparse.CCSMatrix;
import org.la4j.matrix.sparse.CRSMatrix;

import com.google.common.collect.*;
import com.google.common.collect.Multiset.Entry;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Ints;



public class slimfast {

	//TODO, can we keep a representative sequence in each signature and compare some sequences to that signature
	// How do we choose the representative sequence?
	
	// TODO: make sure no sequences are left with -1. so far, removeSingletons can remove all the signatures of a sequence if it is alone in there.
	
	
	
	
	static boolean computeRandIndex = true;
	
	static final int MIN_NUM_COMPS = 40;
	static final int MIN_NB_MATCHES = 10;
	
	static final int MIN_NB_KMER_MATCHES = 220;
	static final int MIN_NB_COM_SIGS=8;

	
	//static String inFileName = "2M.fa";
	//static String inFileName = "1M.fa";

	//static String inFileName = "500k.fa";
	static String inFileName = "/Users/mahdi/Desktop/backup/80k.fa";
	//static String inFileName = "/tmp/2_seqs";
	
	
	//static String inFileName = "/Users/mahdi/Desktop/backup/2k.fa";
	//tatic String inFileName = "/Users/mahdi/Desktop/backup/430.fa";
	
	//static int nbSequences = 1975443;
	//static int nbSequences = 1000000;
	//static int nbSequences = 492322;
	static int nbSequences = 80001;
	//static int nbSequences = 7657;
	//static int nbSequences = 2;
	
	static int nbSubsets = 128;
	static byte subsetSize = 33;
	//static String correctClusteringFile;

	//static String correctClusteringFile = "/Users/mahdi/Desktop/backup/correctClust_vsearch_2k";
	//static String correctClusteringFile = "/tmp/correst_2.out";
	//static String correctClusteringFile = "/Users/mahdi/Desktop/backup/correctClust_vsearch_2k_93";
	static String correctClusteringFile = "/Users/mahdi/Desktop/backup/correctClust_vsearch_80k_93";
	//static String correctClusteringFile = "/Users/mahdi/Desktop/backup/out";
	
	


	
	
//	public static int nbKmerMatches(Sequence seq1, Sequence[] otherSeq){
//		int nbSeqMatches=0;
//		for (Sequence seq2: otherSeq){
//			if (seq1.seqId == seq2.seqId)
//				continue;
//			int nbKmerMatches = 0;
//			for (int i=0; i< seq1.kmersVec.length; i++){
//				if ((seq1.kmersVec[i] == seq2.kmersVec[i])){
//					// TODO make sure that the signature is not mostly 0's for those sequences in that cluster
//					nbKmerMatches++;
//				}
//			}
////			System.out.println(" Seq1:"+seq1.seqId+" Seq2:"+seq2.seqId+" have "+nbKmerMatches+" signature sizes in common");
////			System.out.println("Seq1"+Arrays.toString(seq1.kmersVec));					
////			System.out.println("Seq2"+Arrays.toString(seq2.kmersVec));
//			
//			if (nbKmerMatches >= MIN_NB_KMER_MATCHES){
//				nbSeqMatches++;
//			}
//		}
//		return nbSeqMatches;
//	}
	
	
	public static void droppingLargeClusters(List<Sequence> seqs, Map<Double, Signature> signatures, int[][] bands){
		List<Signature> largeSigs = new ArrayList<Signature>();
		for (Signature sig: signatures.values()){
			//TODO find a better way for this:  Why Dropping 500???
			if (sig.sequences.size() > 500){	
				//System.out.println("Signature "+sig.signature+" is large, it has "+sig.sequences.size()+" Sequences");
				Sequence seq= returnSequences(seqs, sig.sequences.subList(0, 1))[0];
				String vals="";
				int sum=0;
				for (int pos: bands[sig.bandId]){
					vals += seq.kmersVec[pos];
					if (seq.kmersVec[pos]!=0)
						sum++;
				}
				largeSigs.add(sig);
			}
		}
		for (Signature sig: largeSigs){
			removeSignature(sig, seqs, signatures);
		}
	}
	
	
	public static void removeSingleValSigs(Map<Double, Signature> signatures, Utils utils){
		int nbDropped=0;
		for (int prime : utils.primes_257){
			if (signatures.containsKey(prime)){
				signatures.remove(prime);
				nbDropped++;
			}	
		}
	}
	
	public static int numberOfFinalMutations(int[] positions){
		//System.out.println("myPOsitions"+Arrays.toString(positions));
		int start=0;
		int nbvalsHolding=0;
		int nbMutations=0;
		int pos;
		for(pos=1; pos<positions.length;pos++){
			if ((positions[pos]- positions[start]) <= 3){
				nbvalsHolding++;
			}else{
				if (nbvalsHolding > 0){
					//TODO test that that kmer pos-1 and kmer start overlap by at most one
					//TODO test that both kmers have overlapping values (0,1)(1,0)
					nbMutations++;
					nbvalsHolding=0;
				}
				start= pos;
			}
		}
		
		//System.out.println("Pos "+pos+" start "+start+" nbvalsHolding "+nbvalsHolding);
		if ((pos<positions.length) && ((positions[pos]-positions[start]) <= 3) && (nbvalsHolding>=1)){
			//(positions[start],positions[pos])
			//TODO test that that kmer pos-1 and kmer start overlap by at most one
			//TODO test that both kmers have overlapping values (0,1)(1,0)
			nbMutations++;
		} 
		return nbMutations;
	}		
	

		
		
		
	public static boolean nbKmerMatches2(int seq1Id, int seq2Id, List<Sequence> seqs){
		
		if (seq1Id == seq2Id)
			return true;
		Sequence mySeq1 = seqs.get(seq1Id);
		Sequence mySeq2 = seqs.get(seq2Id);

		int nbKmerMatches = 0;
		int nonMatches =0;
		int sumNonMatches=0;
		int sumNonMatches_6k=0;
			
		for (int i=0; i< mySeq1.kmersVec.length; i++){			
			if ((mySeq1.kmersVec[i] == mySeq2.kmersVec[i])){
				// TODO make sure that the signature is not mostly 0's for those sequences in that cluster
				nbKmerMatches++;
			}else{
				nonMatches++;
				sumNonMatches += Math.abs(mySeq1.kmersVec[i] - mySeq2.kmersVec[i]);
			}
		}
			
			

		if ((nbKmerMatches >= MIN_NB_KMER_MATCHES) && (nonMatches < 31) && (sumNonMatches < 31) ){
				//TODO:keep this in memory to speed things up?	
				//System.out.println(seq1.seqId+","+seq2.seqId+"\tseq1Kmers: "+seq1Kmers+" seq2Kmers: "+seq2Kmers+" nonZeroMatch:"+nonZeroMatch);
//				System.out.println("MATCH FOUND "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
//						+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
				return true;
		}else{
//			System.out.println("NO MATCHES "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
//					+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
			return false;
		}
	}
	

	
	
	public static int nbKmerMatches2(int seq1, int[] tempList,  List<Sequence> seqs){
		int nbSeqMatches=0;
		Sequence mySeq1 = seqs.get(seq1);		
		
		for (int seq2: tempList){
			if (seq1 == seq2)
				continue;
			
			Sequence mySeq2 = seqs.get(seq2);

			
			int nbKmerMatches = 0;
			int nonMatches =0;
			int sumNonMatches=0;
			int sumNonMatches_6k=0;
			List<Integer> mismatches= new ArrayList<Integer>();

		
			for (int i=0; i< mySeq1.kmersVec.length; i++){
				
				if ((mySeq1.kmersVec[i] == mySeq2.kmersVec[i])){
					// TODO make sure that the signature is not mostly 0's for those sequences in that cluster
					nbKmerMatches++;
				}else{
					nonMatches++;
					sumNonMatches += Math.abs(mySeq1.kmersVec[i] - mySeq2.kmersVec[i]);
					mismatches.add(i);
				}
			}
			
			
//			for (int i=0; i< seq1.kmersVec_k6.length; i++){
//				if ((seq1.kmersVec_k6[i] != seq2.kmersVec_k6[i])){
//					sumNonMatches_6k += Math.abs(seq1.kmersVec_k6[i] - seq2.kmersVec_k6[i]);
//					mismatches.add(i);
//				}
//			}
//			int minNbMutations = numberOfFinalMutations(Ints.toArray(mismatches));
			
//			System.out.println("Seq1"+Arrays.toString(seq1.kmersVec));					
//			System.out.println("Seq2"+Arrays.toString(seq2.kmersVec));
			//TODO modify this hardcoded requirement to  a dynamic one
			//if ((nbKmerMatches >= MIN_NB_KMER_MATCHES) && (nonMatches < 31) && (sumNonMatches < 31) && (sumNonMatches_6k < 40) && (minNbMutations<5)){
			if ((nbKmerMatches >= MIN_NB_KMER_MATCHES) && (nonMatches < 31) && (sumNonMatches < 31) ){
				//TODO:keep this in memory to speed things up?	
				//System.out.println(seq1.seqId+","+seq2.seqId+"\tseq1Kmers: "+seq1Kmers+" seq2Kmers: "+seq2Kmers+" nonZeroMatch:"+nonZeroMatch);

				nbSeqMatches++;
//				System.out.println("MATCH FOUND "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
//						+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
			}else{
//				System.out.println("NO MATCHES "+nbKmerMatches+"- Seq1:"+mySeq1.seqId+" Seq2:"+mySeq2.seqId
//						+" have "+nonMatches+" nonMatches "+sumNonMatches+" sumNonMatches"+ " and "+sumNonMatches_6k+" sumNonMatches_6k");
				
			}
		}
		return nbSeqMatches;
	}

	
	
	
	
	
	public static boolean areSimilar(int seq1Id, int seq2Id, int minSimilarSigs, List<Sequence> seqs){
		Sequence seq1 = seqs.get(seq1Id);
		Sequence seq2 = seqs.get(seq2Id);
		int nbSimilarSigs = 0;
		for (int i=0; i<seq1.signatures.length; i++){
			if ((seq1.signatures[i]!=null) && (seq1.signatures[i]==seq2.signatures[i])){
				nbSimilarSigs++;
			}
		}
		//TODO Move this in the For loop

		if (nbSimilarSigs >= minSimilarSigs){
			System.out.println("Seq1 "+seq1.seqId+" seq2 "+seq2.seqId+" have "+nbSimilarSigs+" shared signatures");
			return true;
		}

		//System.out.println("Seq1 "+seq1.seqId+" seq2 "+seq2.seqId+" are not similar and have only "+nbSimilarSigs+" shared signatures");
		return false;
	}
	

	public static void findAbherrentSigsUsingKmers(int[][] bands, List<Sequence> seqs, Map<Double, Signature> signatures, int[] correct_clusters){
		System.out.println("Total number of signatures to process for abherrations : "+signatures.size());

		int totalProcessed =0;

		
		Set<Integer> clusters = new HashSet<Integer>();		
		List<Integer> outliersList = new ArrayList<Integer>();
		List<Double> sigsToRemove=new ArrayList<Double>();

		for (Signature sig: signatures.values()){
			//System.out.println("----Working on sig:"+sig.signature+"\t"+Arrays.toString(sig.sequences.toArray()));
			
			clusters.clear();
			outliersList.clear();
	
			
			String line = "";
			
			// if sequence in that signature greater than MIN_NUM_COMPS
			// then randomly select a number of sequences to compare sequences against
			
			List<Integer> reps =new ArrayList(sig.sequences);
			Multiset<Integer>  counts = HashMultiset.create();
			
			
			if (sig.sequences.size() > MIN_NUM_COMPS){
				Collections.shuffle(reps);
				reps = reps.subList(0, MIN_NUM_COMPS);				
				for (int seq1id: sig.sequences){
					for(int rep: reps){
						if (areSimilar(seq1id, rep, MIN_NB_COM_SIGS, seqs) == true){
							counts.add(seq1id);
							counts.add(rep);
						}
					}
				}
				
			}else{
				int[] seqList = Ints.toArray(sig.sequences);
				for (int i=0; i<seqList.length-1; i++){
					for (int j=i+1; j<seqList.length; j++){
						if (areSimilar(seqList[i], seqList[j], MIN_NB_COM_SIGS, seqs) == true){
							counts.add(seqList[i]);
							counts.add(seqList[j]);
						}
					}					
					
				}
				reps =sig.sequences;
			}
			
			
			//Process the sequences in the the signature
			for (int seqId: sig.sequences){
				if (computeRandIndex==true){
					line +="\t"+seqId+":"+correct_clusters[seqId]+":"+counts.count(seqId)+"_";
					clusters.add(correct_clusters[seqId]);
				}
				//System.out.println("Sequence: "+seqId+" has "+counts.count(seqId)+" similar sequences in "+Arrays.toString(reps.toArray()));
				//TODO make the percentage a user def. value Instead of hard code it
				// TODO should update the counts but cannot do it here to now slow things down
				if(	counts.count(seqId) < 1){	
					//System.out.println("Sequence:"+seqId+" Does not have enough matches: REMOVED from sig: "+sig.signature);
					outliersList.add(seqId);
				}
			}

			if(outliersList.size()>0){
				removeSignatureFromSeqs(sig, Ints.toArray(outliersList), seqs, signatures);
			}
			//if we removed all sequences form the signature then add it to get rid of it
			if(sig.sequences.size()<=1){
				sigsToRemove.add(sig.signature);
				//System.out.println("TO REMOVE _SIG_ "+sig.signature+":"+sig.sequences.size()+" has duplicates "+clusters.size()+"\t"+line);

			}
			
			if (clusters.size() >1){
				System.out.println("myCounts for this signature"+counts.toString());
				if (outliersList.size()>=1){
					System.out.println("DETECTED:CORRECT "+Arrays.toString(outliersList.toArray())+"\t_SIG_:"+sig.signature+":"+sig.sequences.size()+" has duplicates "+clusters.size()+"\t"+line);
				}else{
					System.out.println("NOT_DETECTED _SIG_ "+sig.signature+":"+sig.sequences.size()+" has duplicates "+clusters.size()+"\t"+line);
				}				
			}else{
				if (outliersList.size()>1){
					//System.out.println("FALSE_POS"+Arrays.toString(outliersList.toArray())+"\t_SIG_:"+sig.signature+":"+sig.sequences.size()+" has duplicates"+line);
				}else{
					//System.out.println("NOTHING TO REPORT"+Arrays.toString(outliersList.toArray())+"\t_SIG_:"+sig.signature+":"+sig.sequences.size()+" has duplicates"+line);	
				}				
			}
			totalProcessed++;
			if (totalProcessed % 100000 ==0){
				System.out.println("Number of Signatures Processed = "+totalProcessed);
			}
		}
		
		for (Double sig: sigsToRemove){
			signatures.remove(sig);
			
		}
	}

	
	
	
	
	
	public static Sequence[] returnSequences(List<Sequence> seqs, List<Integer> tempList){
		Sequence[] otherList = new Sequence[tempList.size()];
		int nextPos = 0;
		for (int seqId: tempList){
			otherList[nextPos] = seqs.get(seqId);	
			nextPos++;
		}
		return otherList;
	}
	
	
	public static int getNumberValidKmers(Sequence seq1, Sequence seq2, int[] signatures, int[][] bands){
		int nbValidKmers=0;
		Set<Integer> kmers = new HashSet<Integer>();
		for (int sigId: signatures){
			kmers.addAll(Ints.asList(bands[sigId]));
		}
		
		for (int kmer: kmers){
			if ((seq1.kmersVec[kmer]==seq2.kmersVec[kmer]) && (seq2.kmersVec[kmer]!=0)){
				nbValidKmers++;
			}
		}
		return nbValidKmers;
	}
	
	
	/***********
	 * for a given Seq1, find the best match among otherSeqs then computes the number of 
	 * non-null shared sequences between seq1 and its best hit in otherSeqs
	 * @param seq1
	 * @param otherSeq
	 * @return
	 */
	public static int nbKmerMatches(Sequence seq1, Sequence[] otherSeqs, int[][] bands){
		int nbKmerMatches=0;
		List<Integer> colMatches = new ArrayList<Integer>();		
		Sequence bestHitSequence=null;
		List<Integer> bestColMatches = new ArrayList<Integer>();
		
		for (Sequence seq2: otherSeqs){
			if (seq1.seqId == seq2.seqId)
				continue;
			colMatches.clear();
			for (int i=0; i< nbSubsets; i++){
				if ((seq1.sigSizes[i] == seq2.sigSizes[i]) && (seq1.sigSizes[i] != 1)){
					colMatches.add(i);
				}
			}

			if ((colMatches.size() > MIN_NB_MATCHES) && (colMatches.size() > bestColMatches.size())){
				bestHitSequence = seq2;
				bestColMatches = new ArrayList<Integer>(colMatches);
			}
		}
		// Get number non empy kmers sor seq1 and seq2 on cols in bestColMatches
		if (bestHitSequence != null){
			System.out.println("seq1: "+seq1.seqId+"\tbestHitId: "+bestHitSequence.seqId+"\t"+Arrays.toString(bestColMatches.toArray())+
					"\tscore is:"+getNumberValidKmers(seq1, bestHitSequence, Ints.toArray(bestColMatches), bands));
		}else{
			System.out.println("seq1: "+seq1.seqId+"Has no best hit");
		}
		return getNumberValidKmers(seq1, bestHitSequence, Ints.toArray(bestColMatches), bands);
	}	
	



	
	public static int getPercentile(Sequence seq, int percentile){
		List<Integer> sigSizes = new ArrayList<Integer>();
		for (Signature sig : seq.signatures){
			if (sig == null)
				sigSizes.add(0);
			else
				sigSizes.add(sig.sequences.size());

		}
		Collections.sort(sigSizes);		
		return sigSizes.get(percentile);
	}
	
	
	
	public static int getBiggestSigSize(Sequence seq, double notThisSig){
	/**
	 * This method takes as input a sequence, and finds the the size 
	 * of its largest signature(number of sequences in that sig.), excluding
	 * the signature passed as a param (notThisSig)
	 */
		int maxSize = 0;
		for (Signature sig: seq.signatures){
			if (sig != null && sig.signature != notThisSig){
				if (sig.sequences.size() > maxSize){
					maxSize = sig.sequences.size(); 
				}
			}
		}
		return maxSize;
	}
	
    public static int[] getCorrectClustering(String fileName, int nbSequences) throws IOException{
    	int[] clusters = new int[nbSequences];
    	BufferedReader br = new BufferedReader(new FileReader(fileName));
        try {
            String line = "";
            int clust = -1;
            line = br.readLine();
            while (line != null) {
                if ( (line.length() > 7)  && line.substring(0, 7).equals("Cluster")){
                	clust = Integer.parseInt(line.split(" ")[1]);
                	
                }else{
                	clusters[Integer.parseInt(line)] = clust;
                }
                line = br.readLine();
            }
        } finally {
            br.close();
        }
        return clusters;
    }

    public static float computeModRandIndex(int [] correct_clustering, int [] new_clustering){
        int T_T = 0; // Together in A and together in B
        int	S_S = 0; // Separate in A and Seaparate in B 
        int S_T = 0; // SeparaTte in A and Seaparate in B 
        int T_S = 0; // Together in A and Separate in B
		
//        System.out.println("*******");
//        System.out.println(Arrays.toString(correct_clustering));
//        System.out.println(Arrays.toString(new_clustering));
//        System.out.println("*******");
        
        int[] pair;
		CombinationGenerator x = new CombinationGenerator (nbSequences, 2);
		while (x.hasMore()) {
		  pair = x.getNext ();
		  
		  if (correct_clustering[pair[0]] ==  correct_clustering[pair[1]]){ 
			  if (new_clustering[pair[0]] == new_clustering[pair[1]]){
				  T_T +=1;
			  }else{
				  System.out.println(pair[0]+" and "+pair[1]+" were together in clust: "+ 
						  correct_clustering[pair[0]]+" now are separated in clusts: "+ 
						  new_clustering[pair[0]]+" and "+new_clustering[pair[1]]);
				  T_S += 1;
			  }
		  }else if (correct_clustering[pair[0]] !=  correct_clustering[pair[1]]){
			  if (new_clustering[pair[0]] != new_clustering[pair[1]]){
				  S_S += 1;
			  }else{
				  S_T += 1;
			  }
		  }
		}
		return  ( (float) T_T + S_S +  S_T ) / (T_T + S_S +  S_T + T_S);
    }
    
    
	
	public static void removeSingletons(List<Sequence> seqs, Map<Double, Signature> signatures){
		/* Sets all the signatures of size one of a sequence to null.
		 * */
		List<Signature> suspicousSigs = new ArrayList<Signature>();
		Iterator<Signature> i = signatures.values().iterator();
		for(Signature sig: signatures.values()){
			if ((sig.sequences.size() == 1) || (sig.sequences.size() > 500) ){
				suspicousSigs.add(sig);
			}
		}
		for (Signature sig: suspicousSigs){
			removeSignature(sig, seqs, signatures);
		}
	}
	
	

	
	
	
	
	public static int[] findClustering(List<Sequence> seqs, Map<Double, Signature> signatures){
		int[] clusters = new int[nbSequences];
		LinkedList<Integer> queue = new LinkedList<Integer>();
		Set<Integer> neighborsSet = new HashSet<Integer>();
		// Set all values by default to -1
		for (int i=0; i< clusters.length; i++){
			clusters[i] = -1;
		}
		int nextClustId = -1;
		for (Sequence seq: seqs){	
			if (clusters[seq.seqId]==-1){
				// Assign it to its own cluster
				// Cannot belong to a previous cluster, or else it would have been labeled through another sequence
				nextClustId++;
				clusters[seq.seqId]= nextClustId;
			}else{
				continue;
			}
			// Add the sequence to a queue so as to inspect its children
			queue.push(seq.seqId);
			//System.out.println("Queued : "+ seq.seqId);

			while(queue.isEmpty()==false){
				int newSeqId = queue.removeLast();
				//System.out.println("unqueuing:"+ newSeqId);
				// For each signature that this sequence is involved in
				// get all the other sequences that are in the signature
				for (Signature sig: seqs.get(newSeqId).signatures){
					// once a signature is visited, or if null (had only one read in it)
					// we do not need to process it further
					if (sig == null || sig.visited==true)
							continue;

					// TODO WAS THIS WHAT WAS CAUSING THE ERRORS
					// search for the position where the sequence 
					//int myPosition = sig.getPosition(newSeqId);
					int myPosition = 0;
					int seqSig;
					// We start at the position of the sequence since all sequences
					// before it would have been assigned already
					// since seqIds are incremental
					// we are sure any sequence before it has already a clusterId 
					for (int i = myPosition; i<sig.sequences.size(); i++ ){
						seqSig = sig.sequences.get(i);
						if(clusters[seqSig] == -1){
							// is sequence was not assigned to a cluster
							// then add it to the queue to process other sequences in its cluster
							neighborsSet.add(seqSig);
						}
					}
					sig.visited=true;
				}
				// add the neighbors (sequences sharing this signature) to the queue 
				// so that we can add their neighbors
				for(int seqId: neighborsSet){
					clusters[seqId] = clusters[seq.seqId];
					queue.addFirst(seqId);
					//System.out.println("Queued : "+seqId+" through: "+ seq.seqId);
				}
				neighborsSet.clear();
				//if ((queue.size() % 10000)==0)
				//System.out.println(" queue is: "+queue.size()+" in size");
			}
			
			if ((seq.seqId % 10000)==0)
				System.out.println(" X ");
		}
		System.out.println("");
		
		// System.out.println("seq: "+i+" is in cluster: "+clusters[i]);
		System.out.println("Number of clusters generated is: "+(nextClustId+1));
		return clusters;
	}

	
	public static void removeSignature(Signature sig, List<Sequence> seqs, Map<Double, Signature> signatures){
		for (int seqId: sig.sequences){
			//REMOVE THIS CHECK AFTER VALIDATING THE seqId = seqs.get(seqId).seqId
			if (seqId != seqs.get(seqId).seqId){
				System.out.println("ERROR IN removeSiganture: seqId != seqs.get(seqId).seqId ");
				System.exit(0);
			}
			seqs.get(seqId).signatures[sig.bandId]=null;
		}
		signatures.remove(sig.signature);
	}
	
	public static void removeSignatureFromSeqs(Signature sig, int[] seqIds, List<Sequence> seqs, Map<Double, Signature> signatures){
		for (int seqId: seqIds){
			//REMOVE THIS CHECK AFTER VALIDATING THE seqId = seqs.get(seqId).seqId
			if (seqId != seqs.get(seqId).seqId){
				System.out.println("ERROR IN removeSiganture: seqId != seqs.get(seqId).seqId ");
				System.exit(0);
			}
			int position = sig.getPosition(seqId);
			sig.sequences.remove(position);
			seqs.get(seqId).signatures[sig.bandId]=null;
		}
		//for thread safete not remove this here. Do another look 
		//		if (sig.sequences.size() ==0)
		//			signatures.remove(sig.signature);
	}
	
	
	
//	public static void removeSignature(Sequence[] outliers, Signature sig){
//		//System.out.println("I AAAAAAAAAAAAMMMMMM HERE"+);
//		// if removing two sequence from signature of size 1, stop after removing the first one.
//		
//		for (Sequence seq: outliers){			
//			for (int i=0; i<seq.signatures.length; i++){
//				if ((seq.signatures[i]!=null) && (seq.signatures[i].signature == sig.signature)){
//					seq.signatures[i]=null;
//				}
//			}
//		}
//	}
	
	
	
	public static void main(String[] args) throws IOException {
		
		Utils utils = new Utils((byte) 4);
		
		//System.exit(0);

		
		
		
		
		int[] correct_clusters =null;
		if (computeRandIndex == true)
			correct_clusters =	getCorrectClustering(correctClusteringFile, nbSequences);
		
		
		long startTime = System.nanoTime();
		
		System.out.println("Starting");
		Fasta myFasta = new Fasta(inFileName, nbSubsets);
		
		myFasta.parse(nbSequences);
		System.out.println("Nb Sequences is: " + myFasta.getNumSequences());
		List<Sequence> seqs = (List<Sequence>) myFasta.getSequences();		
	
		
//		// TODO Change this (tempFix )
//		Sequence [] seqs = new Sequence[seqsList.size()];
//		for (int i=0; i< seqsList.size(); i++){
//			seqs[i] =seqsList.get(i); ;
//		}
//		System.out.println("Done convertine sequences list ot an array");
		
		Map<Double, Signature> signatures = new HashMap<Double, Signature>(nbSequences * nbSubsets); // allocate max number of possible signatures
		int [][] bands = myFasta.getRandomSamples(nbSubsets, subsetSize);

		
		
		
		
		System.out.println("Number of bands is: "+bands.length);
		
		int nbSequencesProcessed = 0;
		long s = System.nanoTime();
		int subsetNumber = 0; // this is the band number
		int bandNumber=0;
		
		for (int [] band: bands){
			System.out.println(subsetNumber+"\t"+Arrays.toString(band));
			for (Sequence seq: seqs){
				// if signature does not exist, create it and add it to seq
				// need to check because of sig collisions
				Double sig = myFasta.returnSignature(seq.kmersVec, band, subsetNumber);
				Signature signature = signatures.get(sig);
				
			
				if (signature == null){
					signature = new Signature(sig);
					signatures.put(sig, signature);
				}
				signature.bandId=bandNumber;
				signature.sequences.add(seq.seqId);
				seq.addSignature(signature, subsetNumber);
				nbSequencesProcessed++;
			}
			subsetNumber++;	
			if (( subsetNumber % 10) == 0){ 
				System.out.print(subsetNumber+"\t nbSigs is: "+signatures.size()+"\t");
				System.out.println((System.nanoTime() - s)/1000000000);
				System.out.println("started cleaning");
				removeSingletons(seqs, signatures);
				System.out.println("Done cleaning: number of remaining signatures is: "+signatures.size());
				s = System.nanoTime();
			}
		bandNumber++;
		}
		System.out.println("--------------------------------------------------------------");
	
		

		
		
		//COMPUTE PERCENTYLE SEQUENCE SIZES
		System.out.println("Computing Percentile Sizes for all sequences "+(System.nanoTime() - startTime)/1000000000);
		for (Sequence seq: seqs){
			// TODO FIX this to get percentile based on number of signature. Here we assume it is 100
			seq.pSize = getPercentile(seq, 94);
		}
		System.out.println("Done computing percentyle sizes "+(System.nanoTime() - startTime)/1000000000);
		
		
		//COMPUTE PERCENTYLE SEQUENCE SIZES
		System.out.println("Compute the number of unique kmers in each sequence "+(System.nanoTime() - startTime)/1000000000);
		for (Sequence seq: seqs){
				seq.nbKmers = Sets.newHashSet(Bytes.asList(seq.kmersVec)).size() -1;
		}
		System.out.println("Done computing nbUnique Kmers for each sequence "+(System.nanoTime() - startTime)/1000000000);

		
		
		
		
		
//		System.out.print("SEQUENCES PERCENTILE AND SIGANTURES SIZES ARE : +++++++++++++++++++++++ ");
//		for (Sequence seq: seqs){
//			System.out.print("SEQ: "+seq.seqId+":"+correct_clusters[seq.seqId]+":"+seq.pSize+"\t");
//			for (Signature sig: seq.signatures){
//				if (sig == null)
//					System.out.print(1+" ");
//				else
//					System.out.print(sig.sequences.size()+" ");
//			}
//			System.out.print("\n");
//		}
			
		
			
		

		
		System.out.println("updating signature sizes: "+(System.nanoTime() - startTime)/1000000000);
		for (Sequence seq: seqs){
			for (int i=0; i< nbSubsets; i++){
				if (seq.signatures[i] != null)
					seq.sigSizes[i]= seq.signatures[i].sequences.size();
			}
		}
		System.out.println("Done computing signature sizes: "+(System.nanoTime() - startTime)/1000000000);

		System.out.println("Dropping large signatures: "+(System.nanoTime() - startTime)/1000000000);
		droppingLargeClusters(seqs, signatures, bands); //cluster is a signature!!!!!
		System.out.println("Done Dropping large signatures: "+(System.nanoTime() - startTime)/1000000000);

		System.out.println("Dropping signatures with prime values: "+(System.nanoTime() - startTime)/1000000000);
		removeSingleValSigs(signatures, utils); // if signatue equals the prime number, then remove
		System.out.println("Done Dropping signatures with prime values: "+(System.nanoTime() - startTime)/1000000000);
	
		//findAbherrentSequences(seqs, signatures, correct_clusters);
		findAbherrentSigsUsingKmers(bands, seqs, signatures, correct_clusters);
		System.out.println("Done identifying outliers in all signatures"+(System.nanoTime() - startTime)/1000000000);

		
		
		

		System.out.println("Cleaning graphs");		
		removeSingletons(seqs, signatures);
		System.out.println("Done cleaning graphs");
		System.out.println("Finding Clusters");

		int[] new_clusters = findClustering(seqs, signatures);
		//System.out.println(Arrays.toString(new_clusters));
		System.out.println("Done Finding Clusters\n\n");
		
		if (computeRandIndex == true){			
			System.out.println("Computing Rand index\n");
			System.out.println("Mod rand Index is: "+computeModRandIndex(correct_clusters, new_clusters));
			System.out.println("Done");
		}
		
		long endTime = System.nanoTime();
		long duration = (endTime - startTime)/1000000000; 
		
		System.out.println("Completed Successfully "+duration);
		
		
		
//		Sequence seq = seqs.get(950);
//		Multiset<Integer> seqsCounts = HashMultiset.create();
//		for (Signature sig: seq.signatures){
//				if (sig != null)
//					seqsCounts.addAll(sig.sequences);
//		}
//		System.out.println("Done:---------------");
//
//		for (Entry<Integer> a: seqsCounts.entrySet()){
//			System.out.println("voila:  "+a);			
//		}
	}
}
