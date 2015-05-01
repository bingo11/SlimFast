package edu.hawaii.himb.slimfast;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Chars;


public class Utils {

	
	public final static int[] primes_257 = new int[]{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621};
	private double[] log_primes = new double[this.primes_257.length];

	
	private static final Map<Character, Integer> basesMap = new HashMap<Character, Integer>(5);
	static {
        basesMap.put(Character.valueOf('A'), 0);
        basesMap.put(Character.valueOf('C'), 1);
        basesMap.put(Character.valueOf('G'), 2);
        basesMap.put(Character.valueOf('T'), 3);
        basesMap.put(Character.valueOf('N'), 0);
    }
	byte kmerSize;
	
	public Utils(byte kmerSize){
	    for (int i = 0; i < this.primes_257.length; ++i) {
	        this.log_primes[i] = Math.log10(this.primes_257[i]);
	    }
        this.kmerSize = kmerSize;
	}
	
	public double returnSignature(byte[] kmers, int [] band, int subsetNumber){

		double signature=0;
		for(int col: band){
			signature +=   (kmers[col] * log_primes[col]);
		}
		signature +=   (subsetNumber * log_primes[256]);
		return signature;
	}
	
	public byte[] getKmers(String seq){
		char[] sequence = seq.toCharArray();		
		String base = seq.substring(0, this.kmerSize);
		int baseEncoding = getKmerEncoding(base);
		int n = baseEncoding;
		int nbKmers = (int)Math.pow(4.0, this.kmerSize);
		byte[] counts = new byte[nbKmers];
		counts[n] = (byte)(counts[n] + 1);
        for (int j = this.kmerSize; j < seq.length(); ++j) {
            baseEncoding = getHash(baseEncoding, sequence[j], this.kmerSize);
            counts[baseEncoding] = (byte)(counts[baseEncoding] + 1);
        }	
		return counts;
	}
	public byte[] getKmers_6k(String seq){
		//hardcoded to 6
		char[] sequence = seq.toCharArray();		
		String base = seq.substring(0, 6);
		int baseEncoding = getKmerEncoding(base);
		int n = baseEncoding;
		int nbKmers = (int)Math.pow(4.0, 6);
		byte[] counts = new byte[nbKmers];
		counts[n] = (byte)(counts[n] + 1);
        for (int j = this.kmerSize; j < seq.length(); ++j) {
            baseEncoding = getHash(baseEncoding, sequence[j], 6);
            counts[baseEncoding] = (byte)(counts[baseEncoding] + 1);
        }	
		return counts;
	}	
	
	
	public int[] getBases(String seq){
		int[] basesCounts = new int[4];
		char[] basesList = seq.toCharArray();		
		Multiset<Character> bases = HashMultiset.create(Chars.asList(basesList));
		basesCounts[0] = bases.count('A');
		basesCounts[1] = bases.count('C');
		basesCounts[2] = bases.count('G');
		basesCounts[3] = bases.count('T');
		return basesCounts;
	}
	
	
	
	public int [][] getRandomSamples(int nbSubsets, int subsetSize){
		Integer[] cols = new Integer[(int)Math.pow(4.0, this.kmerSize)];
		for (int i=0; i< cols.length; i++ ){
			cols[i] =i;
		}
		int [][] bands = new int[nbSubsets][subsetSize];
		
		for(int i=0; i< nbSubsets; i++){
			List<Integer> colList = Arrays.asList(cols);
			Collections.shuffle(colList);
			for(int j=0; j< subsetSize; j++){
				bands[i][j] = cols[j];
			}
		}
		return bands;
	}
	
    private static int getKmerEncoding(String kmer) {

        int hashCode = 0;
        for (int i = 0; i < kmer.length(); ++i) {
            char base = kmer.charAt(i);
            hashCode<<=2;
            hashCode|=basesMap.get(Character.valueOf(base)).intValue();
        }
        return hashCode;
    }
    
    private int getHash(int hashCode, char currentChar, int kmerSize) {

        hashCode&=~ (1 << 2 * kmerSize - 1);
        hashCode&=~ (1 << 2 * kmerSize - 2);
        switch (currentChar) {
            case 'A': {
                hashCode<<=2;
                hashCode|=0;
                break;
            }
            case 'C': {
                hashCode<<=2;
                hashCode|=1;
                break;
            }
            case 'G': {
                hashCode<<=2;
                hashCode|=2;
                break;
            }
            case 'T': {
                hashCode<<=2;
                hashCode|=3;
                break;
            }
            case 'N': {
                hashCode<<=2;
		hashCode|=0;
                break;
            }
            default: {
                System.out.println("Unknown char");
            }
        }
        return hashCode;
    }
	
}
