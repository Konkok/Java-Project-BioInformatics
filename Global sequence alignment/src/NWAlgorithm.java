
public class NWAlgorithm {
	public static boolean verbose=false;
	 
	/*
	* prints out the score matrix 
	*/
	public static void dumpMatrix(int[][] matrix, 
	                              String row, 
	                              String column){
	         
	 System.out.print(String.format("%5s",""));
	 for (int j =0; j< row.length(); j++){
	   System.out.print(String.format("%5s", row.charAt(j)+"  "));
	 }          
	 System.out.println();
	          
	for (int i =0; i< column.length(); i++){
	  System.out.print(String.format("%5s",column.charAt(i)+ " "));
	  for (int j =0; j< row.length(); j++){
	     System.out.print(String.format("%5s", matrix[j][i]+" "));
	   }
	  System.out.println();
	 }
	}
	      
	/*
	* Needleman-Wunsch Dynamic Programming Algorithm
	* see http://amrita.vlab.co.in/?sub=3&brch=274&sim=1431&cnt=1
	* Runtime complexity: O(NxM), 
	* Space complexity: O(NxM) where N,M are lengths of the sequences
	*/
	public static AlignmentResult computeNWAlignment(String seq1, String seq2, SimpleAlignment parameters){
	             
	//It is easier to index and faster to use
	//integer matrix for fixed length storage
	int[][] scoreMatrix = new int[seq1.length()+1][seq2.length()+1];
	             
	int gapPenalty=parameters.getGapPenalty();
	int substitutePenalty=parameters.getSubstitutePenalty();
	             
	AlignmentResult result = new AlignmentResult();
	result.setParameters(parameters);
	             
	//Initialize the score matrix
	//the first row and column are for the gap
	//Complexity: O(NxM)
	for (int i =0; i< seq1.length()+1; i++){
	 for (int j =0; j< seq2.length()+1; j++){
	     scoreMatrix[i][j]=0;
	     if (i==0) 
	    scoreMatrix[i][j] = gapPenalty*j;
	     else if (j==0) 
	        scoreMatrix[i][j] = gapPenalty*i;
	  }
	}
	 
	if (verbose){
	   System.out.println("Initial Matrix");
	   dumpMatrix(scoreMatrix, "-"+seq1, "-"+seq2);
	}
	             
	int similarityCost=0;
	int matchCost=1;
	int seq1GapCost=0;
	int seq2GapCost=0;
	             
	//Compute the minimum cost scores between all 
	//possible pairs of prefixes
	//Complexity: O(NxM)
	for (int i =1; i< seq1.length()+1; i++){
	for (int j =1; j< seq2.length()+1; j++){
	                     
	//Case 1: The cost of mismatch between the two prefixes
	similarityCost= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : 
	                                   substitutePenalty;   
	 
	matchCost = scoreMatrix[i-1][j-1] + similarityCost;
	                     
	//Case 2: the cost of adding a gap on sequence 2
	seq2GapCost = scoreMatrix[i-1][j] + gapPenalty;
	                     
	//Case 3: the cost of adding a gap on sequence 1
	seq1GapCost = scoreMatrix[i][j-1] + gapPenalty;
	                     
	scoreMatrix[i][j] = Math.min(Math.min(matchCost,seq1GapCost),
	                             seq2GapCost);
	 }
	}
	             
	if (verbose){
	 System.out.println("\nFilled Matrix");
	 dumpMatrix(scoreMatrix, "-"+seq1, "-"+seq2);
	}
	             
	//Reconstruct the Alignment by backtracking on 
	//the score matrix
	//Complexity O(N+M)
	StringBuilder alignedSequence1= new StringBuilder();
	StringBuilder alignedSequence2= new StringBuilder();
	             
	int j = seq2.length();
	int i = seq1.length();
	             
	while (i >0 || j > 0) {
	 if (i>0 && j > 0) 
	    similarityCost= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : 
	                     substitutePenalty;
	    else
	    similarityCost = Integer.MAX_VALUE;
	                  
	if ( i > 0 && j >0 && 
	     scoreMatrix[i][j] == scoreMatrix[i-1][j-1] + similarityCost) { 
	    alignedSequence1.append(seq1.charAt(i-1));
	        alignedSequence2.append(seq2.charAt(j-1));
	        i=i-1;
	    j=j-1;
	 }
	 else if ( i > 0 && 
	      scoreMatrix[i][j] == scoreMatrix[i-1][j] + gapPenalty){
	    alignedSequence2.append("-");
	    alignedSequence1.append(seq1.charAt(i-1));
	    i=i-1;
	 }
	else if ( j > 0 && 
	       scoreMatrix[i][j] == scoreMatrix[i][j-1] + gapPenalty){
	    alignedSequence1.append("-");
	    alignedSequence2.append(seq2.charAt(j-1));
	    j=j-1;
	}
	} // end of while
	             
	 result.setTotalCost(scoreMatrix[seq1.length()][seq2.length()]);
	 result.setAlignmentLength(alignedSequence1.length());
	 result.setAlignments(new String[] 
	    {alignedSequence1.reverse().toString(), 
	     alignedSequence2.reverse().toString()});
	                 
	 return result;
	 }
	public static void testNWAlignment(String seq1, String seq2, int maxgap){
        
		NWAlgorithm.verbose=true;
		          
		System.out.println("\n\nTest Alignments using Needleman-Wunsch Algorithm");
		             
		     
		System.out.println("Original Sequences:");
		System.out.println(seq1);
		System.out.println(seq2);
		          
		long start=System.nanoTime();
		//SimpleAlignmentParameters par = new SimpleAlignmentParameters(-1,  -1);
		AlignmentResult result=NWAlgorithm.computeNWAlignment(seq1,
		                       seq2,
		                       new SimpleAlignment());
		long stop=System.nanoTime();
		          
		System.out.println( "elapsed time (sec): " + 
		                    String.format("%2.5f",
		                    (float)(stop-start)/Math.pow(10, 9)));
		String[] alignments= result.getAlignments();
		                  
		System.out.println("Alignments:");
		System.out.println(alignments[0] + " size:" + 
		                   alignments[0].length());
		           
		int matches=0;
		int gaps=0;
		int firstalign= 0;
		int lastalign=0;
		
		for (int k=0; k < alignments[0].length(); k++){
		if (alignments[0].charAt(k)==alignments[1].charAt(k)) {
			matches++;
			lastalign=k;
			if(matches==1){firstalign=k;}
		    System.out.print('|');
		} else System.out.print(" ");
		  
		
		  if ( (alignments[0].charAt(k)=='-') ||
		       (alignments[1].charAt(k)=='-')  ) 
		  gaps++;
		}
				
		System.out.println();
		           
		System.out.println(alignments[1] + " size:" + 
		                   alignments[1].length());
		System.out.println("Match Score=" + matches + ", " + 
		                   ((float)matches/alignments[0].length())+ 
		                   " gaps:" + gaps);
		           
		System.out.println("Edit Distance="+ result.getTotalCost());
		System.out.println("Alignment Length="+ 
		                    result.getAlignmentLength());
		System.out.println();
		int c1=0;
		int maxc=0;
		for(int i = firstalign; i<=lastalign;i++){
			if(alignments[0].charAt(i)!=alignments[1].charAt(i)){
				c1++;
				if (maxc<c1) maxc=c1;
			}else c1=0;
		}
		if(maxc > maxgap)
		{
			System.out.println("The sequences cannot be globaly aligned.");
		}
		
		int x = matches - gaps - (alignments[0].length() - matches - gaps); 
		System.out.println("Final score=" + x);
		System.out.println("first: " +firstalign+ "last: "+lastalign + "maxgap: "+maxc);
		}
}
