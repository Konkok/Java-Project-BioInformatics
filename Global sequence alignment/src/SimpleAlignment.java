
public class SimpleAlignment {
	private int gapPenalty;
    
    private int substitutePenalty;
 
    public int getGapPenalty() {
        return gapPenalty;
    }
 
    public void setGapPenalty(int gapPenalty) {
        this.gapPenalty = gapPenalty;
    }
 
    public int getSubstitutePenalty() {
        return substitutePenalty;
    }
 
    public void setSubstitutePenalty(int substitutePenalty) {
        this.substitutePenalty = substitutePenalty;
    }
 
    public SimpleAlignment(int gapPenalty,int substitutePenalty) {
        super();
        this.gapPenalty = gapPenalty;
        this.substitutePenalty = substitutePenalty;
    }
     
    public SimpleAlignment() {
        super();
        this.gapPenalty = 1;
        this.substitutePenalty = 1;
    }

}
