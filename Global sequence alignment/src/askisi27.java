import java.util.Scanner;

public class askisi27 {

	public static void main(String[] args) {
		// TODO Auto-generated method stub

		Scanner reader = new Scanner(System.in);  // Reading from System.in
		System.out.println("Enter sequence1 : ");
		String seq1 = reader.nextLine(); // Scans the next token of the input as an string.
		
		Scanner reader2 = new Scanner(System.in);  // Reading from System.in
		System.out.println("Enter sequence2 : ");
		String seq2 = reader2.nextLine(); // Scans the next token of the input as an string.

		Scanner reader3 = new Scanner(System.in);  // Reading from System.in
		System.out.println("Enter a number: ");
		int maxgap = reader3.nextInt(); // Scans the next token of the input as an int.

		//System.out.println(seq1);
		//System.out.println(seq2);
		
		NWAlgorithm.testNWAlignment(seq1,seq2,maxgap);
	}

}
