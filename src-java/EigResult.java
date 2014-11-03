/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 * 							Class EigResult
 *
 *  each EigResult represents a mutation that was done to the wild typ sequence.
 *
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

import java.util.*;

class EigResult{
	
	public String name; //the mutation like 15-U , means: U is the new point mutation in place 15.
	public String full_name; //the mutation like C15U , means: C was changed to U in place 15.
	public String eigValue; //the Second eigenValue sequence
	public String numOfVertices;
	public String repr;//the representation : Shapiro style. 
	public String WT; //the wild-type sequence
	public int numOfAppiar; //how many  mutations with the same EigenValue exist 
	public double distance; //the distance from the wild type sequence - refers to shapiro representation
	public double dbDistance;//the distance from the wild type sequence - refers to "dot-bracket" representation
	public String deltaG; // minimum free Energy in Kcal per mol
	public String dotBracket;//the representation : dot-bracket style. 


	
	public EigResult(){
		name = "";
		full_name="WT";
		eigValue = "";
		numOfVertices = "";
		numOfAppiar = 0;
		repr = "";
		WT = " - ";
		distance = 0 ;
		dbDistance = 0;
		deltaG = " 0 ";
		dotBracket = "";
		

	}
	
	public String toString(){
		return "|    "+eigValue+"      |            "+numOfVertices+
		       "             |         "+WT+"        |    "+numOfAppiar;
	}
	
	public void setShapiroDistance (double toSet){
		distance = toSet;
	}
	
	public void setDbDistance (double toSet){
		dbDistance = toSet;
	}
	
	public void set_full_name(String new_name){
		full_name = new_name;
	}
	
	public void set_new_Energy(){
		String temp = deltaG.substring(1,deltaG.length()-1);
		deltaG=temp;
	}
}
