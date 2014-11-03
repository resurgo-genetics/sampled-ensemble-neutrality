/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *
 *									Class MuteCell
 *
 *  Each mute cell has a field - 'eigValue' - that represents a specific Second eigen value.
 *	The muteCell Object has a vector - eigResults, that holds all  EigResults that have
 *  this eigen Value.
 * 
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
import java.util.*;

class MuteCell{
	
	public String eigValue; // the eigResult.
	public Vector eigResults; // the vector that holds all the mutations...
	public boolean  containsWT; //true if this muteCell contains the wild type
	
	public MuteCell(String eigValue){
		this.eigValue = eigValue;
		eigResults = new Vector();
		containsWT = false;
	}
	
	
	public String toString(){
		return "|    "+ ((EigResult)eigResults.elementAt(0)).eigValue+"      |            "+((EigResult)eigResults.elementAt(0)).numOfVertices+
			       "             |         "+((EigResult)eigResults.elementAt(0)).WT+"        |    "+eigResults.size();
	}
}