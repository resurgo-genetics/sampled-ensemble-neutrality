/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*                         
*                        class ProcessResult
*
*  this class is responsible for processing the "result" file, that is produced 
*  after running the  mute_single proceduer.
*  the class containd s buffers that reads the file and data structurs that 
*  accumulates the necessey data.
*  this class is also, produce the HTML files that are used in the interface, 
*  and prints to them the required information.
*
*
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


import java.util.*;
import java.io.*;
import java.lang.*;

class ProcessResult{
	
//	public Vector mutations;;
        	
	public ProcessResult(){
//		mutations = new Vector();
            Double dotBracketdist;
            Double shapiroDist;
	}
    
    public static void go () {
        go("RESULT_TABLE.html");
    }

	/* the main function that runs it all...*/ 
	public static void go(String outfile){
		 	// cleaning all the  leftovers from the last executaion	
	    try{
			Process p = Runtime.getRuntime().exec(new String[] {  "/bin/sh", "-c", "rm -r -f ./htmlDir/ "});	
	        p.waitFor();
		}
		catch (Exception e){System.out.println(e);}
		
		/***** needed fields that holds the data *******/
		
		//creats a directory that will hold all the html files.
    	String adress="htmlDir/"; 	
    	File htmlDir = new File(adress);
    	htmlDir.mkdir();

	    

        BufferedReader reader = null;
        BufferedReader seqReader=null;//reads the wild type file
        BufferedReader jumpsReader = null;
        String our_seq="";//holds the wild type sequence
        String inputData="", inputLine="";
        String the_wild_type_representation="";//shapiro repr.
        String the_wild_type_repr_dot_b="";//dot-bracket repr.
        String the_wild_type_eigenValue="";
        String the_wild_type_freeEnergy="";
        Vector finalResults = new Vector(); //the vector that will hold the mutations
        Vector finalResults_backup = new Vector();//contains the mutations after we delete them from the finalResults vector
		Vector allEigResults = new Vector();
	    Vector allEigResults2 = new Vector();
	    
	    int shapiroJump = 6;
	    int dbJump = 4;
        int mut_counter=0;//counts the number of the mutations
        double sum = 0; //sums the Shapiro distances of the mutations
        double db_sum = 0; //sums the Shapiro distances of the mutations
        double average_distance = 0; //Shapiro average distance
        double db_average_distance = 0; //dot barket average distance
        String string_average_distance=""; //the string representation of the distance
        String db_string_average_distance="";


        try{
            reader = new BufferedReader(new FileReader("result"));
        }
        catch(FileNotFoundException e){
            System.err.println("Error: file result"+/*args[0]+*/" not found!");
            System.exit(1);
        }
        
        try{
            seqReader = new BufferedReader(new FileReader("our_sequence"));
        }
        catch(FileNotFoundException e){
            System.err.println("Error: file our_sequence"+/*args[0]+*/" not found!");
            System.exit(1);
        }
        
        try{//to take the values of DB and Shpiro Clustering Resolution  for the 3rd &4th tables printers
            jumpsReader = new BufferedReader(new FileReader("Clustering Resolution.txt")); 
        }
        catch(FileNotFoundException e){
            System.err.println("Error: file Clustering Resolution "+/*args[0]+*/" not found!");
            System.exit(1);
        }
        try{ 
            String jump_str = jumpsReader.readLine();
            int pass = jump_str.indexOf(" ");
            String sJump = jump_str.substring(0,pass);
            String dJump = jump_str.substring(pass+1);
            shapiroJump = stringToInt(sJump,6);
            dbJump = stringToInt(dJump,4);
            
        }
        catch(IOException e){
            System.err.println("Error: error while reading input file!");
            System.exit(1);
        } 
        
        /* starting to read the first line (the wild type).
           saving the sequence and manipulating it:
           1.replacing all its characters to upper case.
           2. replacing Ts with Us. */
        try{ 
            String temp = seqReader.readLine();
            String temp2 =temp.toUpperCase();
            our_seq= temp2.replace('T','U');
                  
        }
        catch(IOException e){
            System.err.println("Error: error while reading input file!");
            System.exit(1);
        } 
        // continue reading the "result" file
        try{
            while ((inputLine = reader.readLine())!=null)
                inputData+= inputLine+"\n";
        }
        catch(IOException e){
            System.err.println("Error: error while reading input file!");
            System.exit(1);
        }
        
        StringTokenizer tokens = new StringTokenizer(inputData);
        String numOfVerticesOfWT = "0";

        // extracting all the data about the "WILD TYPE"
        // creating a "EigResult"(a mutation)
        // and initialise all its fields
        try{	
      		EigResult result = new EigResult();//new EigResult item
        	if (tokens.hasMoreTokens()){
        		result.dotBracket= tokens.nextToken();
        		result.deltaG = tokens.nextToken();
        		result.set_new_Energy();
        		
        	    result.eigValue = tokens.nextToken();
        	    result.numOfVertices = tokens.nextToken();
        	 	result.repr = tokens.nextToken();
        	 	result.numOfAppiar++;
        	 	result.WT = "WT";
        	 	numOfVerticesOfWT = result.numOfVertices;
        	 	the_wild_type_eigenValue= result.eigValue;
        	 	the_wild_type_freeEnergy = result.deltaG;
        	 	the_wild_type_representation = result.repr;
        	 	the_wild_type_repr_dot_b = result.dotBracket;	
        	}

        	//creates a new muteCell
        	MuteCell mc = new MuteCell(result.eigValue); 
        	mc.containsWT = true; //the w.t is in this MuteCell
        	mc.eigResults.addElement(result); /// adding this mutation to the MuteCell
        	finalResults.addElement(mc); //adding the MuteCell to the final mutations vector
            allEigResults.addElement(result);
			allEigResults2.addElement(result);
			


        	//continue with all the other mutations...
        	while(tokens.hasMoreTokens()){
        		result = new EigResult();
        		boolean flag = true; //indicates if we need to open a new MuteCell
        							// or if it is the same EigenValue, so we do'nt 
        		result.name = tokens.nextToken();
        	    // finds the full name of the mutaion like: U15G 
        		String original_letter = get_source(our_seq,result.name);
        		int place = place_in_sequence(result.name);
        		String mute_letter = get_mute_letter(result.name);
        		String new_name = original_letter + place + mute_letter;
        		result.set_full_name(new_name);
        
          		result.dotBracket= tokens.nextToken();
        		result.deltaG = tokens.nextToken();
        		result.set_new_Energy();
        		
        		result.eigValue = tokens.nextToken();
        		result.numOfVertices = tokens.nextToken();
        		result.repr = tokens.nextToken();


			    // calculate distance 
			    double newDistance;
			    double newDbDistance ;
			    
			    // Shapiro's distance
			    if (the_wild_type_representation.equals(result.repr))
			         newDistance = 0;     	
			    else 
			    	 newDistance = calculateDistance_byShapiro(the_wild_type_representation,result.repr);
		
           		result.setShapiroDistance(newDistance);
	
           		mut_counter++;
           		sum = sum+newDistance;
           		
           		// DotBracket distance
           		if (the_wild_type_repr_dot_b.equals(result.dotBracket))
           			newDbDistance = 0;
           		else
        			newDbDistance = calculateDistance_byDotB(the_wild_type_repr_dot_b,result.dotBracket);
		
                result.setDbDistance(newDbDistance);
	
                db_sum = db_sum+newDbDistance;
            
    		///the e.r. is ready to be insert to allEigResults
                allEigResults.addElement(result);
		allEigResults2.addElement(result);
             
                		
        		for (int i=0; i<finalResults.size(); i++){
        			if (((MuteCell) finalResults.elementAt(i)).eigValue.equals(result.eigValue)
        				&& (((EigResult)((MuteCell) finalResults.elementAt(i)).eigResults.elementAt(0)).numOfVertices.equals(result.numOfVertices)  )){
        				flag = false;
           				((MuteCell) finalResults.elementAt(i)).eigResults.addElement(result);
	           			if (result.numOfVertices.equals(numOfVerticesOfWT)) {
	        			   result.WT = "*";
	        			}
        			}
        		}
        		if (flag){ //we didn't find MuteCell with this EigValue
        			if (result.numOfVertices.equals(numOfVerticesOfWT)){
        			   result.WT = "*";
           			}
	        		mc = new MuteCell(result.eigValue);
	        		mc.eigResults.addElement(result);
	        		finalResults.addElement(mc);
        		}	
        	}
        	
     
   		    

		    //after reading the file and axtracting all neccesery data
        	//all the mutations are in the vector "allEigResults".
        	//now I need to sort it by distance
        	
        	//SORTING..
       		Vector allERshapiro = Sort.mergesort_by_distance(allEigResults); 
       		Vector allERdb = Sort.mergesort_by_db_distance(allEigResults2); 


     
        	//after reading the file and axtracting all neccesery data
        	//all the mutations are in the vector "finalResults".
        	average_distance = sum/mut_counter;
        	db_average_distance = db_sum/mut_counter;

                		
		string_average_distance = ((new Double (average_distance)).toString());
        	db_string_average_distance = ((new Double (db_average_distance)).toString());
		if (string_average_distance.length() > 6)
        		string_average_distance = string_average_distance.substring(0,6);
		if (db_string_average_distance.length() > 6)
        		db_string_average_distance = db_string_average_distance.substring(0,6);

        	// writing the firs HTML page: RESULT_TABLE.html
        	try {    	
       	 		BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
                out.write(  "<HTML><HEAD><TITLE> RESULT TABLE </TITLE></HEAD><BODY BGCOLOR=#cccff><H3> The results for your query:\n"+"</H3><p>"+
                			"The wild type sequence : " + our_seq + 
		                	"<br><br><div align=\"center\">"+
		                	"<p><p>");
            	out.write(  "<H4> For more information about a specific eigenValue, press the required link<H4>"+
            				"<p><p>"+
				 			"<TABLE BORDER=1 WIDTH=80% rules=all><THEAD>"+
				 			"<TR ALIGN=CENTER BGCOLOR=#ff66ff><TH>Second eigenvalue</TH>"+
				 			"<TH>Number of vertices</TH>"+
							"<TH>wild type</TH>"+
							"<TH>Frequency</TH>"+
							"</TR>");
	        
	        	int counter=1;
	        	double minDistance=100000;
	        	double maxDistance=0;
			double minDbDistance=100000;
	        	double maxDbDistance=0;
	        	while(finalResults.size()!=0){
	        		MuteCell minResult = ((MuteCell) finalResults.elementAt(0));
	       			for (int i=1; i<finalResults.size(); i++)
	        			if (Double.parseDouble(minResult.eigValue)>
	        				Double.parseDouble(((MuteCell)finalResults.elementAt(i)).eigValue))
					minResult = ((MuteCell)finalResults.elementAt(i));
					//sorts the muteCell by distance.
	    	        Vector t = Sort.mergesort_by_distance(minResult.eigResults);
	    	        minResult.eigResults = t;

	                
	        		printLineInMainTable(out,adress+"eigen" + counter +".html",minResult,numOfVerticesOfWT);
	        
	        		createSubTables(adress +"eigen" + counter +".html", minResult,minResult.eigValue,our_seq,
	        						((EigResult)minResult.eigResults.elementAt(0)).numOfVertices,
	        						 the_wild_type_eigenValue, the_wild_type_freeEnergy, 
	        						 string_average_distance,db_string_average_distance,the_wild_type_representation , the_wild_type_repr_dot_b);
	        	
	        	   //updates the minDistance value, and the maxDistanceValue.
	            	double new_min = ((EigResult)minResult.eigResults.elementAt(0)).distance;
					double new_db_min = ((EigResult)minResult.eigResults.elementAt(0)).dbDistance;
	        		if (new_min < minDistance) minDistance = new_min;
					if (new_db_min < minDbDistance) minDbDistance = new_min;
	           		int length  = minResult.eigResults.size(); //the length of the eigResults vector of this MuteCell
	        		double new_max = ((EigResult)minResult.eigResults.elementAt(length-1)).distance;
					double new_db_max = ((EigResult)minResult.eigResults.elementAt(length-1)).dbDistance;
			       	if (new_max > maxDistance) maxDistance = new_max;
	        		if (new_db_max > maxDbDistance) maxDbDistance = new_db_max;
	        	
	        	    finalResults_backup.addElement(minResult);
	        		finalResults.remove(minResult);
	        		counter++;
	        	}
	        	out.write("</TABLE><p><p>");
	       		out.write("<br><p><H4> For more information about the mutations in a specific distance range, press the required link</H4>");
	  			out.write("</TABLE> </HTML>" );
	      		printFourthTable(out,allERdb,minDbDistance,maxDbDistance,adress,db_string_average_distance,dbJump);//<--TO ADD THESE PARAMETERS
	       		out.write("</TABLE>  * the average of Dot-Bracket Distances is: "+ db_string_average_distance  );
       			out.write("<br>Clustering Resolution : "+dbJump+"                          </HTML>");
	       		printThirdTable(out,allERshapiro,minDistance,maxDistance,adress,string_average_distance,shapiroJump);
	      		out.write("</TABLE>  * the average of Shapiro's Distances is: "+ string_average_distance  );
	    		out.write("<br>Clustering Resolution : "+shapiroJump+"                          </HTML>");
	        	out.close();
	        	

	        }
	        catch(IOException e){System.out.println(e);}        	
	    }
	    catch(Exception e){System.out.println("Error "  + e);} 
	       
	    try{		
			Process p = Runtime.getRuntime().exec(new String[] {  "/bin/sh", "-c", "mv ./psPictures ./htmlDir/ "});
	        p.waitFor();      
		}
		catch (Exception e){System.out.println(e);}
		try{		
			Process p = Runtime.getRuntime().exec(new String[] {  "/bin/sh", "-c", "mv ./jpgPictures ./htmlDir/ "});
	        p.waitFor();      
		}
		catch (Exception e){System.out.println(e);}



	
	} //END OF GO
	
/*****************************************************/
	
	//print a line in the main table in RESULT_TABLE.html
	public  static void printLineInMainTable(BufferedWriter outFile,String path,MuteCell mc,String numOfVerticesOfWT){
		try{
			String temp1 = "*";
			String temp2= ((EigResult)mc.eigResults.elementAt(0)).WT;
			String temp = "-";
			if (mc.containsWT == true) temp = "WT";
				else if (temp1 == temp2){
	    			temp="*";
	    		}
				outFile.write(	"<TBODY><TR ALIGN=CENTER >"+	
								"<TD>"+ "<A HREF ="+ path +">"+ mc.eigValue + " </A>"+"</TD>"+
								"<TD>"+((EigResult)mc.eigResults.elementAt(0)).numOfVertices +"</TD>"+
		                    	"<TD>" + temp + "</TD>"+
								"<TD>" + mc.eigResults.size() + "</TD>"+"</TR>");
				outFile.flush();
		}
		catch(IOException e){System.out.println("I HAVE A PROBLEM  - please help me");}
	}
	
/*****************************************************/	
	
	public static void createSubTables(String path, MuteCell mc,String eigVal,String our_seq, String numOfNods, 
									   String  the_wild_type_eigenValue, String  the_wild_type_freeEnergy, String string_average_distance,
									   String db_string_average_distance, String the_wild_type_representation, String the_wild_type_repr_dot_b){
	    String adress="htmlDir/";
     	String picturesDir="psPictures/";
     	String picturesDir2="jpgPictures/";
     	
		try{
	 		BufferedWriter out = new BufferedWriter(new FileWriter(path));
            out.write(  "<HTML><HEAD><TITLE> MUTATIONS TABLE </TITLE></HEAD><BODY BGCOLOR=#99ccff>"+
            			"<H4>The mutations with eigenvalue = " + "<strong>"+
            			  eigVal + " and " + numOfNods + " vertices </strong> </H4>"+"<p>"+
            		    "(The mutations are sorted according to their Shapiro-Distance)<p>"+
            		    "* the average of Shapiro's Distances is: "+ string_average_distance + "<p>"+
			 			"<TABLE BORDER=1 WIDTH=100% rules=all><THEAD>"+
			 			"<TR ALIGN=CENTER BGCOLOR=#9966ff>"+
			 			"<TH>Mutation</TH>"+
						"<TH>Shapiro-Distance*</TH>"+
						"<TH>Minimum Energy Kcal/mol </TH>"+
						"<TH width=300>Shapiro representation</TH>"+
						"</TR></HTML>");
			for (int i = 0; i< mc.eigResults.size();i++){	
				EigResult current_mut =	((EigResult)mc.eigResults.elementAt(i));
				if (current_mut.full_name.equals("WT")){
					out.write(	"<TBODY><TR ALIGN=CENTER >"+	
								"<TD><strong>WT</strong></TD>"+
								"<TD>0.0</TD>"+
		                    	"<TD>" + current_mut.deltaG + "</TD>"+
		                    	"<TD>" +  current_mut.repr + "</TD>" );
				}
				else
      			out.write(	"<TBODY><TR ALIGN=CENTER >"+	
						    "<TD>"+ "<A HREF = "+ 	
							 current_mut.full_name +".html" +">"+
							 current_mut.full_name + " </A>"+"</TD>"+
							"<TD>"+	current_mut.distance +"</TD>"+
	                    	"<TD>" + current_mut.deltaG + "</TD>"+
	                    	"<TD>" +  current_mut.repr + "</TD>"+ "</TR>");			
				BufferedWriter out2 = new BufferedWriter(new FileWriter( adress + current_mut.full_name+".html"));
				out2.write( "<HTML><HEAD><TITLE>MUTATION  " +  current_mut.full_name + "</TITLE></HEAD><BODY BGCOLOR=#CC99CC>"+
                            "<H2>" + current_mut.full_name + "</H2>");
				out2.write( "<div align=\"left\">" + "   <H4><B> wild type "+
							"&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; "+
							"&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; "+
							"&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; "+
							"&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; "+
							
						
							"mutation </H4></B>" + "</div>" + 	
					
					
					    	"<IMG src=" +  picturesDir2 + "wild_typ.jpg  width=350 height=350 align=left border=3>  "+
					    	"<IMG src="+ picturesDir2 +  current_mut.name + ".jpg width=350 height=350 align=right border=3> " );
			   out2.write( "<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>"+
					    	" If you can't see the pictures click here:<br> "+
					    	"<A HREF =" + picturesDir + "wild_typ.ps target=_blank> Wild-Type </A><br>" +
					    	"<A HREF ="  + picturesDir + current_mut.name + ".ps target=_blank> Mutation &nbsp;" + current_mut.full_name + " </A><br>");
				String temp  = 	set_new_sequence(current_mut.name,our_seq);	   
				
				out2.write( "<H3><B><U> Information About This Mutation:</H3></B></U>") ;
				out2.write( "<li> Wild-Type Sequence:&nbsp;" + "<br>  &nbsp;&nbsp;"+our_seq+"<P>" ); 
				out2.write( "<li> Mutation &nbsp;Sequence:  &nbsp;&nbsp;" + "<br>  &nbsp;&nbsp;" + set_new_sequence(current_mut.name,our_seq) + "<P>");
				out2.write( "<li> Wild-Type EigenValue :&nbsp;&nbsp;&nbsp;" + the_wild_type_eigenValue + 
							"<li> Mutation EigenValue :&nbsp;&nbsp;&nbsp;&nbsp;" + current_mut.eigValue + "<P>"+
							"<li> Wild-Type Free Energy :&nbsp;&nbsp;&nbsp;" + the_wild_type_freeEnergy + " Kcal/mol"+
							"<li> Mutation Free Energy :&nbsp;&nbsp;&nbsp;&nbsp;" + current_mut.deltaG + " Kcal/mol"+"<P>"+
							"<li> Wild-Type Shapiro Representation:&nbsp;&nbsp; " + the_wild_type_representation +
							"<li> Mutation Shapiro Representation:&nbsp;&nbsp;&nbsp; " + current_mut.repr  + 
							"<li> Mutation Shapiro Distance:&nbsp;&nbsp;&nbsp;" + current_mut.distance + "&nbsp;&nbsp;" + "( avg:  "+string_average_distance+")<P>"+
							"<li> Wild-Type Dot-Bracket Representation:&nbsp;&nbsp; " + the_wild_type_repr_dot_b +
							"<li> Mutation Dot-Bracket Representation:&nbsp;&nbsp;&nbsp;&nbsp;" + current_mut.dotBracket +
						    "<li> Mutation Dot-Bracket Distance:&nbsp;&nbsp;&nbsp;" + current_mut.dbDistance+ "&nbsp;&nbsp;" + "( avg:  "+db_string_average_distance+")<P>");
							out2.close();
			}
			out.close();
		}
		catch(IOException e){System.out.println(e);}		 
	}

   /*****************************************************/
   
    //finds the place of the mutaion in the sequence	
	public static int place_in_sequence(String mute_name){
		int border = mute_name.lastIndexOf("-");
		String newLetter = mute_name.substring(0,border);
			return Integer.parseInt(newLetter);
	}
	
   /*****************************************************/
   
   // create the muted sequence 	
	public static String  set_new_sequence(String mute_name, String wild_type_seq){
		if (mute_name.length()>0){
			int border = mute_name.lastIndexOf("-");
			String newLetter = mute_name.substring(border+1);
			String prefix = wild_type_seq.substring(0,place_in_sequence(mute_name)-1);
			String suffix = wild_type_seq.substring(place_in_sequence(mute_name));
			return ( prefix +"<font color=red><B>"+ newLetter + "</B></font>" + suffix);
		}
		else return "";		
	}
	
	/*****************************************************/
	
	 //returns the point mutation letter 
	 public static String get_mute_letter(String mute_name){
		int border = mute_name.lastIndexOf("-");
		String newLetter = mute_name.substring(border+1,mute_name.length());
		return newLetter;
	}

   	/*****************************************************/
   	
	//returns the original letter in the wild type.
	public static String get_source(String wild_type, String mute_name){
		int place = place_in_sequence(mute_name) - 1;
		char a = wild_type.charAt(place);
		Character cr = new Character (a);
		String ans = cr.toString();
		return  ans;
	}
	
 
	
   	/*****************************************************/
	
	
	  public static void printThirdTable(BufferedWriter outFile,Vector allERshapiro,
										double minDistance,double maxDistance,String path,String string_average_distance,int shapiroJump ){
				double [] 	values2;										
				Vector ans = new Vector ();//in the end it will hold vectors which are the gruops by shapiro dis
				Vector current_gruop = new Vector ();
				int size = allERshapiro.size();
				values2 = new double [size];// holds for each gruop it's min distance
				double minDis = minDistance;
				double curr_distance = minDistance;
				double hefresh = 0;
				int max_jump = shapiroJump;//the jamp can be changed. it's for the gruoping by distance
				int group_counter = 1;
				values2[group_counter] = curr_distance;
				for (int i = 0;i< allERshapiro.size();i++){
						EigResult current_eigResult = ((EigResult)allERshapiro.elementAt(i));
						minDistance = curr_distance;
						curr_distance = current_eigResult.distance;
						
						hefresh = curr_distance - minDistance;
						if(hefresh<max_jump){
						current_gruop.addElement(current_eigResult);
						}
						else{//the jump is too big, open a new group
							ans.addElement(current_gruop);//insert to ans the first group
							current_gruop = new Vector ();//open a new group
							current_gruop.addElement(current_eigResult);//insert the eigResult  to the new group	
							group_counter++;// the num of groups increaces
							values2[group_counter] = curr_distance;
							}
				
				}// now all the mutations divided to groups acording to the given max_jump 
				 //and now we can print the third table
				 ans.addElement(current_gruop);
			////////////////////////////////////////////////	 
				 
				 
					try{
			outFile.write (	"<p><TABLE BORDER=1 WIDTH=50% rules=all><THEAD>"+
			 				"<TR ALIGN=CENTER BGCOLOR=#ff66ff>"+
			 				"<TH >Shapiro-Distances range</TH>"+
			 				"<TH>Frequency</TH>"+
							"</TR>" +			
							"<TBODY><TR ALIGN=CENTER >");
			for(int j = 0; j<ans.size();j++){		  		
				Vector v = (Vector)ans.elementAt(j);//takes the j-th vector
				outFile.write(  "<TBODY><TR ALIGN=CENTER >"+
							   	"<TD>"+ "<A HREF =" + path + "distance" + j +
								".html>" + 	((EigResult)v.elementAt(0)).distance + " - " + ((EigResult)v.elementAt(v.size() - 1)).distance + "</A> </TD>"+
								"<TD>" + v.size() + "</TD>"	);
				BufferedWriter out2 = new BufferedWriter(new FileWriter(path+ "distance" + (j) +".html"));//starting to write the page that shel open
			    out2.write( "<HTML><HEAD><TITLE> RESULT TABLE </TITLE></HEAD><BODY BGCOLOR=#99ccff><H4> Mutations in the Shapiro-Distance range of "+
			    			"&nbsp;" + 	((EigResult)v.elementAt(0)).distance+"-"+ ((EigResult)v.elementAt(v.size() - 1)).distance + 
			    			":\n"+"</H4><p>"+
			    			"(THE mutations are sorted according to their Shapiro-Distance)<p>"+
			    			"* the average of Shapiro's Distances is: "+ string_average_distance +
	                        "<br><br><div align=\"center\">"+
			                "<TABLE BORDER=1 WIDTH=100% rules=all><THEAD>"+
			            	"<TR ALIGN=CENTER BGCOLOR=#9966ff>"+
			            	"<TH>mutation</TH>"+	
			            	"<TH>Second eigenvalue</TH>"+
			            	"<TH>vertices num.</TH>"+
			            	"<TH>Shpiro-Distance* </TH>"+
			            	"<TH>minimum Energy Kcal/mol</TH>"+ 
			            	"<TH>Shapiro representation</TH>"+
					    	"</TR></HTML>");
				 for (int i = 0; i< v.size();i++){	
					EigResult current_mut = (EigResult)v.elementAt(i);
					if (current_mut.full_name.equals("WT")){
						out2.write(	"<TBODY><TR ALIGN=CENTER >"+	
									"<TD><strong>WT</strong></TD>"+
									"<TD>"+ current_mut.eigValue +"</TD>"+
				                    "<TD>" + current_mut.numOfVertices + "</TD>"+
				                    "<TD>0.0</TD>"+
								    "<TD>" + current_mut.deltaG + "</TD>"+
				                    "<TD>" +  current_mut.repr + "</TD>"+	
									"</TR>");
					}
					else
					out2.write(	"<TBODY><TR ALIGN=CENTER >"+	
								"<TD>"+ "<A HREF = "+
								 current_mut.full_name +".html" +">"+
								 current_mut.full_name + " </A>" + "</TD>"+
								"<TD>"+ current_mut.eigValue +"</TD>"+
								"<TD>"+ current_mut.numOfVertices +"</TD>"+
			                    "<TD>" + current_mut.distance + "</TD>"+
							    "<TD>" + current_mut.deltaG + "</TD>"+
			                    "<TD>" +  current_mut.repr + "</TD>"+	
								"</TR>");
				}
				out2.close();
			}					
		}catch(Exception e){System.out.println("Error while trying to do the THIRD table " +e);}						 
				 
				 
				
				}		
	
	/*****************************************************/
		
	  public static void printFourthTable(BufferedWriter outFile,Vector allERdb,
										double minDistance,double maxDistance,String path,String string_average_distance,int dbJump ){
				double []   values3;										
				Vector ans = new Vector ();//in the end it will hold vectors which are the gruops by shapiro distance
				Vector current_gruop = new Vector ();
				int size = allERdb.size();
				values3 = new double [size];// holds for each gruop it's min distance
				double minDis = minDistance;
				double curr_distance = minDistance;
				double hefresh = 0;
				int max_jump = dbJump;//the jamp can be changed. it's for the gruoping by distance
				int group_counter = 1;
				values3[group_counter] = curr_distance;
				for (int i = 0;i< allERdb.size();i++){
						EigResult current_eigResult = ((EigResult)allERdb.elementAt(i));
						minDistance = curr_distance;
						curr_distance = current_eigResult.dbDistance;
						
						hefresh = curr_distance - minDistance;
						if(hefresh<max_jump){
						current_gruop.addElement(current_eigResult);
						}
						else{//the jump is too big, open a new group
							ans.addElement(current_gruop);//insert to ans the first group
							current_gruop = new Vector ();//open a new group
							current_gruop.addElement(current_eigResult);//insert the eigResult  to the new group	
							group_counter++;// the num of groups increaces
							values3[group_counter] = curr_distance;
							}
				
				}// now all the mutations divided to groups acording to the given max_jump 
				 //and now we can print the third table
				 ans.addElement(current_gruop);
			////////////////////////////////////////////////	 
				 
				 
					try{
			outFile.write (	"<p><TABLE BORDER=1 WIDTH=50% rules=all><THEAD>"+
			 				"<TR ALIGN=CENTER BGCOLOR=#ff66ff>"+
			 				"<TH >Dot-Bracket Distance range</TH>"+
			 				"<TH>Frequency</TH>"+
							"</TR>" +			
							"<TBODY><TR ALIGN=CENTER >");
			for(int j = 0; j<ans.size();j++){		  		
				Vector v = (Vector)ans.elementAt(j);//takes the j-th vector
				outFile.write(  "<TBODY><TR ALIGN=CENTER >"+
							   	"<TD>"+ "<A HREF =" + path + "DBdistance" + j +
								".html>" + 	((EigResult)v.elementAt(0)).dbDistance + " - " + ((EigResult)v.elementAt(v.size() - 1)).dbDistance + "</A> </TD>"+
								"<TD>" + v.size() + "</TD>"	);
				BufferedWriter out2 = new BufferedWriter(new FileWriter(path+ "DBdistance" + (j) +".html"));//starting to write the page that shell open
			    out2.write( "<HTML><HEAD><TITLE> RESULT TABLE </TITLE></HEAD><BODY BGCOLOR=#99ccff><H4> Mutations in the Dot-Bracket Distance range of "+
			    			"&nbsp;" + 	((EigResult)v.elementAt(0)).dbDistance+"-"+ ((EigResult)v.elementAt(v.size() - 1)).dbDistance + 
			    			":\n"+"</H4><p>"+
			    			"(THE mutations are sorted according to their Dot-Bracket Distance)<p>"+
			    			"* the average of Dot-Bracket Distance  is: "+ string_average_distance +
	                        "<br><br><div align=\"center\">"+
			                "<TABLE BORDER=1 WIDTH=100% rules=all><THEAD>"+
			            	"<TR ALIGN=CENTER BGCOLOR=#9966ff>"+
			            	"<TH>mutation</TH>"+	
			            	"<TH>Second eigenvalue</TH>"+
			            	"<TH>vertices num.</TH>"+
			            	"<TH>Dot-Bracket Distance* </TH>"+
			            	"<TH>minimum Energy Kcal/mol</TH>"+ 
			            	"<TH>Dot-Bracket representation</TH>"+
					    	"</TR></HTML>");
				 for (int i = 0; i< v.size();i++){	
					EigResult current_mut = (EigResult)v.elementAt(i);
					if (current_mut.full_name.equals("WT")){
						out2.write(	"<TBODY><TR ALIGN=CENTER >"+	
									"<TD><strong>WT</strong></TD>"+
									"<TD>"+ current_mut.eigValue +"</TD>"+
				                    "<TD>" + current_mut.numOfVertices + "</TD>"+
				                    "<TD>0.0</TD>"+
								    "<TD>" + current_mut.deltaG + "</TD>"+
				                    "<TD>" +  current_mut.dotBracket + "</TD>"+	
									"</TR>");
					}
					else
					out2.write(	"<TBODY><TR ALIGN=CENTER >"+	
								"<TD>"+ "<A HREF = "+
								 current_mut.full_name +".html" +">"+
								 current_mut.full_name + " </A>" + "</TD>"+
								"<TD>"+ current_mut.eigValue +"</TD>"+
								"<TD>"+ current_mut.numOfVertices +"</TD>"+
			                    "<TD>" + current_mut.dbDistance + "</TD>"+
							    "<TD>" + current_mut.deltaG + "</TD>"+
			                    "<TD>" +  current_mut.dotBracket + "</TD>"+	
								"</TR>");
				}
				out2.close();
			}					
		}catch(Exception e){System.out.println("Error while trying to do the Fourth table " +e);}						 
				 
				 
				 
				 
				
	
				
				}		
	
	/*****************************************************/
	
	
	
 
 	//removes "R"s from the shapirro representation	
	public static String removeR (String s){
		int length = s.length();
		String ans = s.substring(1,length-2);
		return ans;
	}
	
	/*****************************************************/
		
	public static double calculateDistance_byShapiro(	String Original_WT, String Original_mute){
		try{
			String WT = removeR(Original_WT);
			String mute = removeR(Original_mute);
                        System.out.println("wt:" + WT);
                        System.out.println("mute:" + mute);
			BufferedWriter out = new BufferedWriter(new FileWriter("temporary.txt"));
			out.write(WT);
			out.write("\n");
			out.write(mute);
			out.close();
			Process p = Runtime.getRuntime().exec(new String[] {  "/bin/sh", "-c", "bin/RNAdistance  -Dw < temporary.txt "});
		    p.waitFor();
	        BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
	        BufferedReader br1 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	        String line="";
	        String ans =  br.readLine();
	        String toReturn = ans.substring(3);
                System.out.println("ans:" + ans);
		    return Double.parseDouble(toReturn);
		}
		catch (Exception e){System.out.println(e+" while calculating Shapiro Distance");}
		return -5;
	}
	
	/*****************************************************/
	
	public static double calculateDistance_byDotB(	String WT, String mute){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter("temporary1.txt"));
			out.write(WT);
			out.write("\n");
			out.write(mute);
			out.close();
			Process p = Runtime.getRuntime().exec(new String[] {  "/bin/sh", "-c", "bin/RNAdistance   < temporary1.txt "});
		    p.waitFor();
	       
	        BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
	        BufferedReader br1 = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	        String line="";
	        String ans =  br.readLine();
	        String toReturn = ans.substring(3);
		    return Double.parseDouble(toReturn);
		}
		catch (Exception e){System.out.println(e +"while calculating Dot-Bracket Distance");}
		return -5;
	}
	
        public static int stringToInt(String str, int defaultValue) {
        try {
            //return Integer.parseInt(str);
            Integer i = new Integer(str);
            return i.intValue();
        } catch (NumberFormatException nfe) {
        	System.out.println("in here");
            return defaultValue;
        }
    }


}// END OF CLASS






