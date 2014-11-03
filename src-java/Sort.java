/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 *
 *									class Sort
 *
 *  this class have static function for sorting our vectors.
 *  using mergeSort algorithm.
 *	there are three different types sort: 
 *	1.by the Shapiro distance value of the mutations - mergesort_by_distance(Vector v)
 *	2.by Dot Bracket distance value of the mutations-  mergesort_by_db_distance(Vector v)
 *	3.by_eigen value of the mutations-  mergesort_by_eigen(Vector v)
 *
 *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/ 

import java.util.*;

public class Sort{
	
	public static Vector mergesort_by_distance(Vector v){
	    if (v.isEmpty() || v.size()==1) {
	        return  v; 
	    }
	    else{
	       Vector [] pair = split(v);
	       Vector ans = 
	       merge_by_distance(mergesort_by_distance(pair[0]), mergesort_by_distance(pair[1]));
	       return ans;
	    }
    }
    
    public static Vector mergesort_by_db_distance(Vector v){
	    if (v.isEmpty() || v.size()==1) {
	        return  v; 
	    }
	    else{
	       Vector [] pair = split(v);
	       Vector ans = 
	       merge_by_db_distance(mergesort_by_db_distance(pair[0]), mergesort_by_db_distance(pair[1]));
	       return ans;
	    }
    }
    
            
    public static Vector merge_by_distance(Vector v1, Vector v2){
        Vector vec;
        EigResult head;
        if (v1.isEmpty()) 
            return v2;
        else if (v2.isEmpty()) 
            return v1;
        else {
            double data1 = (double)((EigResult) v1.elementAt(0)).distance;
            double data2 = (double)((EigResult) v2.elementAt(0)).distance;
          
            if (data1< data2) 
                head = (EigResult)v1.remove(0);
            else 
                head = (EigResult)v2.remove(0);
             
            vec = merge_by_distance(v1,v2);

            vec.add(0,head);
            return vec;
        }
    }
    
     public static Vector merge_by_db_distance(Vector v1, Vector v2){
        Vector vec;
        EigResult head;
        if (v1.isEmpty()) 
            return v2;
        else if (v2.isEmpty()) 
            return v1;
        else {
            double data1 = (double)((EigResult) v1.elementAt(0)).dbDistance;
            double data2 = (double)((EigResult) v2.elementAt(0)).dbDistance;
          
            if (data1< data2) 
                head = (EigResult)v1.remove(0);
            else 
                head = (EigResult)v2.remove(0);
             
            vec = merge_by_db_distance(v1,v2);

            vec.add(0,head);
            return vec;
        }
    }


    public static Vector[] split (Vector v){
        Vector[] pair = {new Vector(),new Vector()};
        int i = 0;
        while ( !v.isEmpty() ){
            pair[i%2].add(0,v.remove(0));
            i = i+1;
        }
        return pair;
    }
    
    
    public static Vector mergesort_by_eigen(Vector v){
        if (v.isEmpty() || v.size()==1) {
            return  v; 
        }
        else{
            Vector [] pair = split(v);
            Vector ans = 
           merge_by_eigen(mergesort_by_eigen(pair[0]), mergesort_by_eigen(pair[1]));
        
           return ans;
        }
    }
            
    public static Vector merge_by_eigen(Vector v1, Vector v2){
        Vector vec;
        EigResult head;
        if (v1.isEmpty()) 
            return v2;
        else if (v2.isEmpty()) 
            return v1;
        else {
            double data1 = Double.parseDouble(((EigResult) v1.elementAt(0)).eigValue);
            double data2 = Double.parseDouble(((EigResult) v2.elementAt(0)).eigValue);
          
            if (data1< data2) 
                head = (EigResult)v1.remove(0);
            else 
                head = (EigResult)v2.remove(0);
             
            vec = merge_by_eigen(v1,v2);
            vec.add(0,head);
            return vec;
        }
    }

	
	public static void sort_by_eigenValue(Vector v){
		int min;
		EigResult temp;
		for(int index=0;index<v.size()-1;index++){
			min=index;
			for(int scan=index+1;scan<v.size();scan++){
				double eig1 = Double.parseDouble(((EigResult)v.elementAt(scan)).eigValue);
				double eig2 = Double.parseDouble(((EigResult)v.elementAt(min)).eigValue);
				if (eig1<eig2)
					min=scan;
			}
			temp= (EigResult)v.elementAt(min);
			v.setElementAt((EigResult)v.elementAt(index),min);
			v.setElementAt(temp,index);
		}
	}

}