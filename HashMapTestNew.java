import java.util.HashMap;

import mosek.Env;
import mosek.Task;

import java.util.*;

public class HashMapTestNew {
	
	public static void main(String[] args) {
		HashMap<Integer, Double> hm1 = new <Integer, Double> HashMap();
		HashMap<Integer, Double> hmresult1 = new <Integer, Double> HashMap();
		
		HashMapTestNew hmt = new HashMapTestNew();
	    hm1.put(1, 10.0);
	    hm1.put(2, 10.0);
	    hm1.put(3, 10.0);
	    hm1.put(4, 10.0);
	    hm1.put(5, -40.0);
	    
		
		hmresult1 = hmt.mister(hm1);
		System.out.println("The expected result is ");
		for (int q = 0; q < hmresult1.size(); q++) {
			  System.out.println(hmresult1.get(q+1));
		  }
		
		
	}
	
	public HashMap mister(HashMap m)
	  {
	    // Since the value infinity is never used, we define
	    // 'infinity' symbolic purposes only
		  
		HashMap<Integer, Double> hm1 = m;
		HashMap<Integer, Double> hmresult = new <Integer, Double> HashMap(); 
		
		List<Integer> excess = new ArrayList<Integer>();
		List<Integer> deficit = new ArrayList<Integer>();
		List<Integer> neutral = new ArrayList<Integer>();
		
		int i = 0;
		int j = 0;
		int poscount = 0;
		int negcount = 0;
		int zerocount = 0;
		
		for (j = 1; j < hm1.size() + 1 ; j++) {
			double check = ((Double)hm1.get(j)).doubleValue();
			if (check > 0 ) {
				poscount += 1;
				//excess.put(j, check);
				excess.add(j);
			}
			else if (check < 0) {
				negcount += 1;
				//deficit.put(j, check);
				deficit.add(j);
			}
			else {
				zerocount += 1;
				//neutral.put(j, check);
				neutral.add(j);
			}
		}
		
		//int numvar = hm1.size();
		final int numvar = excess.size();
		final int numcon = deficit.size();
		int asub[][] = new int [numcon][excess.size()];
		double aval[][] = new double [numcon][excess.size()];
		int d = 0;
		int e = 0;
		int f = 0;
		double blcvalue = 0.0;
		double bucvalue = 0.0;
		
		//double c1[] = new double[numvar];
		double c[] = new double[numvar];
		
		for (j = 0; j < numvar; j++) {
			d = excess.get(j);
			c[j] = ((Double)hm1.get(d)).doubleValue();
		}
		
		for (i = 0; i < numcon; i++) {
			for (j = 0; j < excess.size(); j++) {
				asub[i][j] = j;
				d = excess.get(j);
				aval[i][j] = ((Double)hm1.get(d)).doubleValue();
			}
		}
		
		double blc[] = new double[numcon];
		double buc[] = new double[numcon];
		double blx[] = new double[numvar];
		double bux[] = new double[numvar];
		
		for (i = 0; i < numcon; i++) {
			e = deficit.get(i);
			blcvalue = (-1.0)*((Double)hm1.get(e)).doubleValue();
			blc[i] = blcvalue;
			buc[i] = blcvalue;
		}
		
		for (j = 0; j < numvar; j++) {
			blx[j] = 0;
			bux[j] = 1;
		}
		
		
		  //final int numcon = 4;
		  //final int numvar = 8;
		  final int NUMANZ = 9;
		  double result[] = new double[8]; 
		  
		  double
	      infinity = 0;
	    
	    //double c[]    = {3.0, 3.0, 1.0, 5.0, 2.0, 0.0, 0.0, 3.0};
	    /*int    asub[][] = { {0},
	                        {0,1,2,3,5},
	                        {0,3,4,6},
	                        {1,5,6,7} };*/
	    /*double aval[][] = { {1.0},
	                        {1.0, -1.0, -1.0, 1.0, 1.0},
	                        {1.0, -1.0, -1.0, 1.0},
	                        {1.0, -1.0, -1.0, -1.0}};*/
	    
		  mosek.boundkey[] bkc = new mosek.boundkey[numcon];
		    
		  for (i = 0; i < numcon; i++) {
		    	bkc[i] = mosek.boundkey.fx;
		  }
		  
		 /*mosek.boundkey[]
	                 bkc    = {mosek.boundkey.fx,
	                           mosek.boundkey.fx,
	                           mosek.boundkey.fx,
	                           mosek.boundkey.fx};*/
	    /*double  blc[]  = {1.0,
	                      0.0,
	                      0,0,
	                      0.0};
	    double  buc[]  = {1.0,
	                      0.0,
	                      0.0,
	                      0.0};*/
	    
		
		mosek.boundkey[] bkx = new mosek.boundkey[numvar];
		    
		for (i = 0; i < numvar; i++) {
		    	bkx[i] = mosek.boundkey.ra;
		}  
		  
		/*mosek.boundkey
	            bkx[]  = {mosek.boundkey.ra,
	                      mosek.boundkey.ra,
	                      mosek.boundkey.lo,
	                      mosek.boundkey.ra,
	                      mosek.boundkey.lo,
	                      mosek.boundkey.lo,
	                      mosek.boundkey.ra,
	                      mosek.boundkey.lo,};*/
	    /*double  blx[]  = {0.0,
	                      0.0,
	                      0.0,
	                      0.0,
	                      0.0,
	                      0.0,
	                      0.0,
	                      0.0};
	    double  bux[]  = {1.0,
	                      2.0,
	                      +infinity,
	                      2.0,
	                      +infinity,
	                      +infinity,
	                      2.0,
	                      +infinity};*/
	    double[] xx  = new double[numvar];

	    try (Env  env  = new Env(); 
	         Task task = new Task(env,0,0))
	    {
	      // Directs the log task stream to the user specified
	      // method task_msg_obj.stream
	      task.set_Stream(
	        mosek.streamtype.log,
	        new mosek.Stream() 
	          { public void stream(String msg) { System.out.print(msg); }});

	      /* Give MOSEK an estimate of the size of the input data. 
	     This is done to increase the speed of inputting data. 
	     However, it is optional. */
	      /* Append 'numcon' empty constraints.
	     The constraints will initially have no bounds. */
	      task.appendcons(numcon);
	      
	      /* Append 'numvar' variables.
	     The variables will initially be fixed at zero (x=0). */
	      task.appendvars(numvar);
	        
	      for(j=0; j<numvar; ++j)
	      {
	        /* Set the linear term c_j in the objective.*/  
	        task.putcj(j,c[j]);
	        /* Set the bounds on variable j.
	           blx[j] <= x_j <= bux[j] */
	        task.putbound(mosek.accmode.var,j,bkx[j],blx[j],bux[j]);
	      }
	      /* Set the bounds on constraints.
	       for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
	      for(i=0; i<numcon; ++i)
	      {  
	        task.putbound(mosek.accmode.con,i,bkc[i],blc[i],buc[i]);

	          /* Input row i of A */   
	          task.putarow(i,                     /* Row index.*/
	                       asub[i],               /* Column indexes of non-zeros in row i.*/
	                       aval[i]);              /* Non-zero Values of row i. */
	      }

	      /* A maximization problem */ 
	      task.putobjsense(mosek.objsense.maximize);

	      /* Solve the problem */
	      mosek.rescode r = task.optimize();
	                
	      // Print a summary containing information
	      //   about the solution for debugging purposes
	      task.solutionsummary(mosek.streamtype.msg);     
	       
	      mosek.solsta solsta[] = new mosek.solsta[1];
	      /* Get status information about the solution */ 
	      task.getsolsta(mosek.soltype.bas,solsta);

	      task.getxx(mosek.soltype.bas, // Basic solution.     
	                 xx);
	      System.out.println("The solsta is not " + solsta[0]);
	      switch(solsta[0])
	      {
	      case optimal:
	    	  //return null;
	      case near_optimal:
	        //System.out.println("Optimal primal solution\n");
	        //for(int j = 0; j < numvar; ++j)
	            //System.out.println ("x[" + j + "]:" + xx[j]);
	        	//result[j] = xx[j];
	        //return result[];
	    	for (j = 0; j < numvar; j++) {
	    		hmresult.put(j+1, xx[j]);
	    	}  
	    	return hmresult;  
	        //return xx;
	        //break;
	      case dual_infeas_cer:
	    	  return null;
	      case prim_infeas_cer:
	    	  return null;
	      case near_dual_infeas_cer:
	    	  return null;
	      case near_prim_infeas_cer:  
	        System.out.println("Primal or dual infeasibility.\n");
	        return null;
	        //break;
	      case unknown:
	        System.out.println("Unknown solution status.\n");
	        return null;
	        //break;
	      default:
	        System.out.println("Other solution status");
	        return null;
	        //break;
	      }
	    }
	    catch (mosek.Exception e1)
	    {
	      System.out.println ("An error/warning was encountered");
	      System.out.println (e1.toString());
	      throw e1;
	    }
	  }

}
