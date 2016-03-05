import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import dijkstra.engine.DijkstraAlgorithm;
import dijkstra.model.Edge;
import dijkstra.model.Graph;
import dijkstra.model.Vertex;
import jsc.distributions.Binomial;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
/*
* This is the main class.
* Use 'runSimulation("graphname")' to run a single simulation on the graph 'graphname' from the input folder.
*
* 'spsa("graphname",numberOfIterations)' runs a single simulation on the graph 'graphname' from the input folder where
* the penalty function is optimized using SPSA. 'numberOfIterations' is the total number of iterations to optimize the
* penalty function for each graph. The more is better but takes time...
*
* Use 'runFromTxt("textfilename",a)' to bulk run a list of graphs provided in textfilename.txt located in the input
* folder. 'a' is a binary variable. a = 0 uses dt penalty function (or whatever function you designate in 'getPenalty'.
* a = 1 uses an SPSA-obtained penalty function.
*
* You can use GenerateWeather.generateWeather("graphname",n) to generate weather (probability realization) for a graph
* file. 'n' denotes the total number of weather files that will be generated.
*
*
*
* NOTES:
* 1. Start node is the first indexed vertex and terminal node is the last indexed vertex.
* 2. Edge probabilities are interpreted as P(traversable).
* 3. Edge probabilities will be obtained: 
*   a. by reading the corresponding weather files IF obtainWeatherFromFile = true (for a.graph it's a.prb),
*   b. by randomly generating using Binomial Distribution IF obtainWeatherFromFile = false.
* 4. iterationNumber is by default 1. If obtainWeatherFromFile = false you can set it to whatever you like,
*   but if obtainWeatherFromFile = true, you will need the corresponding weather files to run different instances.
* 5. This code computes two paths:
*   a. PBA method: Using a penalty function from getPenalty(), weighted graph is obtained and then solved according to
*       the navigate-disambiguate-reset (NDR) strategy.
*   b. Optimism method: Without using any probabilistic aid, this algorithm assumes all edges are traversable and then
*       solves the graph according to the NDR strategy. Referred to as "Optimistic Algorithm" in Eyerich et. al. 2009.
* 6. Results will be printed to 'SimulationResults.csv' for standard runs. For SPSA runs you can find the results in
* 'SPSAResults.csv'.
*/
/*
* furkansahin@std.sehir.edu.tr
*/
public class SampleBasedCTP {
    //initializing...
    private List<Vertex> nodes;
    private List<Edge> edges;
    private List<Edge> edgesForOptimistic;
    private List<Vertex> nodesForOptimistic;
    private List<Edge> weightedEdges;
    private List<Vertex> weightedNodes;
    private List<Edge> prunedEdges;
    private List<Vertex> prunedNodes;
    double total = 0;
    double expected = 0;
    double expectedDistance = 0;
    String thePath = "1";//initializing the solution path.
    //initializing...

    int cost = 0;//cost of disambiguating an edge. default: 0.
	boolean obtainWeatherFromFile = true;//use '.prb' files or random probability realizations?
	//public static double slope = 10.0;

    public static double[] parameters = {4.9334, 0.6839, -29.5583, 6.0710, -0.02, 34.4602};
    //SPSA parameters. please refer to http://www.jhuapl.edu/spsa/ for more information about SPSA.
    //close-to-dt: {4.9334, 0.6839, -29.5583, 6.0710, -0.02, 34.4602}

	public static void main(String[] args) throws IOException {
		//GenerateWeather.generateWeather("delaunay50_100b_2",2);
        //runFromTxt("PBA",0);

        //SampleBasedCTP ctp = new SampleBasedCTP();
        //ctp.runSPSASimulation("grid10_100b_1",5);
        //spsa("delaunay50_100b_2",25);
        runFromTxt("PBA",0);

    }

    public static void runFromTxt(String filename, int a) throws IOException {
        int runNo = 1;
        String[] nextLine1;
        SampleBasedCTP ctp = new SampleBasedCTP();
        CSVReader reader1 = new CSVReader(new FileReader("input"+File.separator+filename+".txt"), (char)32);
        if (a==0) {
            while ((nextLine1 = reader1.readNext()) != null) {
                ctp.runSimulation(nextLine1[0]);
                System.out.println("-------------------------------------------------------"+runNo+"----------------------------------------------------------------");
                runNo++;
            }
        }

        if (a==1) {
            while ((nextLine1 = reader1.readNext()) != null) {
                spsa(nextLine1[0],25); // default: 100.
            }
        }

        reader1.close();
        System.out.println("Done.");
    }

    public static void writeToSPSAFile(double adaptive, double dt, long rta, long rtdt) throws IOException {
        CSVWriter writeSPSA = new CSVWriter(new FileWriter("SPSAResults.csv",true));
        String[] outputColumn = new String[4];

        outputColumn[0] = Double.toString(adaptive);
        outputColumn[1] = Double.toString(dt);
        outputColumn[2] = Long.toString(rta);
        outputColumn[3] = Long.toString(rtdt);

        writeSPSA.writeNext(outputColumn);

        writeSPSA.close();
    }

    public static double[] spsa(String filename, int numberOfIterations) throws IOException{
        int n = 100;
        int p = parameters.length;
        double[] delta = new double[p];
        double A = n*0.05; //default: n*0.1;
        double alpha = 0.602; //default: 0.602
        double gamma = 0.101; //default: 0.101
        //double a = ((A/n)*Math.pow(A+1,alpha))/5; //why 10? mean magnitude of elements of g0(theta0)??? scaling matrix: different ak for different coefficients of the penalty function?? default: ((A/n)*Math.pow(A+1,alpha))/10
        //double c = 20; //stdev of the noise in y(theta)??? default: 50
        double[] theta = parameters; //should be given
        double[] thetaplus = new double[p];
        double[] thetaminus = new double[p];
        double yplus = 0;
        double yminus = 0;
        double[] ghat = new double[p];


        //double[] a = {500,100,16,500,16,500};
        //double[] c = {0.5,10,50,0.5,50,0.5};
        double[] a = {500,1,1,1,1,500};
        double[] c = {0.1,1,1,0.1,1,0.1};
        double[] ak = new double[p];
        double[] ck = new double[p];

        SampleBasedCTP ctp = new SampleBasedCTP();

        GenerateWeather.generateWeather(filename,numberOfIterations); //generateweather burada olursa hep aynı instancelar üzerinden çalışır.

        long startTime = System.currentTimeMillis();

        for (int k = 0; k < n; k++) {
            for (int d = 0; d < p; d++) {
                ak[d] = a[d] / Math.pow((k+1+A),alpha);
                ck[d] = c[d] / Math.pow((k+1),gamma);
            }

            for (int d = 0; d < p; d++) {
                delta[d] = (int) (1000*Math.round(Math.random())-500);
                //delta[d] = Math.round(Math.random())-0.5;
                //delta[d] = (Math.round(Math.random())/10)-0.05;
            }

            for (int d = 0; d < p; d++) {
                thetaplus[d] = theta[d] + ck[d]*delta[d];
                thetaminus[d] = theta[d] - ck[d]*delta[d];
                //System.out.println(Arrays.toString(thetaplus));
            }

            yplus = ctp.lossFunction(filename, numberOfIterations, thetaplus);
            yminus = ctp.lossFunction(filename, numberOfIterations, thetaminus);

            for (int d = 0; d < p; d++) {
                ghat[d] = (yplus - yminus) / (2*ck[d]*delta[d]);
            }

            for (int d = 0; d < p; d++) {
                theta[d] = theta[d] - (ak[d]*ghat[d]);
            }

            System.out.println(k+": "+ctp.lossFunction(filename,1,theta)+" (dt) "+ctp.lossFunction(filename,5,theta)+" (adap.) "+Arrays.toString(theta));
            System.out.println(Arrays.toString(ghat));
        }
        long endTime = System.currentTimeMillis();
        long runTime = (endTime - startTime);
        System.out.println("SPSA done!-----------------------------------------------------------------------"+filename+"--------------------------------------------");
        //System.out.println(Arrays.toString(theta));
        //System.out.print(filename+" ");
        //System.out.print((-1)*ctp.lossFunction(filename,-1,theta)+" ");
        //System.out.println(ctp.lossFunction(filename, 1, theta)+" "+runTime);

        double dt = (-1)*ctp.lossFunction(filename,-1,theta);
        long endTime2 = System.currentTimeMillis();
        long runTime2 = (endTime2 - endTime);

        //writeToSPSAFile(ctp.lossFunction(filename,1,theta),dt,runTime,runTime2);

        return theta;
    }

    public double lossFunction (String filename, int numberOfIterations, double[] a) throws IOException {
        //Loss function is an SPSA term. Please refer to SPSA documentation for more information.
        double total = 0;

        //BURADA DA OLABİLİR GENERATEWEATHER

        nodes = new ArrayList<Vertex>();
        edges = new ArrayList<Edge>();

        CSVReader reader = new CSVReader(new FileReader("input"+File.separator+filename+".graph"), (char)32);
        String [] nextLine;
        nextLine = reader.readNext();
        int numberOfNodes = Integer.parseInt(nextLine[1]);
        int numberOfEdges = Integer.parseInt(nextLine[2]);
        for (int i = 1; i <= numberOfNodes; i++) {
            Vertex location = new Vertex(i-1, "Node_" + i);
            nodes.add(location);
        }
        int i=1;
        while ((nextLine = reader.readNext()) != null) {
            //the edges in the input file are being added into the graph.
            addLane(i-1, Integer.parseInt(nextLine[1])-1, Integer.parseInt(nextLine[2])-1, Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 0);
            //opposite direction of the edges above are being added into the graph (the graph has to be undirected).
            addLane(i, Integer.parseInt(nextLine[2])-1, Integer.parseInt(nextLine[1])-1, Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 0);
            i++;
        }
        reader.close();

        Graph graph = new Graph(nodes, edges, nodes.get(0), nodes.get(numberOfNodes-1));

        if (numberOfIterations==-1) {
            CSVReader reader2 = new CSVReader(new FileReader("input"+File.separator+filename+".prb"), (char)32);
            String [] nextLine2;
            i=0;
            while ((nextLine2 = reader2.readNext()) != null) {
                int s = Integer.parseInt(nextLine2[2]);
                weightedEdges.get(i).setActualStatus(s);
                weightedEdges.get(i+1).setActualStatus(s);
                weightedEdges.get(i).setWeight( edges.get(i).getLength() + getPenalty(graph, edges.get(i)) );
                weightedEdges.get(i+1).setWeight( weightedEdges.get(i).getWeight() );
                i=i+2;
            }
            reader2.close();
            Graph graphWeighted = new Graph(weightedNodes, weightedEdges, weightedNodes.get(0), weightedNodes.get(numberOfNodes - 1));
            total = computeExpectedLength(graphWeighted, graphWeighted.getStartNode());
        }
        for (int it = 0; it < numberOfIterations; it++) {
            double result;
            weightedEdges = new ArrayList<Edge>();
            weightedNodes = new ArrayList<Vertex>();
            weightedEdges = cloneEdges(edges);
            weightedNodes = cloneVertices(nodes);


            if (numberOfIterations==1) {
                CSVReader reader2 = new CSVReader(new FileReader("input"+File.separator+filename+".prb"), (char)32);
                String [] nextLine2;
                i=0;
                while ((nextLine2 = reader2.readNext()) != null) {
                    int s = Integer.parseInt(nextLine2[2]);
                    weightedEdges.get(i).setActualStatus(s);
                    weightedEdges.get(i+1).setActualStatus(s);
                    weightedEdges.get(i).setWeight( edges.get(i).getLength() + getSPSAPenalty(graph, edges.get(i), a) );
                    weightedEdges.get(i+1).setWeight( weightedEdges.get(i).getWeight() );
                    i=i+2;
                }
                reader2.close();
            }
            else {
                CSVReader reader2 = new CSVReader(new FileReader("input"+File.separator+"prb"+File.separator+filename+"_"+(it+1)+".prb"), (char)32);
                String [] nextLine2;
                i=0;
                while ((nextLine2 = reader2.readNext()) != null) {
                    int s = Integer.parseInt(nextLine2[2]);
                    weightedEdges.get(i).setActualStatus(s);
                    weightedEdges.get(i+1).setActualStatus(s);
                    weightedEdges.get(i).setWeight( edges.get(i).getLength() + getSPSAPenalty(graph, edges.get(i), a) );
                    weightedEdges.get(i+1).setWeight( weightedEdges.get(i).getWeight() );
                    i=i+2;
                }
                reader2.close();
            }

            Graph graphWeighted = new Graph(weightedNodes, weightedEdges, weightedNodes.get(0), weightedNodes.get(numberOfNodes - 1));
            result = computeExpectedLength(graphWeighted, graphWeighted.getStartNode());

            total = total + result;
        }
        return total/numberOfIterations;
    }

	public double getPenalty(Graph graph, Edge e) {
        //Weight function.
		double penalty = 0;
		double dT;
		double traversableProb = e.getProb();//traversableProb is (1-rho).
		if(traversableProb==0) {
			traversableProb = 0.00001;
		} else if(traversableProb!=1) {
			//instead of obtaining the Euclidian distance to the termination node like solving an SOSP, in classic CTP
			//we calculate DT by executing Dijkstra to get the shortest DETERMINISTIC distance to the termination
			//point from the end of the edge (if it's the final edge then it is calculated from the start).
			DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
			if (e.getDestination().getId()+1==graph.getVertexes().size()){
				dijkstra.execute(e.getSource());
			}
			else {
				dijkstra.execute(e.getDestination());
			}
			dT = dijkstra.getShostestDistance(graph.getTerminalNode());

			penalty = cost + Math.pow(dT/traversableProb,-Math.log(traversableProb));//This is the weight function for "distance-to-termination" (DT) algorithm. Please refer to Aksakalli and Ari, 2014 for more information.
            //penalty = cost + Math.exp(Math.log(2)*Math.log(10)) + (Math.exp(Math.log(2)*Math.log(10))*Math.log(2)*(dT - 5))/5 - 2*Math.exp(Math.log(2)*Math.log(10))*Math.log(20)*(traversableProb - 1/2) + (Math.exp(Math.log(2)*Math.log(10))*Math.pow((2 * traversableProb - 1), 2)*(Math.log(400) + 2*Math.pow(Math.log(20),2) + 4))/4 - Math.exp(Math.log(2)*Math.log(10))*((2*Math.log(2)*Math.log(20))/5 + 2/5)*(dT - 5)*(traversableProb - 1/2) + (Math.exp(Math.log(2)*Math.log(10))*Math.log(2)*Math.pow((dT - 5),2)*(Math.log(2) - 1))/50;
			//penalty = dT/e.getLength()*traversableProb;
			//penalty = dT/traversableProb;
			//penalty = dT/e.getLength()/traversableProb;

            //penalty = cost + (1 - traversableProb) * slope;
		}
		return penalty;
	}

    public double getSPSAPenalty(Graph graph, Edge e, double[] a) {
        double penalty = 0;
        double r = e.getProb();
        if (r == 0) {
            r = 0.00001;
        }
        else if (r < 0.9999) {
            DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
            if (e.getDestination().getId()+1==graph.getVertexes().size()){
                dijkstra.execute(e.getSource());
            }
            else {
                dijkstra.execute(e.getDestination());
            }
            double dT = dijkstra.getShostestDistance(graph.getTerminalNode());
            penalty = cost + a[0] + a[1]*(dT-500) + a[2]*(r-0.5) + a[3]*(dT-500)*(r-0.5) + a[4]*Math.pow((dT-500),2) + a[5]*Math.pow((r-0.5),2); //Bivariate quadratic function.
            if(e.getId()==20) {
                //System.out.println(penalty);
            }

            if (penalty < 0) {
                penalty = 0;
            }
        }
        return penalty;
    }
	
	//----------------------------------------------------------------------------------------------------------------

    public double runSimulation(String filename) throws IOException {
		nodes = new ArrayList<Vertex>();
		edges = new ArrayList<Edge>();
		double omt = 0.0;
        int performance = 0;
		int iterationNumber = 1; //times current graph will be iterated (beware: if you selected "obtainWeatherFromFile" to be true, you'll need .prb files for all of the corresponding seeds) default is 1. keep it 1.
		
		System.out.println("Solving "+filename+".graph...");
		System.out.println("-------------------------");

		CSVReader reader = new CSVReader(new FileReader("input"+File.separator+filename+".graph"), (char)32);
		String [] nextLine;
		nextLine = reader.readNext();
		int numberOfNodes = Integer.parseInt(nextLine[1]);
		int numberOfEdges = Integer.parseInt(nextLine[2]);
		for (int i = 1; i <= numberOfNodes; i++) {
			Vertex location = new Vertex(i-1, "Node_" + i);
			nodes.add(location);
		}	
		int i=1;
		while ((nextLine = reader.readNext()) != null) {
			//the edges in the input file are being added into the graph.
	    	addLane(i-1, Integer.parseInt(nextLine[1])-1, Integer.parseInt(nextLine[2])-1, Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 0);
	    	//opposite direction of the edges above are being added into the graph (the graph has to be undirected).
	    	addLane(i, Integer.parseInt(nextLine[2])-1, Integer.parseInt(nextLine[1])-1, Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 0);
	        i++;
	    }
	    reader.close();
	    
	    Graph graph = new Graph(nodes, edges, nodes.get(0), nodes.get(numberOfNodes-1));
	    
	    DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
		dijkstra.execute(nodes.get(0));
		System.out.println("DIJ: "+dijkstra.getShostestDistance(nodes.get(numberOfNodes-1))+" (perfect weather)");

		//simulation results is written into SimulationResuls.csv
		CSVWriter writer_results = new CSVWriter(new FileWriter("SimulationResults.csv",true));
		for(int counter=1; counter<=iterationNumber; counter++) {
            double t1 = System.currentTimeMillis();

			String[] outputRow = new String[6];
			outputRow[0]=Integer.toString(counter);//What is this?
			
			edgesForOptimistic = new ArrayList<Edge>();
			nodesForOptimistic = new ArrayList<Vertex>();
			edgesForOptimistic = cloneEdges(edges);
			nodesForOptimistic = cloneVertices(nodes);
			
			weightedEdges = new ArrayList<Edge>();
			weightedNodes = new ArrayList<Vertex>();
			weightedEdges = cloneEdges(edges);
			weightedNodes = cloneVertices(nodes);
			
			prunedEdges = new ArrayList<Edge>();
			prunedNodes = new ArrayList<Vertex>();
			prunedEdges = cloneEdges(edges);
			prunedNodes = cloneVertices(nodes);

			if(obtainWeatherFromFile = true) {
				//obtainWeatherFromFile = true; EDGE STATUS IS OBTAINED FROM PRB FILES
				CSVReader reader2 = new CSVReader(new FileReader("input/"+File.separator+filename+".prb"), (char)32);
				String [] nextLine2;
				i=0;
				while ((nextLine2 = reader2.readNext()) != null) {
					int s = Integer.parseInt(nextLine2[2]);
					edges.get(i).setActualStatus(s);
					edges.get(i+1).setActualStatus(s);
					edgesForOptimistic.get(i).setActualStatus(s);
					edgesForOptimistic.get(i+1).setActualStatus(s);
					weightedEdges.get(i).setActualStatus(s);
					weightedEdges.get(i+1).setActualStatus(s);
					weightedEdges.get(i).setWeight( edges.get(i).getLength() + getPenalty(graph,edges.get(i)) );
					weightedEdges.get(i+1).setWeight( weightedEdges.get(i).getWeight() );
					prunedEdges.get(i).setActualStatus(s);
					prunedEdges.get(i+1).setActualStatus(s);
					prunedEdges.get(i).setWeight(prunedEdges.get(i).getLength()+(999999999*(1-prunedEdges.get(i).getActualStatus())));
					prunedEdges.get(i+1).setWeight(prunedEdges.get(i).getLength()+(999999999*(1-prunedEdges.get(i).getActualStatus())));
					i=i+2;
			    } 
			    reader2.close();
			}  
			else {
			   //obtainWeatherFromFile = false; EDGE STATUS IS DETERMINED PROBABILISTICALLY
				for (i=0; i<2*numberOfEdges; i=i+2) {
					Binomial binomialDist = new Binomial(1,edges.get(i).getProb());
					int randomIntBinomial = (int) binomialDist.random();					
					edges.get(i).setActualStatus(randomIntBinomial);
					edges.get(i+1).setActualStatus(randomIntBinomial);
					edgesForOptimistic.get(i).setActualStatus(randomIntBinomial);
					edgesForOptimistic.get(i+1).setActualStatus(randomIntBinomial);
					weightedEdges.get(i).setActualStatus(randomIntBinomial);
					weightedEdges.get(i+1).setActualStatus(randomIntBinomial);
					weightedEdges.get(i).setWeight( edges.get(i).getLength() + getPenalty(graph,edges.get(i)));
					weightedEdges.get(i+1).setWeight(weightedEdges.get(i).getWeight());
					prunedEdges.get(i).setActualStatus(randomIntBinomial);
					prunedEdges.get(i+1).setActualStatus(randomIntBinomial);
					prunedEdges.get(i).setWeight(prunedEdges.get(i).getLength()+(999999999*(1-prunedEdges.get(i).getActualStatus())));
					prunedEdges.get(i+1).setWeight(prunedEdges.get(i).getLength()+(999999999*(1-prunedEdges.get(i).getActualStatus())));
				}
			}

            double t2 = System.currentTimeMillis();
            double initialize = (t2 - t1) / 1000;

			Graph graphForOptimal = new Graph(prunedNodes, prunedEdges, prunedNodes.get(0), prunedNodes.get(numberOfNodes-1));
		    DijkstraAlgorithm dijkstraForOptimal = new DijkstraAlgorithm(graphForOptimal);
			dijkstraForOptimal.execute(prunedNodes.get(0));
			double dijCurrent = dijkstraForOptimal.getShostestDistance(prunedNodes.get(numberOfNodes-1));
			System.out.println("DIJ: "+dijCurrent+" (current weather)");
            thePath = "1";
			
			//Optimistic approach (dijkstra with actual statuses)
            double omtStart = System.currentTimeMillis();
			Graph graphForOptimistic = new Graph(nodesForOptimistic, edgesForOptimistic, nodesForOptimistic.get(0), nodesForOptimistic.get(numberOfNodes-1));
			omt = computeExpectedLength(graphForOptimistic,graphForOptimistic.getStartNode());
            double omtEnd = System.currentTimeMillis();
            double omtTime = initialize + ((omtEnd - omtStart) / 1000);

			//If the expectedDistance is less than zero, then problem is unsolvable under the current actual statuses of the edges.
			if(expectedDistance > 0) {
				outputRow[1]=Double.toString(expectedDistance);
			} else {
				outputRow[1]="-";
			}
			System.out.println("-------------------------");
			System.out.println("OMT: "+omt+" ("+omtTime+")");
			thePath = thePath.replaceAll("Node_", "");
			System.out.println(thePath);
			System.out.println("-------------------------");
			thePath = "1";
			
			//PBA approach (assigning weights using the selected penalty function)
            double pbaStart = System.currentTimeMillis();
			Graph graphWeighted = new Graph(weightedNodes, weightedEdges, weightedNodes.get(0), weightedNodes.get(numberOfNodes-1));
			expectedDistance = computeExpectedLength(graphWeighted, graphWeighted.getStartNode());
			thePath = thePath.replaceAll("Node_", "");
            double pbaEnd = System.currentTimeMillis();
			double pbaTime = initialize + ((pbaEnd - pbaStart)/1000);
			
			if (expectedDistance-dijCurrent<=omt-dijCurrent) {
				performance = 1;
			}
			
			//If the expectedDistance is less than zero, then problem is unsolvable under the current actual statuses of the edges.
			if(expectedDistance >= 0) {
				outputRow[0] = String.valueOf(filename);
                outputRow[1] = String.valueOf(thePath);
				outputRow[2] = Double.toString(omtTime);
				outputRow[3] = Double.toString(pbaTime);
				outputRow[4] = Double.toString(omt);
				outputRow[5] = String.valueOf(expectedDistance);
				expected = expected + expectedDistance;
				total = total + 1;
			} else {
				outputRow[1]="-";
			}
			System.out.println("PBA: "+expectedDistance+" ("+pbaTime+")");
			System.out.println(thePath);
			writer_results.writeNext(outputRow);
	
			System.out.println(performance);
		}
		writer_results.close();
		/*Printing the path
		LinkedList<Vertex> path = dijkstra.getPath(nodes.get(numberOfNodes-1));		
		for (Vertex vertex : path) {
			System.out.println(vertex);
		}
		*/
        return omt - expectedDistance;
	}

	public double computeExpectedLength(Graph graph, Vertex source) {
		double trueWeight = 0, falseWeight=0;
		double expectedLength = 0;
				
		if(hasNeighbor(graph.getEdges(), source)) {
			DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
			dijkstra.execute(source);
			LinkedList<Vertex> path = dijkstra.getPath(graph.getTerminalNode());
			
			if (path!=null && path.size() > 2) {
				Edge currentEdge = graph.getEdge(path.get(0), path.get(1));
				if(currentEdge.getActualStatus()==1) {
					graph.getEdges().get(graph.getEdges().indexOf(currentEdge)).setProb(1.0);
					thePath = thePath + "-" + currentEdge.getDestination();
					trueWeight = computeExpectedLength(graph, path.get(1));
					expectedLength = currentEdge.getLength() + trueWeight;
				}
                else {
					for(Edge e: graph.getEdges()) {
						if( (e.getSource().equals(path.get(0)) && e.getDestination().equals(path.get(1))) 
								|| (e.getSource().equals(path.get(1)) && e.getDestination().equals(path.get(0))) ) {
							graph.getEdges().get(graph.getEdges().indexOf(e)).setWeight(Double.POSITIVE_INFINITY);
							graph.getEdges().get(graph.getEdges().indexOf(e)).setLength(Double.POSITIVE_INFINITY);
						}
					}
					if(hasNeighbor(graph.getEdges(), path.get(0))) {
						falseWeight = computeExpectedLength(graph, path.get(0));
						expectedLength = falseWeight;
					}
                    else {
						System.out.println("Problem is unsolvable under the current actual statuses (1)");
						expectedLength = -9999999;
					}
				}	
					
			}
            else {
				if (path!=null) {
					Edge currentEdge = graph.getEdge(path.get(0), path.get(1));
					if (currentEdge.getActualStatus()!=0) {
						expectedLength = currentEdge.getLength();

						thePath = thePath + "-" + currentEdge.getDestination();
					}
					else {
						for(Edge e: graph.getEdges()) {
							if( (e.getSource().equals(path.get(0)) && e.getDestination().equals(path.get(1))) 
									|| (e.getSource().equals(path.get(1)) && e.getDestination().equals(path.get(0))) ) {
								graph.getEdges().get(graph.getEdges().indexOf(e)).setWeight(Double.POSITIVE_INFINITY);
								graph.getEdges().get(graph.getEdges().indexOf(e)).setLength(Double.POSITIVE_INFINITY);
							}
						}
						trueWeight = computeExpectedLength(graph, path.get(0));
						expectedLength = trueWeight;
					}
				}
                else {
					System.out.println("The path is null !!!");
					expectedLength = -9999999;
				}
			}
		}
        else {
			System.out.println("Problem is unsolvable under the current actual statuses (2)");
			expectedLength = -9999999;
		}
		return expectedLength;
	}
	
	public List<Edge> cloneEdges(List<Edge> initialEdges) {
		List<Edge> newEdges = new ArrayList<Edge>();
		for(Edge e : initialEdges) {
			Edge newEdge = new Edge(e.getId(), e.getSource(), e.getDestination(), e.getLength(), e.getProb(), e.getActualStatus());
			newEdges.add(newEdge);
		}
		return newEdges;
	}
	
	public List<Vertex> cloneVertices(List<Vertex> initialVertices) {
		List<Vertex> newVertices = new ArrayList<Vertex>();
		for(Vertex v : initialVertices) {
			Vertex newVertex = new Vertex(v.getId(), v.getName());
			newVertices.add(newVertex);
		}
		return newVertices;
	}
	
	private boolean hasNeighbor(List<Edge> edges, Vertex node) {
		boolean returnValue=false;
		for (Edge edge : edges) {
			if (edge.getSource().equals(node)) {
				returnValue = true;
				break;
			}
		}
		return returnValue;
	}
	
	private void addLane(Integer laneId, int sourceLocNo, int destLocNo,
			double duration, double probability, int actualStatus) {
		Edge lane = new Edge(laneId,nodes.get(sourceLocNo), nodes.get(destLocNo), duration, probability, actualStatus);
		edges.add(lane);
	}
}