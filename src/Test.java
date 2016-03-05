import au.com.bytecode.opencsv.CSVReader;
import dijkstra.engine.DijkstraAlgorithm;
import dijkstra.model.Edge;
import dijkstra.model.Graph;
import dijkstra.model.Vertex;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;


public class Test {
	public static void runTest(String filename) throws IOException {
		List<Vertex> nodes = new ArrayList<Vertex>();
		List<Edge> edges = new ArrayList<Edge>();
		
		System.out.println("Initializing...");
		
		//Graph information (Edges,Vertices,transition probabilities and edge weights) is loading .....
		//Input file is opening...
		CSVReader reader = new CSVReader(new FileReader("input"+File.separator+filename+".graph"), (char)32);
		String [] nextLine;
		nextLine = reader.readNext();
		int numberOfNodes = Integer.parseInt(nextLine[1]);
		for (int i = 1; i <= numberOfNodes; i++) {
			Vertex location = new Vertex(i-1, "Node_" + i);
			nodes.add(location);
		}
		
		int i=1;
		while ((nextLine = reader.readNext()) != null) {
			//The edge in the input file is being added into graph . . .
			
			Edge lane = new Edge(i-1, nodes.get(Integer.parseInt(nextLine[1])-1), nodes.get(Integer.parseInt(nextLine[2])-1), Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 1);
	    	edges.add(lane);
	    	//Opposite direction of the edge above is being added into graph . . . (Because, the graph has to be undirected)
	    	Edge laneOpposite = new Edge(i, nodes.get(Integer.parseInt(nextLine[2])-1), nodes.get(Integer.parseInt(nextLine[1])-1), Integer.parseInt(nextLine[4]), Double.parseDouble(nextLine[3]), 1);
	    	edges.add(laneOpposite);
	        i++;
	    } 
	    reader.close();
	    
	    CSVReader reader2 = new CSVReader(new FileReader("input/"+File.separator+filename+".prb"), (char)32);
		String [] nextLine2;
		i=0;
		while ((nextLine2 = reader2.readNext()) != null) {
			int s = Integer.parseInt(nextLine2[2]);
			//System.out.println(s);
			edges.get(i).setActualStatus(s);
			edges.get(i+1).setActualStatus(s);
			i=i+2;
	    } 
	    reader2.close();
	    
	    Graph graph = new Graph(nodes, edges, nodes.get(0), nodes.get(numberOfNodes-1));
	    
	    DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
		dijkstra.execute(nodes.get(0));
		double deterministicLength = dijkstra.getShostestDistance(nodes.get(numberOfNodes-1));
		//System.out.println("Deterministic Shortest Path: "+dijkstra.getShostestDistance(nodes.get(numberOfNodes-1)));
			
		//Optimistic approach
		double optimisticNDR = computeExpectedLength(graph,graph.getStartNode());
		//If the expectedDistance is less than zero, then problem is unsolvable under the current actual statuses of the edges.
		if(deterministicLength == optimisticNDR) {
			System.out.println("Test is OK! ("+deterministicLength+" = "+optimisticNDR+")");
		} else {
			System.out.println("("+deterministicLength+" != "+optimisticNDR+")");
			//System.exit(0);
		}
		System.out.println("Done.");
	}
	
	public static double computeExpectedLength(Graph graph, Vertex source) {
		double trueWeight = 0, falseWeight=0;//9999999 is just a big number.
		double expectedDistance = 0;//-9999999 is just a big negative number, and it makes the ecpectedDistance negative when problem is unsolvable.
		
		if(hasNeighbor(graph.getEdges(), source)) {
			
			DijkstraAlgorithm dijkstra = new DijkstraAlgorithm(graph);
			dijkstra.execute(source);
			LinkedList<Vertex> path = dijkstra.getPath(graph.getTerminalNode());
			//dijkstra.getShostestDistance(nodes.get(target.getId()));;
			
			if (path!=null && path.size()>2) {
				Edge currentEdge = graph.getEdge(path.get(0), path.get(1));
				
				if(currentEdge.getActualStatus()==1) {
					graph.getEdges().get(graph.getEdges().indexOf(currentEdge)).setProb(1.0);
					trueWeight = computeExpectedLength(graph, path.get(1));
					expectedDistance = currentEdge.getLength() + trueWeight;
				} else {
					//graph.getEdges().remove(graph.getEdges().indexOf(currentEdge));
					
					for(Edge e: graph.getEdges()) {
						if( (e.getSource().equals(path.get(0)) && e.getDestination().equals(path.get(1))) 
								|| (e.getSource().equals(path.get(1)) && e.getDestination().equals(path.get(0))) ) {
							graph.getEdges().get(graph.getEdges().indexOf(e)).setWeight(Double.POSITIVE_INFINITY);
							graph.getEdges().get(graph.getEdges().indexOf(e)).setLength(Double.POSITIVE_INFINITY);
						}
					}

					
					if(hasNeighbor(graph.getEdges(), path.get(0))) {
						falseWeight = computeExpectedLength(graph, path.get(0));
						expectedDistance = falseWeight;
					} else {
						System.out.println("Problem is unsolvable!");
						//System.exit(0);
						expectedDistance = -9999999;
					}
				}	
					
				
			} else {
				if (path!=null) {
					Edge currentEdge = graph.getEdge(path.get(0), path.get(1));
					if (currentEdge.getActualStatus()!=0) {
						expectedDistance = currentEdge.getLength();
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
						expectedDistance = trueWeight;
					}
				} else {
					System.out.println("No feasible paths under the current set of actual statuses!");
					expectedDistance = -9999999;
				}
			}

		} else {
			System.out.println("Problem is unsolvable!");
			//System.exit(0);
			expectedDistance = -9999999;
		}
		
		return expectedDistance;

	}
	
	private static boolean hasNeighbor(List<Edge> edges, Vertex node) {
		boolean returnValue=false;
		for (Edge edge : edges) {
			if (edge.getSource().equals(node)) {
				returnValue = true;
				break;
			}
		}
		return returnValue;
	}
}
