package dijkstra.model;

import java.util.List;

public class Graph {
	private final List<Vertex> vertexes;
	private final List<Edge> edges;
	private final Vertex startNode;
	private final Vertex terminalNode;

	public Graph(List<Vertex> vertexes, List<Edge> edges, Vertex startVertex, Vertex terminalVertex) {
		this.vertexes = vertexes;
		this.edges = edges;
		this.startNode = startVertex;
		this.terminalNode = terminalVertex;
	}

	public List<Vertex> getVertexes() {
		return vertexes;
	}

	public List<Edge> getEdges() {
		return edges;
	}

	public Vertex getStartNode() {
		return startNode;
	}

	public Vertex getTerminalNode() {
		return terminalNode;
	}
	
	public Edge getEdge(Vertex source, Vertex end) {
		Edge edge = null;
		for (Edge e : edges) {
			if (e.getSource().equals(source) && e.getDestination().equals(end)) {
				edge =  e;
				break;
			}
		}
		return edge;
	}
	
}