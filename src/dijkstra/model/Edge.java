package dijkstra.model;

public class Edge  {
	private final Integer id; 
	private final Vertex source;
	private final Vertex destination;
	private double length; 
	private double weight; 
	private double probability;
	private int actualStatus;
	
	public Edge(Integer id, Vertex source, Vertex destination, double length, double probability, int actualStatus) {
		this.id = id;
		this.source = source;
		this.destination = destination;
		this.length = this.weight = length;
		this.probability = probability;
		this.actualStatus = actualStatus;
	}
	
	public Integer getId() {
		return id;
	}
	public Vertex getDestination() {
		return destination;
	}

	public Vertex getSource() {
		return source;
	}
	public double getWeight() {
		return weight;
	}
	public double getLength() {
		return length;
	}
	
	public void setWeight(double weight) {
		this.weight = weight;
	}
	public void setLength(double length) {
		this.length = length;
	}
	
	public double getProb() {
		return probability;
	}
	
	public void setProb(double probability) {
		this.probability = probability;
	}
	
	public int getActualStatus() {
		return actualStatus;
	}
	
	public void setActualStatus(int actualStatus) {
		this.actualStatus = actualStatus;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + actualStatus;
		result = prime * result
				+ ((destination == null) ? 0 : destination.hashCode());
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		long temp;
		temp = Double.doubleToLongBits(length);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(probability);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((source == null) ? 0 : source.hashCode());
		temp = Double.doubleToLongBits(weight);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Edge other = (Edge) obj;
		if (actualStatus != other.actualStatus)
			return false;
		if (destination == null) {
			if (other.destination != null)
				return false;
		} else if (!destination.equals(other.destination))
			return false;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		if (Double.doubleToLongBits(length) != Double
				.doubleToLongBits(other.length))
			return false;
		if (Double.doubleToLongBits(probability) != Double
				.doubleToLongBits(other.probability))
			return false;
		if (source == null) {
			if (other.source != null)
				return false;
		} else if (!source.equals(other.source))
			return false;
		if (Double.doubleToLongBits(weight) != Double
				.doubleToLongBits(other.weight))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return source + " " + destination;
	}
	
	
}
