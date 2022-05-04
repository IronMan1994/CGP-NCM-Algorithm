package GA_Resources;

public class Chromosome implements Comparable
{
	//一条染色体
	public int[] chromosome;

	

	//几个对应的适应度(第一个为总和的适应度)
	public double[] fitness;
	
	public Chromosome(int[] chromosome, double[] fitness)
	{
		this.chromosome = chromosome;
		this.fitness = fitness;
	}
	
	public Chromosome()
	{
		
	}

	public double getFitness(){
		return fitness[0];
	}
	
	@Override
	public String toString()
	{
		String str = "";
		for(int i = 0; i < chromosome.length; i++)
		{
			str += chromosome[i] + "\t";
		}
		for(int i = 0; i < fitness.length; i++)
		{
			str += fitness[i] + "\t";
		}
		return str;
	}

	@Override
	public int compareTo(Object o)
	{
		double fitness0 = ((Chromosome) o).fitness[0];
		if(this.fitness[0] > fitness0){
			return -1;
		}
		if(this.fitness[0] == fitness0){
			return 0;
		}
		else {
			return 1;
		}
	}	
}
