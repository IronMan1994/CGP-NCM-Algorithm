package PSO_Resources;

public class Particle implements Comparable
{
	//二进制编码
	//Xi
	public Integer[] Xi;
	//Vi
	public Double[] Vi;
	//几个对应的适应度(第一个为总和的适应度)
	public Double[] fitness;
	//对应的Xi中为1的下标，方便找基因名
	public Integer[] chromosome;
	
	//粒子编号
	public int index;
	
	public Double[] s;
	

	public Particle(Integer[] Xi, Double[] Vi, Integer[] chromosome, Double[] fitness, int index, Double[] s){
		this.Xi = Xi;
		this.Vi = Vi;
		this.chromosome = chromosome;
		this.fitness = fitness;	
		this.index = index;
		this.s = s;
	}
	
	public Particle(){
		
	}

	@Override
	public String toString(){
		String str = "";
		//位置数组
		for(int i = 0; i < Xi.length; i++){
			str += Xi[i] + "\t";
		}
		str += "\n";
		
		//速度数组
		for(int i = 0; i < Vi.length; i++){
			str += Vi[i] + "\t";
		}
		str += "\n";
		
		//对应的基因序号和适应值
		for(int i = 0; i < chromosome.length; i++){
			str += chromosome[i] + "\t";
		}
		for(int i = 0; i < fitness.length; i++){
			str += fitness[i] + "\t";
		}
		str += "\n";
		return str;
	}

	@Override
	public int compareTo(Object o)
	{
		double fitness0 = ((Particle) o).fitness[0];
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
