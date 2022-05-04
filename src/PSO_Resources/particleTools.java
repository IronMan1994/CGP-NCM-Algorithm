package PSO_Resources;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import GA_Resources.Chromosome;

@SuppressWarnings("all")
public class particleTools<T extends Comparable>
{
	
	//Rank the population according to the fitness
	public static void sort_pop_list(List<Particle> pop)
	{
		Collections.sort(pop);
	}
	
	
	//particle to particle
	public static Particle copy_particle(Particle a){
		Integer[] Xi = new Integer[a.Xi.length];
		Double[] Vi = new Double[a.Vi.length];
		Integer[] chromosome = new Integer[a.chromosome.length];
		Double[] fitness = new Double[a.fitness.length];
		Double[] s = new Double[a.s.length];
		
		copy_int(a.Xi, Xi);
		copy_int(a.Vi, Vi);
		copy_int(a.chromosome, chromosome);
		copy_int(a.fitness, fitness);
		copy_int(a.s, s);
		
		Particle temp = new Particle(Xi, Vi, chromosome, fitness, a.index, s);
		
		return temp;
	}
	
	
	//chromosome to particle
	public static Particle copy_individual_To_particle(Chromosome a, int n, double Vm){
		Random random = new Random();
		Integer[] Xi = new Integer[n];
		Double[] Vi = new Double[n];
		Integer[] chromosome = new Integer[a.chromosome.length];
		Double[] fitness = new Double[a.fitness.length];
		Double[] s = new Double[n];
		
		for(int i = 0; i < Xi.length; i++){
			Xi[i] = 0;
			Vi[i] = random.nextInt(65535) % (Vm+1);
		}
		
		for(int i = 0; i < a.chromosome.length; i++){
			Xi[a.chromosome[i]] = 1;
			Vi[a.chromosome[i]] = Math.random() * 2;
			chromosome[i] = a.chromosome[i];
		}
		
		for(int i = 0; i < n; i++){
			s[i] = Math.abs(1 - 2 / (1 + Math.pow(Math.E, -Vi[i])));
		}

		for(int i = 0; i < a.fitness.length; i++){
			fitness[i] = a.fitness[i];
		}
		
		Particle temp = new Particle(Xi, Vi, chromosome, fitness, -1, s);
		return temp;
	}
	
	
	//Copy particle populations
	public static void copy_list(List<Particle> a, List<Particle> b){
		b.clear();
		for(int i = 0; i < a.size(); i++){
			Integer[] Xi = new Integer[a.get(i).Xi.length];
			Double[] Vi = new Double[a.get(i).Vi.length];
			Integer[] chromosome = new Integer[a.get(i).chromosome.length];
			Double[] fitness = new Double[a.get(i).fitness.length];
			Double[] s = new Double[a.get(i).s.length];
			
			copy_int(a.get(i).Xi, Xi);
			copy_int(a.get(i).Vi, Vi);
			copy_int(a.get(i).chromosome, chromosome);
			copy_int(a.get(i).fitness, fitness);
			copy_int(a.get(i).s, s);
			
			Particle temp = new Particle(Xi, Vi, chromosome, fitness, a.get(i).index, s);
			b.add(temp);
		}
	}
	
	
	//Copy array
	public static <T> void copy_int(T[] a, T[] b){
		for(int i = 0; i < a.length; i++){
			b[i] = a[i];
		}
	}
	
	
	//sort(Bubble)
	public static void sort_list_maxgene(List<Particle_Index> maxGene)
	{
		Particle_Index temp = new Particle_Index();
        int size = maxGene.size();
        for(int i = 0 ; i < size-1; i ++){
	        for(int j = 0 ;j < size-1-i ; j++){
	        	//<表示把最小值移到最后
	            if(maxGene.get(j).speciesIndividual.fitness[0] < maxGene.get(j+1).speciesIndividual.fitness[0])  //交换两数位置
	            {
		            temp = maxGene.get(j);
		            maxGene.set(j, maxGene.get(j+1));
		            maxGene.set(j+1, temp);
	            }
	        }
        }
	}
	
	
	//Calculating algorithm accuracy
	public static double calAccuracy(double standard, List<Particle_Index> maxGene){
		
		double result = 0.0;
		
		for (Particle_Index item : maxGene){
			result += item.speciesIndividual.fitness[0] / standard;
		}
		
		result /= maxGene.size();
		BigDecimal bd = new BigDecimal(result*100).setScale(2, RoundingMode.UP);
		return bd.doubleValue();
	}

	
}
