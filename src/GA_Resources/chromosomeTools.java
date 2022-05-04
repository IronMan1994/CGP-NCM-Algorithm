package GA_Resources;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Collections;
import java.util.List;

import PSO_Resources.Particle;

@SuppressWarnings("all")
public class chromosomeTools
{
	//Sort pop
	public static void sort_pop_list(List<Chromosome> pop)
	{
		Collections.sort(pop);
	}
	
	//Calculating roulette policy probability
	public static double[] crossover_rate(List<Chromosome> Cross_pool)
	{
		int geneSize = Cross_pool.size();
		int fuzhiNum = 0;
		double[] rate_array = new double[geneSize - fuzhiNum + 1];
		double rate = 0.0;
		
		double sum = sum_fitness(fuzhiNum, geneSize, Cross_pool);
		rate_array[0] = 0.0;
		
		for(int i = 1; i < rate_array.length; i++){
			rate += Cross_pool.get(fuzhiNum - 1 + i).fitness[0] / sum;
			rate_array[i] = rate;
		}
		return rate_array;
	}
	
	//Calculate the sum of fitness for all individuals in a population
	public static double sum_fitness(int fuzhiNum, int size, List<Chromosome> pop)
	{
		double count = 0.0;
		for(int i = fuzhiNum; i < size; i++)
		{
			count += pop.get(i).fitness[0];
		}
		return count;
	}
	
	//copy pop
	public static void copy_list(List<Chromosome> a, List<Chromosome> b)
	{
		b.clear();
		for(int i = 0; i < a.size(); i++){
			Chromosome temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}
	
	//copy chromosome
	public static Chromosome copy_SpeciesIndividual(Chromosome from)
	{
		int[] chromosome = new int[from.chromosome.length];
		double[] fitness = new double[from.fitness.length];
		
		for(int i = 0; i < chromosome.length; i++){
			chromosome[i] = from.chromosome[i];
		}
		
		for(int i = 0; i < fitness.length; i++){
			fitness[i] = from.fitness[i];
		}
		Chromosome result = new Chromosome(chromosome, fitness);

		return result;
	}
	
	
	public static void printPOP(List<Chromosome> tempPOP)
	{
		for(int i = 0; i < tempPOP.size(); i++){
			for(int j = 0; j < tempPOP.get(i).chromosome.length; j++){			
				System.out.print(String.format("%-12s", tempPOP.get(i).chromosome[j]));
			}
			for(int j = 0; j < tempPOP.get(i).fitness.length; j++){
				System.out.print(String.format("%-12s", tempPOP.get(i).fitness[j]));
			}
			System.out.println();
		}
	}
	
	
	//The results of each run of the algorithm, sorted by fitness from high to low (bubble sort)
	public static void sort_list_maxgene(List<Chromosome_Index> temppop)
	{
		Chromosome_Index temp = new Chromosome_Index();
        int size = temppop.size();
        for(int i = 0 ; i < size-1; i ++){
	        for(int j = 0 ;j < size-1-i ; j++){
	            if(temppop.get(j).speciesIndividual.fitness[0] < temppop.get(j+1).speciesIndividual.fitness[0])  //交换两数位置
	            {
		            temp = temppop.get(j);
		            temppop.set(j, temppop.get(j+1));
		            temppop.set(j+1, temp);
	            }
	        }
        }
	}
	
}
