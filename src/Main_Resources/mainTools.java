package Main_Resources;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import GA_Resources.Chromosome;
import GA_Resources.Chromosome_Index;
import GA_Resources.GA_Main;
import GA_Resources.chromosomeTools;
import PSO_Resources.PSO_Main;
import PSO_Resources.Particle;

@SuppressWarnings("all")
public class mainTools
{
	//record and output to file
	public static void writeResult(List temp, String initStr) throws IOException
	{
		File f = new File("");
		String pathname = f.getCanonicalPath();
		File ff = new File(pathname + "/src/results/Results.txt");
		OutputStream os = new FileOutputStream(ff, true);
		try{
			os.write((initStr + "：" ).getBytes());
			os.write("\r\n".getBytes());
			for(int i = 0; i < temp.size(); i++){
				os.write(String.format("%-12s",  " ").getBytes());
				os.write((String.valueOf(temp.get(i))).getBytes());
				os.write("\r\n".getBytes());
			}
			os.write("\r\n\r\n\r\n".getBytes());
		}
		catch (Exception e){
			e.printStackTrace();
		}
		finally{
			os.close();
		}
	}
	
	
	//sort(Bubble)
	public static void sort_list_maxgene(List<Chromosome_Index> maxGene)
	{
		Chromosome_Index temp = new Chromosome_Index();
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
	
	
	//Particle to Chromosome
	public static Chromosome copy_SpeciesIndividual(Particle from){
		int[] chromosome = new int[from.chromosome.length];
		double[] fitness = new double[from.fitness.length];
		
		for(int i = 0; i < chromosome.length; i++){
			chromosome[i] = from.chromosome[i];
		}
		
		for(int i = 0; i < fitness.length; i++){
			if(from.fitness[i] != null){
				fitness[i] = from.fitness[i];
			}
		}
		
		Chromosome result = new Chromosome(chromosome, fitness);
		
		return result;
	}
	
	
	//copy Chromosome
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
	
	
	//append list to list
	public static void copy_list(List<Chromosome> a, List<Chromosome> b)
	{
		for(int i = 0; i < a.size(); i++){
			Chromosome temp = copy_SpeciesIndividual(a.get(i));
			b.add(temp);
		}
	}

	
	//copy list
	public static void copy_list2(List<Chromosome> from, List<Chromosome> to)
	{
		to.clear();
		for(int i = 0; i < from.size(); i++){
			Chromosome temp = copy_SpeciesIndividual(from.get(i));
			to.add(temp);
		}
	}
	
	
	//append list of Particle to list of Chromosome
	public static void copy_list3(List<Particle> a, List<Chromosome> b)
	{
		for(int i = 0; i < a.size(); i++){
			int length = a.get(i).chromosome.length;
			int[] chromosome = new int[length];
			for(int j = 0; j < length; j++){
				chromosome[j] = a.get(i).chromosome[j];
			}
			int fLength = a.get(i).fitness.length;
			double[] fitness = new double[fLength];
			for(int j = 0; j < fLength-1; j++){
				fitness[j] = a.get(i).fitness[j];
			}
			fitness[fLength-1] = 0;
			Chromosome speciesIndividual = new Chromosome(chromosome, fitness);
			b.add(speciesIndividual);
		}
	}
	
	
	//Rank the population according to the fitness
	public static void sort_pop_list(List<Chromosome> pop)
	{
		Collections.sort(pop);
	}
	
	
	//update cooperative pool
	public static void update_Cross_pool(GA_Main t, int size, List<Chromosome> Cross_pool)
	{
		double[] rate_array = chromosomeTools.crossover_rate(Cross_pool);
		List<Chromosome> tempList = new ArrayList<Chromosome>();
	
		for(int i = 0; i < size; i++){
			double r1 = Math.random();
			for(int j = 0; j < rate_array.length - 1; j++)
			{
				if(r1 >= rate_array[j] && r1 < rate_array[j + 1])
				{
					Chromosome temp = copy_SpeciesIndividual(Cross_pool.get(j));
					tempList.add(temp);
				}
			}
		}
		copy_list2(tempList, Cross_pool);	
	}

	
	//Cooperative operator [second]
	public static void share(double index, GA_Main t, GA_Main t2, PSO_Main t3){
		if(index == 3){
			t.pop.add(copy_SpeciesIndividual(t3.Pg));
			t2.pop.add(copy_SpeciesIndividual(t3.Pg));	
		}
	}
	
	
	//Another way to jump out is to end the iteration when all groups reach the maxt
	public static double[] bestFitnessAndIndex123(GA_Main t, GA_Main t2, PSO_Main t3, int[] breakArray, double[] preFitnessArray){
		
		double[] results = new double[2];
		
		double pre_fitness = t.pop.get(0).fitness[0] > t2.pop.get(0).fitness[0] ? t.pop.get(0).fitness[0] : t2.pop.get(0).fitness[0];
		int index123 = t.pop.get(0).fitness[0] > t2.pop.get(0).fitness[0] ? 1:2;
		pre_fitness = pre_fitness > t3.Pg.fitness[0] ? pre_fitness : t3.Pg.fitness[0];
		index123 = pre_fitness > t3.Pg.fitness[0] ? index123 : 3;
		
		if(t.pop.get(0).fitness[0] > preFitnessArray[0]){
			preFitnessArray[0] = t.pop.get(0).fitness[0];
			breakArray[0] = 0;
		}else{
			breakArray[0]++;
		}
		if(t2.pop.get(0).fitness[0] > preFitnessArray[1]){
			preFitnessArray[1] = t2.pop.get(0).fitness[0];
			breakArray[1] = 0;
		}else{
			breakArray[1]++;
		}
		
		if(t3.Pg.fitness[0] > preFitnessArray[2]){
			preFitnessArray[2] = t3.Pg.fitness[0];
			breakArray[2] = 0;
		}else{
			breakArray[2]++;
		}
		
		results[0] = pre_fitness;
		results[1] = index123;
		
		return results;
	}
	
	
	//Another way to jump out is to end the iteration when all groups reach the 'maxt'
	public static int ifBreak(int[] breakArray){
		int result = Math.min(breakArray[0], breakArray[1]);
		result = Math.min(result, breakArray[2]);
		return result;
	}
	
	
	//Calculating algorithm accuracy
	public static double calAccuracy(double standard, List<Chromosome_Index> maxgene){
		
		double result = 0.0;
		
		for (Chromosome_Index item : maxgene){
			result += item.speciesIndividual.fitness[0] / standard;
		}
		
		result /= maxgene.size();
		BigDecimal bd = new BigDecimal(result*100).setScale(2, RoundingMode.UP);
		return bd.doubleValue();
	}
	
}
