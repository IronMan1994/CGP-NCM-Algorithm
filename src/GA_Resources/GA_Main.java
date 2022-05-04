package GA_Resources;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.lang.reflect.Method;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class GA_Main
{
	//List of File
	public List<FileToList> A;
	//List of network
	public ArrayList<HashMap<String, ArrayList<String>>> netList;
	//gene name
	public String[] name;
	public int[] name_index;
	//Random()
	public Random random;
	//Model name(Adapting using a reflection mechanism)
	public Method method;
	
	//Current Chromosome Population
	public List<Chromosome> pop;
	//Next Chromosome Population
	public List<Chromosome> nextPOP;
	
	//Size of Chromosome Population
	public int sizeOfChromosome;
	//Number of gene
	public int G;
	//K
	public int K;	

	//Mutation probability
	public double pm;

	//Record time
	public long start;
	public long end;
	

	
	
	//initialization
	public void initData(String[] paths, String[] netPath, int G, int sizeOfChromosome, int K, 
			double pmChromosome, String modelName) throws Exception
	{
		this.G = G;
		this.K = K;
		this.sizeOfChromosome = sizeOfChromosome;
		this.pm = pmChromosome;
		this.name = new String[this.G];
		this.name_index = new int[this.G];
		this.pop = new ArrayList<Chromosome>();
		this.nextPOP = new ArrayList<Chromosome>();
		this.random = new Random(System.currentTimeMillis());
		this.A = new ArrayList<FileToList>();
		
		//File to list
		File f1 = new File("");
		for(int i = 0; i < paths.length;i++)
		{
			File f = new File(this.getClass().getResource("/" + paths[i]).getPath());
			f1 = f;
			List<float[]> B = new ArrayList<float[]>();
			readData(f, B);
			FileToList a = new FileToList(B);
			A.add(a);
		}
		
		//Loading gene name from file
		Reader reader = new FileReader(f1);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		try{
			while((temp = bfr.readLine()) != null){
				String[] ppp = temp.split("\t");
				if(ppp.length == this.G + 1){
					for(int i = 0; i < this.G; i++){
						name[i] = ppp[i + 1];
					}
					break;
				}
			}
			for(int i = 0; i < this.G; i++)
			{
				name_index[i] = i;
			}
		} catch (Exception e){
			e.printStackTrace();
		}finally {
			bfr.close();
			reader.close();
		}


		//Using reflection mechanism to find the required model
		Method[] methods = this.getClass().getMethods();
		for(Method m : methods)
		{
			if(modelName.equals(m.getName()))
			{
				this.method = m;
				break;
			}
		}
		
		//network to List
		netList = new ArrayList<HashMap<String, ArrayList<String>>>();
		for(int j = 0; j < paths.length; j++){
			recordNet(netPath[j]);
		}
	}	
	

	public void readData(File f, List<float[]> B) throws Exception
	{
		InputStream fis = new FileInputStream(f);
		Reader isr = new InputStreamReader(fis);
		BufferedReader bfr = new BufferedReader(isr);
		String tempstr;
		int line = 1;
		try
		{
			while((tempstr = bfr.readLine()) != null){
				String[] ppp = tempstr.split("\t");
				if(ppp.length == this.G + 1){	
					if(line == 2){
						float[] temp = new float[this.G];
						for(int i = 0; i < this.G; i++){
							//提取数字		
							temp[i] = Float.valueOf(ppp[i + 1].trim());		
						}
						B.add(temp);
					}						
				}
				line = 2;
			}		
		} 
		catch (Exception e){
			e.printStackTrace();
		}finally {
			bfr.close();
			fis.close();
		}
	}
	
	//loading network to list
	public void recordNet(String fileName) throws IOException{
		File f = new File(this.getClass().getResource("/" + fileName).getPath());
		HashMap<String, ArrayList<String>> tempMap = new HashMap<String, ArrayList<String>>();
		Reader reader = new FileReader(f);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		try
		{
			while((temp = bfr.readLine()) != null){
				String[] ppp = temp.split("\t");
				if(ppp.length > 1){
					ArrayList<String> tempArray = new ArrayList<String>();
					for (int i = 1; i < ppp.length; i++){
						tempArray.add(ppp[i]);
					}
					tempMap.put(ppp[0], tempArray);
				}
			}
			netList.add(tempMap);
		} catch (Exception e){
			e.printStackTrace();
		}finally{
			bfr.close();
			reader.close();
		}
	}
	
	
	//Counting the number of edges between pairs of K genes
	public int calEdgeCount(int[] chromosome){
		int result = 0;
		for(int i = 0; i < netList.size(); i++){
			int count = 0;
			HashMap<String, ArrayList<String>> tempMap = netList.get(i);
			for(int j = 0; j < chromosome.length; j++){
				if(tempMap.containsKey(name[chromosome[j]])){
					ArrayList<String> edges = tempMap.get(name[chromosome[j]]);
					for(int k = 0; k < chromosome.length; k++){
						if(j != k){
							if(edges.contains(name[chromosome[k]])){
								count++;
							}
						}
					}
				}
			}
			if(result < count){
				result = count;
			}
		}
		return result;
	}
	
	
	//Calculating Fitness(moudle of NCM)
	public double[] model_NCM(int[] chromosome, FileToList B)
	{	
		double[] result = new double[3];
		double fitness = 0.0;
		double red = 0.0,ifcancer = 0.0;
		int x = 0;
		int y = 0;
		List<float[]> sample_one = B.A;	
		for(int i = 0; i < sample_one.size();i++)
		{
			float[] temp = B.A.get(i);
			double[] temp_mean = new double[chromosome.length];
			double sum_line = 0.0; 
			int sum_line_half = 0;
			double max = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] > max)
				{
					max = temp[chromosome[j]];
				}	
				if(temp[chromosome[j]] == 0){
					sum_line_half++;
				}
				if(temp[chromosome[j]] >= 1){
					y++;
				}
				sum_line += temp[chromosome[j]];
				temp_mean[j] = temp[chromosome[j]];
			}
			if(max >= 1){
				x++;
			}
			if(sum_line != 0){
				double mean = sum_line / chromosome.length;
				double fangcha_fenzi = 0.0;
				for(int j = 0; j < temp_mean.length; j++){
					fangcha_fenzi += Math.pow(temp_mean[j] - mean, 2);
				}
				double cov = Math.sqrt(fangcha_fenzi / (this.K-1)) / mean;
				
				if(max <= 0.5){
					cov = cov / (Math.sqrt(this.K) * 2);
				}
				
				cov /= (Math.sqrt(this.K)/1.5);
				
				red += cov;
			}

			ifcancer += max;	
		}
		
		BigDecimal bd = new BigDecimal(ifcancer).setScale(4,RoundingMode.UP);
		ifcancer = bd.doubleValue();
		BigDecimal bd2 = new BigDecimal(red).setScale(4,RoundingMode.UP);
		red = bd2.doubleValue();
		fitness = (ifcancer + red);
		
		//Calculate the number of edges
		double edgeCount = calEdgeCount(chromosome);
		int maxEdgeCount = chromosome.length * (chromosome.length - 1) / 2;
		
		//Calculate value of CM
		BigDecimal bd3 = new BigDecimal(2 / (1.0/ifcancer + 1.0/red)).setScale(4, RoundingMode.UP);
		fitness = bd3.doubleValue();
		result[2] = fitness;
		
		//Calculate value of NCM
		if(edgeCount != 0){
			edgeCount = 1 +  (edgeCount)/(maxEdgeCount) ;
			BigDecimal eDecimal = new BigDecimal(edgeCount).setScale(4, RoundingMode.UP);
			edgeCount = eDecimal.doubleValue();
			BigDecimal bd4 = new BigDecimal(fitness * edgeCount).setScale(4, RoundingMode.UP);
			fitness = bd4.doubleValue();
		}
		
		result[0] = fitness;
		result[1] = edgeCount;
		
		return result;	
	}
	
	
	//Calculating Fitness(model of GA)
	public double[] model_GA(int[] chromosome, FileToList B)
	{
		double[] result = new double[3];
		double fitness = 0.0;
		double red = 0.0,ifcancer = 0.0, overlap = 0.0;
		List<float[]> sample_one = B.A;	
		for(int i = 0; i < sample_one.size();i++)
		{
			float[] temp = B.A.get(i);
			int tempifcancer = 0;
			for(int j = 0; j < chromosome.length;j++)
			{
				if(temp[chromosome[j]] != 0)
				{
					red += temp[chromosome[j]];
					tempifcancer = 1;
				}
			}
			ifcancer += tempifcancer;		
		}
		
		fitness = 2 * ifcancer - red;
		
		BigDecimal bd = new BigDecimal(fitness).setScale(4, RoundingMode.UP);
		
		overlap = red - ifcancer;

		result[0] = bd.doubleValue();
		result[1] = Double.valueOf(ifcancer);
		result[2] = Double.valueOf(red);
		
		return result;
	}
	
	
	//Construct initial chromosome population(Random Initialization)
	public void createBeginningSpeciesRandom() throws Exception
	{
		pop.clear();
		int count = 0;
		List<Integer> index = new ArrayList<Integer>();
		for(int i = 0; i < name_index.length; i++){
			index.add(i);
		}
		while(count < this.sizeOfChromosome){
			Collections.shuffle(index);
			int[] chromosome = new int[K];
			for(int m = 0; m < K; m++){
				chromosome[m] = index.get(m);
			}
			
			double[] fitness_result;
			Object[] args = {chromosome, this.A.get(0)};
			fitness_result = (double[]) method.invoke(this, args);
			Chromosome speciesIndividual = new Chromosome(chromosome, fitness_result);
			pop.add(speciesIndividual);
			count++;
		}//while		
	}
	
	
	//Select operation
	public void select()
	{
		nextPOP.clear();
		chromosomeTools.sort_pop_list(pop);
		for(int i = 0; i < sizeOfChromosome; i++){
			nextPOP.add(pop.get(i));
		}
		chromosomeTools.copy_list(nextPOP, pop);
	}

	
	//Crossover operation
	public void crossover(List<Chromosome> Cross_pool, double[] rate_array) throws Exception
	{
		double r1 = Math.random();
		int index1 = 0;
		
		Chromosome parent1 = null;
		for(int j = 0; j < rate_array.length - 1; j++){
			if(r1 >= rate_array[j] && r1 < rate_array[j + 1]){
				parent1 = Cross_pool.get(j);
				index1 = j;
			}
		}
		
		Chromosome parent2 = null;
		int index2 = index1;
		while(index1 == index2){
			double r2 = Math.random();
			for(int j = 0; j < rate_array.length - 1; j++){
				if(r2 >= rate_array[j] && r2 < rate_array[j + 1]){
					parent2 = Cross_pool.get(j);
					index2 = j;
				}
			}
		}
		
		//Start the crossover operation
		int[] chromosome1 = new int[this.K];
		int[] chromosome2 = new int[this.K];
		List<Integer> remainList1 = new ArrayList<Integer>();
		List<Integer> remainList2 = new ArrayList<Integer>();	
		
		Set<Integer> temp_set = new HashSet<Integer>();
		for(int i = 0; i < parent1.chromosome.length; i++){
			temp_set.add(parent1.chromosome[i]);
			remainList1.add(parent1.chromosome[i]);
		}
		
		int count = 0;
		for(int i = 0; i < parent2.chromosome.length; i++){
			try{
				if(!temp_set.add(parent2.chromosome[i])){
					remainList1.remove(remainList1.indexOf(parent2.chromosome[i]));
					chromosome1[count] = parent2.chromosome[i];
					chromosome2[count] = parent2.chromosome[i];
					count++;
				}
				else{
					remainList2.add(parent2.chromosome[i]);
				}
			}
			catch(Exception e){
				for(int a : remainList1){
					System.out.print(a + "\t");
				}
				System.out.println();
				for(int a : parent2.chromosome){
					System.out.print(a + "\t");
				}
				System.out.println();
			}
		}
		
		
		List<Integer> temp_list = new ArrayList<Integer>();	

		for(Integer a : remainList1){
			temp_list.add(a);
		}
		for(Integer a : remainList2){
			temp_list.add(a);
		}
		
		Collections.shuffle(temp_list);
		
		for(int i = 0; i < temp_list.size(); i+=2)
		{
			boolean wi = random.nextBoolean();
			if(wi == false){
				chromosome1[count] = temp_list.get(i);
				chromosome2[count] = temp_list.get(i+1);
			}
			else{
				chromosome2[count] = temp_list.get(i);
				chromosome1[count] = temp_list.get(i+1);
			}
			count++;
		}
		Object[] args = {chromosome1, this.A.get(0)};
		double[] fitness_result = (double[]) method.invoke(this, args);
		Chromosome speciesIndividual = new Chromosome(chromosome1, fitness_result);
		
		
		Object[] args2 = {chromosome2, this.A.get(0)};
		double[] fitness_result2 = (double[]) method.invoke(this, args2);
		Chromosome speciesIndividual2 = new Chromosome(chromosome2, fitness_result2);
		
		pop.add(speciesIndividual);
		pop.add(speciesIndividual2);
	}
	
	
	//Mutation operation
	public void mutate_SA(int index, double pm) throws Exception
	{
		double rate = Math.random();
		if(rate <= pm){
			mutate(index, false);
		}
	}
	
	
	//Mutation operation
	public void mutate(int Spindex, boolean ifAllIndex) throws Exception
	{
		int[] chromosome = new int[K];
		chromosome = pop.get(Spindex).chromosome;
		double init_fitness = pop.get(Spindex).fitness[0];
		
		int[] max_gene = new int[this.K];
		max_gene[this.K-1] = -1;
		
		int removeIndex = random.nextInt(65535) % K;
		int count2 = 0;
		for(int l = 0; l < chromosome.length; l++)
		{
			if(l != removeIndex){
				max_gene[count2++] = chromosome[l];
			}
		}
		
		//Finding candidate genes
		List<Integer> houxuan = new ArrayList<Integer>();
		boolean ifcunzai;
		for(int j = 0; j < name_index.length; j++){	
			ifcunzai = false;
			for(int k = 0; k < max_gene.length; k++){
				if(name_index[j] == max_gene[k]){
					ifcunzai = true;
					break;
				}
			}
			if(ifcunzai == false){
				houxuan.add(name_index[j]);
			}			
		}
		
		
		double maxfitness = -Double.MAX_VALUE;
		double tempfitness = 0.0;
		int maxIndex = 0;
		double[] max_result = null;
		
		//Randomly disrupt order
		Collections.shuffle(houxuan);
		
		//Jump out of the boundary of the loop
		int maxBreakCount = (int) Math.sqrt(houxuan.size());
//		searchCount = (int) (Math.log(Math.pow(houxuan.size(), k)) / Math.log(2));
		
		//Maximum number of searches
		int maxSearchesCount = houxuan.size();
		int step = 0;
	
		int breakCount = 0;
		while(step < maxSearchesCount){
			if(breakCount > maxBreakCount){
				break;
			}
			max_gene[this.K-1] = houxuan.get(step);
			Object[] args1 = {max_gene, this.A.get(0)};
			double[] tempfitness_result = (double[]) method.invoke(this, args1);
			tempfitness = tempfitness_result[0];
			if(maxfitness < tempfitness){
				maxfitness = tempfitness;
				maxIndex = step;
				max_result = tempfitness_result;
				breakCount = 0;
			}else{
				breakCount++;
			}
			if(maxfitness > init_fitness){
				break;
			}
			step++;
		}
		
		max_gene[this.K - 1] = houxuan.get(maxIndex);
		//update pop
		pop.get(Spindex).chromosome = max_gene;
		pop.get(Spindex).fitness = max_result;
	}
}










