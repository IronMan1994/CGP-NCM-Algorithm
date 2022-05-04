package PSO_Resources;

import java.awt.print.Printable;
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
import java.util.List;
import java.util.Random;

import com.sun.javafx.image.impl.IntArgb;

import GA_Resources.Chromosome;


@SuppressWarnings({"all"})
public class PSO_Main
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
	
	//Current Particle Population
	public List<Particle> pop;
	//List of individual extremes for particles
	public List<Particle> piList;
	//Pg
	public Particle Pg;
	public List<Particle> nextPOP;
	//Statistics results after each algorithm execution
	public List<Particle_Index> maxGene;
	
	//Size of Particle Population
	public int sizeOfParticle;
	//maxg
	public int maxg;
	//gen
	public int gen;
		
	//Number of gene
	public int G;
	//K
	public int K;	
	//W(omega)
	public double W;
	//C1
	public double C1;
	//C2
	public double C2;
	//Mutation probability
	public double pm;
	//Velocity boundary
	public double Vm;
	
	//Record time
	long start;
	long end;

	
	//initialization
	public void initData(int maxg, int maxt, double omega, double C1, double C2, double Vm, double pm,
			String[] paths, String[] netPath, int G, int sizeOfParticle, int K, String modelName) throws Exception{
		
		this.G = G;
		this.K = K;
		this.sizeOfParticle = sizeOfParticle;
		this.maxg = maxg;
		this.W = W;
		this.C1 = C1;
		this.C2 = C2;
		this.gen = 0;
		this.Pg = new Particle();
		this.name = new String[this.G];
		this.name_index = new int[this.G];
		this.pop = new ArrayList<Particle>();
		this.piList = new ArrayList<Particle>();
		this.nextPOP = new ArrayList<Particle>();
		this.random = new Random(System.currentTimeMillis());
		this.A = new ArrayList<FileToList>();
		this.pm = pm;
		this.Vm = Vm;
		
		//File to list
		File f1 = new File("");
		for(int i = 0; i < paths.length;i++){
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
//						name[i] = ppp[i + 1].replaceAll("[^0-9a-zA-Z\u4e00-\u9fa5.，,。？“”]+","");
						name[i] = ppp[i + 1];
					}
					break;
				}
			}
			for(int i = 0; i < this.G; i++){
				name_index[i] = i;
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		bfr.close();
		reader.close();
		
		//Using reflection mechanism to find the required model
		Method[] methods = this.getClass().getMethods();
		for(Method m : methods){
			if(modelName.equals(m.getName())){
				this.method = m;
				break;
			}
		}
		
		//network to List
		netList = new ArrayList<HashMap<String,ArrayList<String>>>();
		for(int j = 0; j < paths.length; j++){
			recordNet(netPath[j]);
		}
	}

	//loading gene name from file
	public void readName(File f) throws Exception{
		//读取基因名字
		Reader reader = new FileReader(f);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		try{
			while((temp = bfr.readLine()) != null){
				String[] ppp = temp.split("\t");
				if(ppp.length == this.G + 1){
					for(int i = 0; i < this.G; i++){
						name[i] = ppp[i + 1].replaceAll("[^0-9a-zA-Z\u4e00-\u9fa5.，,。？“”]+","");
					}
					break;
				}
			}
			for(int i = 0; i < this.G; i++){
				name_index[i] = i;
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		finally {
			bfr.close();
			reader.close();
		}
	}
	
	
	//loading file to list
	public void readData(File f, List<float[]> B) throws Exception
	{
		InputStream fis = new FileInputStream(f);
		Reader isr = new InputStreamReader(fis);
		BufferedReader bfr = new BufferedReader(isr);
		String tempstr;
		int line = 1;
		try{
			while((tempstr = bfr.readLine()) != null){
				String[] ppp = tempstr.split("\t");
				if(ppp.length == this.G + 1){	
					if(line == 2){
						float[] temp = new float[this.G];
						for(int i = 0; i < this.G; i++){	
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
		Reader reader = new FileReader(f);
		BufferedReader bfr = new BufferedReader(reader);
		String temp;
		HashMap<String, ArrayList<String>> tempMap = new HashMap<String, ArrayList<String>>();
		try{
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
	public int calEdgeCount(Integer[] chromosome){
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
								result++;
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
	public Double[] model_NCM(Integer[] chromosome, FileToList B)
	{	
		Double[] result = new Double[3];
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
			edgeCount = 1 +  (edgeCount)/(maxEdgeCount);
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
	public Double[] model_GA(Integer[] chromosome, FileToList B)
	{
		//0-1模型
		Double[] result = new Double[3];
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
	
	
	//Construct initial particle population(Random Initialization)
	public void createBeginningSpeciesRandom() throws Exception{

		int count = 0;
		List<Integer> index = new ArrayList<Integer>();
		for(int i = 0; i < name_index.length; i++){
			index.add(i);
		}
		while(count < sizeOfParticle){
			Integer[] Xi = new Integer[G];			//Array of positions(Xi)
			Double[] Vi = new Double[G];			//Array of velocity(Vi)
			Double[] fitness_result;				
			Integer[] chromosome = new Integer[K];
			Double[] s = new Double[G];
			
			for(int i = 0; i < G; i++){
				Xi[i] = 0;
				Vi[i] = 0.0;
				s[i] = 0.0;
			}
			
			//Random disorder order
			Collections.shuffle(index);
			for(int i = 0; i < K; i++){
				chromosome[i] = index.get(i);
			}
			
			//Initialize velocity with probability of 0.5
			for(int i = 0; i < G; i++){
				double random_1_0 = Math.random() * (Vm);
				if(Math.random() < 0.5){
					Vi[i] = random_1_0;
				}else{
					Vi[i] = -random_1_0;
				}
			}
			
			//Calculating Fitness
			Object[] args = {chromosome, this.A.get(0)};
			fitness_result = (Double[]) method.invoke(this, args);
			Particle individual = new Particle(Xi, Vi, chromosome, fitness_result, count, s);
			pop.add(individual);
			count++;
		}//while
		
		//Sort by fitness from high to low
		particleTools.sort_pop_list(pop);
		//Record global optimum Pg
		Pg = particleTools.copy_particle(pop.get(0));
		//Recording the extreme Pi of an individual particle
		particleTools.copy_list(pop, piList);
	}
	
	
	//Update on particles
	public void Update() throws Exception{
		//1、迭代
		//2、计算适应值
		//3、更新Pi和Pg
		nextPOP.clear();
		for(int i = 0; i < pop.size(); i++){
			Double[] Vi_next = new Double[G];
			Integer[] Xi_next = new Integer[G];
			
			particleTools.copy_int(pop.get(i).Xi, Xi_next);
			
			Integer[] chromosome = new Integer[K];
			Double[] fitness = new Double[6];
			Double[] s = new Double[G];			

			double rand1 = random.nextDouble();
			double rand2 = random.nextDouble();
			
			//Update Omega
			W = this.W - (this.W-0.2)*(this.gen)/1000;
			
			//第一部分
			Double[] one = new Double[G];
			for(int j = 0; j < G; j++){
				one[j] = W * pop.get(i).Vi[j];
			}
			
			//第二部分、第三部分
			double[] Pi_Xi = new double[G];
			double[] Pg_Xi = new double[G];
			double C1_rand1 = C1 * rand1;
			double C2_rand2 = C2 * rand2;
			
			for(int j = 0; j < G; j++){
				Pi_Xi[j] = C1_rand1 * (piList.get(pop.get(i).index).Xi[j] - pop.get(i).Xi[j]);
				Pg_Xi[j] = C2_rand2 * (Pg.Xi[j] - pop.get(i).Xi[j]);
			}
			
			//V之间的加法为或运算
			for(int j = 0; j < G; j++){
				Vi_next[j] = one[j] + Pi_Xi[j] + Pg_Xi[j];
				if(Vi_next[j] > this.Vm){
					Vi_next[j] = this.Vm;
				}else if(Vi_next[j] < -this.Vm){
					Vi_next[j] = -this.Vm;
				}
				
				if(this.gen <= 0){
					s[j] = 1 / (1 + Math.pow(Math.E, -Vi_next[j]));
				}
				else {
					s[j] = Math.abs(1 - 2 / (1 + Math.pow(Math.E, -Vi_next[j])));
				}
			}
			
			//Update Xi
			for(int j = 0; j < G; j++){
				if(this.gen <= 0){
					if(random.nextDouble() <= s[j]){
						Xi_next[j] = 1;
					}else {
						Xi_next[j] = 0;
					}
				}
				else {
					if(Vi_next[j] < 0){
						if(random.nextDouble() <= s[j]){
							Xi_next[j] = 0;
						}
					}
					else if(Vi_next[j] > 0){
						if(random.nextDouble() <= s[j]){
							Xi_next[j] = 1;
						}
					}
				}
			}

			
			
			//Repairing the location of particles
			List<Integer> index_1 = new ArrayList<Integer>();
			
			for(int j = 0; j < G; j++){
				if(Xi_next[j] == 1){
					index_1.add(j);
				}
			}
			
			int complement_num = this.K - index_1.size();
			for(int j = 0; j < chromosome.length; j++){
				chromosome[j] = -1;
			}	
			
			//More than K locations with 1
			if(complement_num < 0){
				
				Collections.shuffle(index_1);
				
				for(int j = 0; j < this.K; j++){
					chromosome[j] = index_1.get(j);
				}
			}
			
			//Less than K locations of 1
			else if(complement_num > 0){
				int j = 0;
				while(j < complement_num){
					int gene_index = random.nextInt(65535) % this.G;
					int k = 0;
					for(; k < index_1.size(); k++){
						if(gene_index == index_1.get(k)){
							break;
						}
					}
					if(k == index_1.size()){
						index_1.add(gene_index);
						j++;
					}
				}//while
				for(int k = 0; k < this.K; k++){
					chromosome[k] = index_1.get(k);
				}
			}
			else {
				for(int j = 0; j < this.K; j++){
					chromosome[j] = index_1.get(j);
				}
			}
			
			//Update fitness
			Object[] args = {chromosome, this.A.get(0)};
			fitness = (Double[]) method.invoke(this, args);
			Particle next = new Particle(Xi_next, Vi_next, chromosome, fitness, pop.get(i).index, s);
			nextPOP.add(next);
		}//for
		
		//Update particle population
		particleTools.copy_list(nextPOP, pop);
			
		//Local search operation
		for(int i = 0; i < pop.size(); i++){
			mutate_SA(i);
		}
		
		//Update Xi、Vi、s、
		for(int j = 0; j < pop.size(); j++){
			Integer[] xi = new Integer[this.G];
			Double[] vi = new Double[this.G];
			Double[] s = new Double[this.G];
			
			for(int k = 0; k < G; k++){
				xi[k] = 0;
				vi[k] = pop.get(j).Vi[k];
				s[k] = pop.get(j).s[k];
			}
			
			for(int k = 0; k < pop.get(j).chromosome.length; k++){
				xi[pop.get(j).chromosome[k]] = 1;
				vi[pop.get(j).chromosome[k]] = 0.0;
				s[pop.get(j).chromosome[k]] = Math.abs(1 - 2 / (1 + Math.pow(Math.E, -vi[pop.get(j).chromosome[k]])));
			}
			
			pop.get(j).Xi = xi;
			pop.get(j).Vi = vi;
			pop.get(j).s = s;
			
			//Update Pi
			int index = pop.get(j).index;
			if(pop.get(j).fitness[0] > piList.get(index).fitness[0]){
				piList.set(index, particleTools.copy_particle(pop.get(j)));
			}
		}
		
		//Update Pg
		particleTools.sort_pop_list(pop);
		if(pop.get(0).fitness[0] > Pg.fitness[0]){
			Pg = particleTools.copy_particle(pop.get(0));	
		}
	}

	
	//Local search operation
	public void mutate_SA(int index) throws Exception
	{
		double rate = Math.random();
		if(rate <= this.pm){
			mutate(index, false);
		}
	}
	
	
	//Local search operation
	public void mutate(int Spindex, boolean ifAllIndex) throws Exception
	{
		Integer[] chromosome = new Integer[K];
		chromosome = piList.get(pop.get(Spindex).index).chromosome;
		double init_fitness = piList.get(pop.get(Spindex).index).fitness[0];

		Integer[] max_gene = new Integer[this.K];
		max_gene[this.K-1] = -1;
		
		int removeIndex = random.nextInt(65535) % K;
		int removeGene = chromosome[removeIndex];
		int count2 = 0;
		for(int l = 0; l < chromosome.length; l++)
		{
			if(l != removeIndex){
				max_gene[count2++] = chromosome[l];
			}
		}
		
		//Finding candidate genes
		List<Integer> houxuan = new ArrayList<Integer>();
		int count = 0;
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
				count++;
			}			
		}
		
		double maxfitness = -Double.MAX_VALUE;
		double tempfitness = 0.0;
		int maxIndex = 0;
		Double[] max_result = null;
		
		//Randomly disrupt order
		Collections.shuffle(houxuan);
		//Number of Searches
		int  searchCount = (int) Math.sqrt(houxuan.size());
		searchCount = (int) (Math.log(Math.pow(houxuan.size(), K)) / Math.log(2));
		
		//Local search operation
		int step = 0;
		while(step < searchCount){
			max_gene[this.K-1] = houxuan.get(step);
			Object[] args1 = {max_gene, this.A.get(0)};
			Double[] tempfitness_result = (Double[]) method.invoke(this, args1);
			tempfitness = tempfitness_result[0];
			if(maxfitness < tempfitness){
				maxfitness = tempfitness;
				maxIndex = step;
				max_result = tempfitness_result;
			}
			step++;
		}
		
		//Update the gene corresponding to Xi
		if(maxfitness > init_fitness){
			max_gene[this.K - 1] = houxuan.get(maxIndex);
			pop.get(Spindex).chromosome = max_gene;
			pop.get(Spindex).fitness = max_result;
			
		}
	}
	
	
	//Calculating Pvalue
	public double calPavlue(double max_fitness, int numberPvalue) throws Exception{
		List<Integer> index = new ArrayList<Integer>();
		for(int i = 0; i < name_index.length; i++){
			index.add(i);
		}
		int correct = 0;
		for(int i = 0; i < numberPvalue; i++){
			//随机算k个基因测适应值
			Integer[] chromosome = new Integer[this.K];
			Collections.shuffle(index);
			for(int j = 0; j < K; j++){
				chromosome[j] = index.get(j);
			}
			Object[] args = {chromosome, this.A.get(0)};
			Double[] fitness_result = (Double[]) this.method.invoke(this, args);
			if(max_fitness > fitness_result[0]){
				correct++;
			}
		}//for
		return 1 - (correct/numberPvalue);
	}
	
	
	//Run
	public void run(int numberPvalue, int maxg, int maxt, double omega, double C1, double C2, double Vm, double pm,
			String[] paths, String[] netPath, int G, int K, String modelName, 
			int numberAlgorithm, double standard) throws Exception{
		try
		{
			maxGene = new ArrayList<Particle_Index>();
			//Record running time
			double one_time = 0.0;
			double time_count = 0.0;
			for(int i = 0; i < numberAlgorithm; i++)
			{
				start = System.currentTimeMillis();	
				System.out.println(String.format("The particle swarm optimization is executed for the '%-2d' time to reconstruct the initial population......", i + 1));
				//initialization
				int sizeOfParticle = (int)Math.sqrt(G/2);
				initData(maxg, maxt, omega, C1, C2, Vm, pm, paths, netPath, G, sizeOfParticle, K, modelName);	
				pop.clear();
				piList.clear();
				nextPOP.clear();
				this.gen = 0;
				
				//Create initial population
				createBeginningSpeciesRandom();	
				
				int countMaxt = 0;
				double pre_fitness = pop.get(0).fitness[0];
				double next_fitness = 0.0;
				System.out.println("The optimal fitness value in the initial population is: " + pre_fitness);
		
				while(this.gen < this.maxg)
				{	
					if(countMaxt == maxt){
						break;
					}
					//Update
					Update();
					next_fitness = pop.get(0).fitness[0];
					if(pre_fitness == next_fitness){
						countMaxt++;
					}
					if(pre_fitness < next_fitness){
						pre_fitness = next_fitness;
						countMaxt = 0;
					}
					this.gen++;
				}
				end = System.currentTimeMillis();
				one_time = Double.valueOf(end - start)/1000;
				time_count += one_time;
				BigDecimal bd = new BigDecimal(one_time).setScale(4, RoundingMode.UP);
				
				System.out.println("The output gene set is: ");
				System.out.print("{ ");
				for(int j = 0; j < Pg.chromosome.length; j++){
					System.out.print(name[Pg.chromosome[j]] + "\t");
				}
				System.out.print("}");
				System.out.print("\nThe partial values are: ");
				for(int j = 0; j < Pg.fitness.length; j++){
					System.out.print(Pg.fitness[j] + "\t");
				}
				
				System.out.println("\nExecution time: " + bd.doubleValue() + "s");
				
				//Calculating Pvalue
				double Pvalue = calPavlue(Pg.fitness[0], numberPvalue);
				if(Pvalue <= 0.001){
					System.out.println("Pvalue <= " + 0.001 + "\n");
				}else{
					System.out.println("Pvalue <= " + Pvalue + "\n");
				}
				
				Particle temp = particleTools.copy_particle(Pg);				
				Particle_Index tempmaxgene = new Particle_Index(temp,i+1);
				maxGene.add(tempmaxgene);			
			}
			
			particleTools.sort_list_maxgene(maxGene);
			Particle_Index max = new Particle_Index();
			max = maxGene.get(0);
			
			//print
			System.out.println("The '" + max.index + "' performing algorithm yields the optimal geneset as: ");
			System.out.print("{ ");
			for(int i =  0; i < max.speciesIndividual.chromosome.length; i++)
			{
				System.out.print(name[max.speciesIndividual.chromosome[i]] + "\t");
			}
			System.out.print("}");
			
			//According to the model corresponding to different fitness values (e.g. NE and CM in model_NCM respectively)
			System.out.print("\nThe partial values are: ");
			for(int i = 0; i < max.speciesIndividual.fitness.length; i++)
			{
				System.out.print(max.speciesIndividual.fitness[i] + "\t");
			}
			
			//average time
			BigDecimal bd = new BigDecimal(time_count/numberAlgorithm).setScale(4, RoundingMode.UP);
			System.out.println("\n\n"+"The average time of "+numberAlgorithm+" executions is: " + bd.doubleValue() + "s");
			System.out.println();
			
			//Print the results of each algorithm execution
			for(int j = 0; j < maxGene.size(); j++){
				System.out.print((maxGene.get(j).index) + "\t");
				for(int i =  0; i < maxGene.get(j).speciesIndividual.chromosome.length; i++){
					System.out.print(name[maxGene.get(j).speciesIndividual.chromosome[i]] + "\t");
				}
				System.out.print(maxGene.get(j).speciesIndividual.fitness[0]);
				System.out.println();
			}	
			
			//Calculating algorithm accuracy
			System.out.println("\nThe accuracy of the algorithm is: " + particleTools.calAccuracy(standard, maxGene) + "%");
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	//Run_Main
	public static void main(String[] args) throws Exception
	{
		PSO_Main psoMain = new PSO_Main();
		
		//Number of times to test P-value
		int numberPvalue = 1000;
		//maxg
		int maxg = 1000;
		//maxt
		int maxt = 10;
		//W(omega)
		double omega = 0.4;
		//C1=C2
		double C1=1.0, C2 = 1.0;
		//Vm
		double Vm = 10;
		//pm
		double pm = 1.0;
		
		//Modifiable parameters
		
		//input file
		String path = "A_cut_off-2-sample_cut_off-0_GBM-geneCount-440_TEST.txt;";
		String[] paths = path.split(";");
		//network file
		String[] netPath = {"mergeNet.txt"};
		
		//gene number
		int G = 440;

		//K
		int K = 10;	
		
		//moudle (model_GA、model_NCM)
		String modelName = "model_GA";
		
		//Number of times algorithm is repeated
		int numberAlgorithm = 10;
		
		//The file's standard optimal fit under the model used to calculate the accuracy of the algorithm
		double standard = 124.7525;
		
		psoMain.run(numberPvalue, maxg, maxt, omega, C1, C2, Vm, pm, paths, netPath, G, K, modelName, numberAlgorithm, standard);
	}
}










