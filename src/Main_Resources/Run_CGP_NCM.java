package Main_Resources;

import java.awt.print.Printable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import GA_Resources.Chromosome;
import GA_Resources.Chromosome_Index;
import GA_Resources.GA_Main;
import GA_Resources.chromosomeTools;
import PSO_Resources.PSO_Main;
import PSO_Resources.Particle;
import PSO_Resources.particleTools;

@SuppressWarnings("unchecked")
public class Run_CGP_NCM
{
	//时间
	public int gen;
	public long start;
	public long end;
	public List<Chromosome> Cross_pool;
	public List<Chromosome_Index> maxGene;
	public int[] breakArray;
	public double[] preFitnessArray;
	
	
	public void run(int numberPvalue, int maxg, int maxt, double omega, double C1, double C2, double Vm, double pmParticle,
			double pmChromosome, String[] paths, String[] netPath, int G, int K, String modelName, 
			int numberAlgorithm, double standard)
	{
		GA_Main t = new GA_Main();
		GA_Main t2 = new GA_Main();
		PSO_Main t3 = new PSO_Main();
		maxGene = new ArrayList<Chromosome_Index>();
		
		String initstr = String.format("%s   %s   %s   %s   %s   %s   %s   %s   %s   %s   %s", paths[0], netPath[0], "modelName = " + modelName, "K = " + K, "Ng = " + (int)G/4, "pmChromosome = " + pmChromosome, "Np = " + (int)Math.sqrt(G/2), "w = " + omega, "C1=C2 = " + C1, "numberAlgorithm = " + numberAlgorithm, "standard = " + standard);
		
		try
		{
			List results = new ArrayList<>();
			double one_time = 0.0;
			double time_count = 0.0;
			for(int i = 0; i < numberAlgorithm; i++)
			{
				Cross_pool = new ArrayList<Chromosome>();
				breakArray = new int[3];
				preFitnessArray = new double[3];
				
				start = System.currentTimeMillis();	
				System.out.println(String.format("The CGP-NCM algorithm is executed for the ===== %d ====== time to reconstruct the initial population......", i + 1));
				
				//initialization
				gen = 0;
				int sizeOfChromosome = (int)G/4;
				int sizeOfParticle = (int)Math.sqrt(G/2);
				t.initData(paths, netPath, G, sizeOfChromosome, K, pmChromosome, modelName);
				t2.initData(paths, netPath, G, sizeOfChromosome, K, pmChromosome, modelName);
				t3.initData(maxg, maxt, omega, C1, C2, Vm, pmParticle, paths, netPath, G, sizeOfParticle, K, modelName);
				
				//ThreadPool
				ExecutorService es = Executors.newFixedThreadPool(3);
				
				//Thread task
				Callable<Boolean> callable1 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							int j = 0;
							double[] rate_array = chromosomeTools.crossover_rate(Cross_pool);
							double[] tRate = chromosomeTools.crossover_rate(t.pop);
				
							while(j < t.sizeOfChromosome / 4){
								t.crossover(t.pop, tRate);
								t.crossover(Cross_pool, rate_array);
								j++;
							}													
							j=0;
							while(j < t.pop.size() - t.sizeOfChromosome)
							{	
								t.mutate_SA(t.sizeOfChromosome + j, pmChromosome);
								j++;
							}
							t.select();
							return true;
						}
						catch(Exception e){
							e.printStackTrace();
							return false;
						}	
					}
				};		
				
				Callable<Boolean> callable2 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							int j = 0;
							double[] rate_array = chromosomeTools.crossover_rate(Cross_pool);
							double[] t2Rate = chromosomeTools.crossover_rate(t2.pop);
							
							while(j < t2.sizeOfChromosome / 4){
								t2.crossover(t2.pop, t2Rate);
								t2.crossover(Cross_pool, rate_array);
								j++;
							}
							j=0;
							while(j < t2.pop.size() - t2.sizeOfChromosome)
							{
								t2.mutate_SA(t2.sizeOfChromosome + j, pmChromosome);
								j++;
							}
							t2.select();
							return true;
						}
						catch(Exception e){
							e.printStackTrace();
							return false;
						}
					}
				};
					
				Callable<Boolean> callable3 = new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception
					{
						try{
							t3.Update();
							return true;
						}
						catch(Exception e){
							e.printStackTrace();
							return false;
						}	
					}
				};		
				
				
				Future<Boolean> r1;
				Future<Boolean> r2;
				Future<Boolean> r3;
				Boolean if_complete1;
				Boolean if_complete2;
				Boolean if_complete3; 
		
				
				//Create initial population(Created by respective thread)
				Future<Boolean> c1 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception{
						t.createBeginningSpeciesRandom();
						return true;
					}
				});
				
				Future<Boolean> c2 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception{
						t2.createBeginningSpeciesRandom();
						return true;
					}
				});
				
				Future<Boolean> c3 = es.submit(new Callable<Boolean>() {
					@Override
					public Boolean call() throws Exception{
						t3.createBeginningSpeciesRandom();
						return true;
					}
				});
				
				
				c1.get();
				c2.get();
				c3.get();
				
				
				chromosomeTools.sort_pop_list(t.pop);
				chromosomeTools.sort_pop_list(t2.pop);
							
				preFitnessArray[0] = t.pop.get(0).fitness[0];
				preFitnessArray[1] = t2.pop.get(0).fitness[0];
				preFitnessArray[2] = t3.pop.get(0).fitness[0];
				
				//The primary population is put into the cooperative pool
				mainTools.copy_list(t.pop, Cross_pool);
				mainTools.copy_list(t2.pop, Cross_pool);
				mainTools.copy_list3(t3.piList, Cross_pool);
				
				//append Particle population into Genetic population
				mainTools.copy_list3(t3.pop, t.pop);
				mainTools.copy_list3(t3.pop, t2.pop);
				
				
				//Take out the best fitness of the three populations
				double[] bestFitnessAndIndex123 = mainTools.bestFitnessAndIndex123(t, t2, t3, breakArray, preFitnessArray);
				double index123 = bestFitnessAndIndex123[1];
				
				Chromosome temp2 = new Chromosome();
				if(index123 == 1){
					temp2 = mainTools.copy_SpeciesIndividual(t.pop.get(0));
				}
				else if(index123 == 2){
					temp2 = mainTools.copy_SpeciesIndividual(t2.pop.get(0));
				}
				else{
					temp2 = mainTools.copy_SpeciesIndividual(t3.Pg);
				}
							
				t3.Pg = particleTools.copy_individual_To_particle(temp2, t3.G, t3.Vm);
				
				System.out.println("The optimal fitness value in the initial population is: "  + temp2.fitness[0]);		
				
				double preFitness = bestFitnessAndIndex123[0];
				int countMaxt = 0;
				
				while(gen < maxg)
				{				
					if(countMaxt == maxt){
						break;
					}
										
					r1 = es.submit(callable1);
					r2 = es.submit(callable2);
					r3 = es.submit(callable3);
					
					if_complete1 = (Boolean) r1.get();
					if_complete2 = (Boolean) r2.get();
					if_complete3 = (Boolean) r3.get();
					
					//The third cooperation strategy
					if(if_complete1 && if_complete2){
						int popcount = t.sizeOfChromosome - 1;
						for(int j = 0; j < 1; j++){
							Chromosome tempIndividual = chromosomeTools.copy_SpeciesIndividual(t.pop.get(j));
							Chromosome tempIndividual2 = chromosomeTools.copy_SpeciesIndividual(t2.pop.get(j));
							
							if(tempIndividual.fitness[0] > t2.pop.get(popcount).fitness[0]){
								t2.pop.set(popcount, tempIndividual);
							}
							if(tempIndividual2.fitness[0] > t.pop.get(popcount).fitness[0]){
								t.pop.set(popcount, tempIndividual2);
							}
							popcount--;
						}
					}
					
					
					chromosomeTools.sort_pop_list(t.pop);
					chromosomeTools.sort_pop_list(t2.pop);
					
					bestFitnessAndIndex123 = mainTools.bestFitnessAndIndex123(t, t2, t3, breakArray, preFitnessArray);
					double Index123 = bestFitnessAndIndex123[1];
					
					mainTools.share(Index123, t, t2, t3);
		
					mainTools.copy_list(t.pop, Cross_pool);
					mainTools.copy_list(t2.pop, Cross_pool);
					mainTools.copy_list3(t3.pop, Cross_pool);

					mainTools.copy_list3(t3.pop, t.pop);
					mainTools.copy_list3(t3.pop, t2.pop);							
					
					//update cooperative pool
					mainTools.update_Cross_pool(t, (sizeOfChromosome * 2 + t3.sizeOfParticle), Cross_pool);
					
					if(preFitness < bestFitnessAndIndex123[0]){
						preFitness = bestFitnessAndIndex123[0];
						countMaxt = 0;
						
					}else{
						countMaxt++;
					}
					gen++;
				}//while
				es.shutdown();
				
				bestFitnessAndIndex123 = mainTools.bestFitnessAndIndex123(t, t2, t3, breakArray, preFitnessArray);
				double Index123 = bestFitnessAndIndex123[1];
				
				//record best chromosome
				Chromosome temp = new Chromosome();
				if(Index123 == 1){
					temp = mainTools.copy_SpeciesIndividual(t.pop.get(0));
				}
				else if(Index123 == 2){
					temp = mainTools.copy_SpeciesIndividual(t2.pop.get(0));
				}
				else{
					temp = mainTools.copy_SpeciesIndividual(t3.Pg);
				}
				
				//print result
				System.out.println("The output gene set is: ");
				System.out.print("{ ");
				for(int j = 0; j < temp.chromosome.length; j++){
					System.out.print(t.name[temp.chromosome[j]] + "\t");
				}
				System.out.print("}" + " Fitness = " + temp.fitness[0]);
				System.out.print("\n" + "The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: ");
				System.out.println(t.pop.get(0).fitness[0] + "\t" + t2.pop.get(0).fitness[0] + "\t" + t3.Pg.fitness[0]);
				System.out.print("The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: ");
				System.out.println(Arrays.toString(breakArray));
				
				//record executon time
				end = System.currentTimeMillis();
				one_time = Double.valueOf(end - start)/1000;
				time_count += one_time;
				BigDecimal bd = new BigDecimal(one_time).setScale(4, RoundingMode.UP);
				System.out.println("Execution time: " + bd.doubleValue() + "s");
				
				//Calculating Pvalue
				double Pvalue = t3.calPavlue(temp.fitness[0], numberPvalue);
				if(Pvalue <= 0.001){
					System.out.println("Pvalue <= " + 0.001 + "\n");
				}else{
					System.out.println("Pvalue <= " + Pvalue + "\n");
				}
				
				//record to file
				String str = "";
				str += String.format("%-20s", (i+1) + "：" + String.format("%.04f", bd.doubleValue()) + "s");;
				for(int j =  0; j < temp.chromosome.length; j++){
					str += String.format("%-10s", t.name[temp.chromosome[j]]);
				}
				str += temp.fitness[0] + " ";
				results.add(str);
				
				//put best to maxGene
				Chromosome_Index tempmaxgene = new Chromosome_Index(temp, i+1, Index123);
				maxGene.add(tempmaxgene);	
			}
			Chromosome_Index max = new Chromosome_Index();
			mainTools.sort_list_maxgene(maxGene);
			max = maxGene.get(0);

			//print the final result
			BigDecimal bd = new BigDecimal(time_count/numberAlgorithm).setScale(4, RoundingMode.UP);
			System.out.println(String.format("The '" + "===== %d ======" + "' performing algorithm yields the optimal geneset as: ", max.index));
			System.out.print("{ ");
			for(int i =  0; i < max.speciesIndividual.chromosome.length; i++)
			{
				System.out.print(t.name[max.speciesIndividual.chromosome[i]] + "\t");
			}
			System.out.print("}");
			
			//According to the model corresponding to different fitness values (e.g. NE and CM in model_NCM respectively)
			System.out.println("\nThe partial values are: ");
			for(int i = 0; i < max.speciesIndividual.fitness.length; i++)
			{
				System.out.print(max.speciesIndividual.fitness[i] + "\t");
			}
			
			String str = "The average time of "+numberAlgorithm+" executions is: " + bd.doubleValue() + "s";
			
			results.add(str);
			results.add("The accuracy of the algorithm is: " + mainTools.calAccuracy(standard, maxGene));
			
			//output to file
			mainTools.writeResult(results, initstr);
			
			System.out.println("\n\n"+"The average time of "+numberAlgorithm+" executions is: " + bd.doubleValue() + "s" + "\n");
			
			for(int j = 0; j < maxGene.size(); j++){
				System.out.print((j+1) + "\t");
				for(int i =  0; i < maxGene.get(j).speciesIndividual.chromosome.length; i++)
				{
					System.out.print(t.name[maxGene.get(j).speciesIndividual.chromosome[i]] + "\t");
				}
				System.out.print(maxGene.get(j).speciesIndividual.fitness[0] + "\t" + maxGene.get(j).Index123);
				System.out.println();
			}
			
			System.out.println("\n" + "The gen is: " + gen);
			
			System.out.println("The accuracy of the algorithm is: " + mainTools.calAccuracy(standard, maxGene) + "%");
			
			//测P值
//			double max_fitness = max.speciesIndividual.fitness[0];
//			int correct = 0;
//			
//			for(int i = 0; i < P_step; i++){
//				//随机算k个基因测适应值
//				int[] chromosome = new int[t.k];
//				int index = t.random.nextInt(65535) % (t.name_index.length);
//				chromosome[0] = index;
//				int j = 1;
//				for(; j < t.k;){
//					index = t.random.nextInt(65535) % (t.name_index.length);
//					int m = 0;
//					for(; m < t.k; m++){
//						if(chromosome[m] == index){
//							break;
//						}
//					}
//					//表示没有重复
//					if(m == t.k){
//						chromosome[j] = index;
//						j++;
//					}
//				}
//				double temp_fitness = t.calfitness(chromosome)[0];
//				if(max_fitness > temp_fitness){
//					correct++;
//				}
//			}//for 1000
//			System.out.println("P值为：" + (double)correct / P_step);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args)
	{
		Run_CGP_NCM run_CGP_NCM = new Run_CGP_NCM();
		
		
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
		//pm of Particle (By default, each particle performs a local search operation)
		double pmParticle = 1.0;
		
		//pm of Chromosome
		double pmChromosome = 0.01;
		
		
		//Modifiable parameters
		
		//input file
		String path = "dataFile//A_cut_off-2-sample_cut_off-0_OV-geneCount-2547.txt;";
		String[] paths = path.split(";");
		//network file
		String[] netPath = {"networkFile//mergeNet.txt"};
		
		//gene number
		int G = 2547;

		//K
		int K = 10;	
		
		//moudle (model_GA、model_NCM)
		String modelName = "model_NCM";
		
		//Number of times algorithm is repeated
		int numberAlgorithm = 10;
		
		//The file's standard optimal fit under the model used to calculate the accuracy of the algorithm
		double standard = 124.7525;
		
		run_CGP_NCM.run(numberPvalue, maxg, maxt, omega, C1, C2, Vm, pmParticle, pmChromosome, paths, netPath, G, K, modelName, numberAlgorithm, standard);		
	}
}







