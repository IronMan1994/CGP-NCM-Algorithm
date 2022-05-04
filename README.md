# CGP-NCM-Algorithm

# `JAVA` source code of `CGP-NCM` algorithm

## A brief description of algorithm

* INPUT: a weighted non-binary mutation matrix `A`, a PPI network `Q`, a parameter `K`;
* OUTPUT: a set of genes corresponding to submatrix `M`;

## Example of `A` input to algorithm

| Gene | TP53 | CDKN2A | CDKN2B| RB1 | CDK4| … |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| Sample_1 | 1 | 1.5 | 0 | 0 | 1 | … |
| Sample_2 | 0.45 | 1 | 0 | 0.4 | 1.5 | … |
| Sample_3 | 0 | 1.5 | 0 | 0.3 | 0 | … |
| … | … | … | … | … | … | … |


## Example of `Q` file input to algorithm, the first item in each row is a gene, and the remaining items in the row are connected to this gene.
**RNF7** UBE2Q1   UBE2D4   UBE2D2   UBE2D1

**RNF26**	UBE2Q1	UBE2D2	UBE2D3	UBE2D1	UBE2U	UBE2W	UBE2E1	UBE2D4	UBE2G2

**EXOSC7**	UBE2Q1	UBE2Q2	IP6K1	EXOSC9

**...**  ...


## The process of executing the project

1. You need to download `CGP-NCM-ALgorithm.zip` first and unzip it.

2. You need to import the `My_Cancer_GA_PSO_GitHub` folder in `CGP-NCM-ALgorithm` folder into `eclipse` or `MyEclipse` and execute them in the `JAVA8` environment whenever it is possible.
   
3. The `main` method in the `Run_CGP_NCM.java` in package `Main_Resources` is the entry to the whole program. The 'A' matrix corresponding to the three cancer species tested in the paper is placed in package `dataFile`, and the PPI network for experiment is placed in package `networkFile`, which has been preprocessed into executable format.
  
4. Enter the relative or absolute path of the `txt` file in the following statement.
    ```java
       //path of `A` matrix file
       String path = "dataFile//A_cut_off-2-sample_cut_off-0_GBM-geneCount-440.txt;";
   ```
   ```java
      //path of network file
      String[] netPath = {"networkFile//mergeNet.txt"};
   ```


   
5. Setting parameters.
   ### * `G`: The number of genes in the input `A` matrix. If you use the three files that come with the dataFile in the project, the number following the last `-'symbol in the file name is the number of genes.
   ### * `K`: Set the scale of driver pathway that need to be looked for.
   ### * `modelName`: The project provides two model choices, `model_GA` stands for the model of method `Dentrix`, used to test algorithm time and accuracy, `model_NCM` represents the model proposed in the paper and is used to find the driver pathway using matrix `A` and PPI network `Q`. By simply setting the model name, the project automatically matches the model using a reflection mechanism.
   ### * `numberAlgorithm`: The project supports repeated execution of `CGP-NCM` algorithm, and this parameter is set to meet the needs of repeated execution.
   ### * `standard`: This parameter represents the fitness value of the accurate solution represented by the input file calculated through the set model, which is used to test the accuracy of the algorithm. If you don't know this value, you don't need to set it, so the accuracy of output has no reference value. The self-contained in the code is the adaptation value of the exact solution under the default setting.
   ### * Other parameters have been set as the optimal value through a large number of experiments and can be kept unchanged for testing.

6. After setting the parameters, `CGP-NCM` algorithm is ready to be executed.
7. This project provides the function of outputting results to file. After running, you can view file `Results.txt` in package `results`.
8. After correct execution, the results are printed in the console as follows:
      ```
      The CGP-NCM algorithm is executed for the ===== 1 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 64.9982
      The output gene set is: 
      { CDKN2A	CDK4	EGFR	ERBB2	TP53	STAT1	RB1	PIK3R1	VAV1	MDM2	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	121.2334
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 4]
      Execution time: 1.851s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 2 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 67.5441
      The output gene set is: 
      { TP53	CDK4	CDKN2A	ERBB2	EGFR	RB1	MDM2	STAT1	PIK3R1	VAV1	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	119.2783
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 4]
      Execution time: 0.863s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 3 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 78.2701
      The output gene set is: 
      { TP53	CDK4	EGFR	PIK3R1	CDKN2A	MDM2	STAT1	ERBB2	VAV1	RB1	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	120.7249
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 6]
      Execution time: 0.7551s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 4 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 66.8354
      The output gene set is: 
      { EGFR	TP53	CDK4	STAT1	PIK3R1	MDM2	RB1	ERBB2	VAV1	CDKN2A	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	121.2334
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 1]
      Execution time: 0.865s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 5 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 69.5855
      The output gene set is: 
      { EGFR	CDK4	TP53	PIK3R1	CDKN2A	RB1	ERBB2	STAT1	MDM2	VAV1	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	118.1122
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 2]
      Execution time: 0.5621s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 6 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 75.6741
      The output gene set is: 
      { CDK4	TP53	RB1	CDKN2A	EGFR	MDM2	MDM4	STAT1	PIK3R1	ERBB2	} Fitness = 120.7249
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 120.7249	120.7249	117.9842
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 4]
      Execution time: 0.6471s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 7 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 71.8921
      The output gene set is: 
      { TP53	MDM2	EGFR	STAT1	ERBB2	CDKN2A	RB1	CDK4	PIK3R1	MDM4	} Fitness = 120.7249
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 120.7249	120.7249	120.7249
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [10, 10, 11]
      Execution time: 0.9111s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 8 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 68.9963
      The output gene set is: 
      { TP53	CDKN2A	PIK3R1	MDM2	EGFR	CDK4	ERBB2	RB1	VAV1	STAT1	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	119.2874
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 4]
      Execution time: 0.864s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 9 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 72.9373
      The output gene set is: 
      { TP53	CDKN2A	EGFR	CDK4	RB1	ERBB2	PIK3R1	STAT1	MDM2	VAV1	} Fitness = 121.2334
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 121.2334	121.2334	120.4014
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 9]
      Execution time: 0.7621s
      Pvalue <= 0.001

      The CGP-NCM algorithm is executed for the ===== 10 ====== time to reconstruct the initial population......
      The optimal fitness value in the initial population is: 70.4337
      The output gene set is: 
      { CDKN2A	TP53	MDM2	MDM4	CDK4	RB1	EGFR	ERBB2	PIK3R1	STAT1	} Fitness = 120.7249
      The optimal 'fitness' of 'GA_1', 'GA_2' and 'PSO_3' populations are: 120.7249	120.7249	118.0887
      The jumping steps of 'GA_1', 'GA_2' and 'PSO_3' populations are: [11, 11, 2]
      Execution time: 0.6461s
      Pvalue <= 0.001

      The '===== 1 ======' performing algorithm yields the optimal geneset as: 
      { CDKN2A	CDK4	EGFR	ERBB2	TP53	STAT1	RB1	PIK3R1	VAV1	MDM2	}
      The partial values are: 
      121.2334	1.3334	90.9205	

      The average time of 10 executions is: 0.8727s

      1	CDKN2A	CDK4	EGFR	ERBB2	TP53	STAT1	RB1	PIK3R1	VAV1	MDM2	121.2334	3.0
      2	TP53	CDK4	CDKN2A	ERBB2	EGFR	RB1	MDM2	STAT1	PIK3R1	VAV1	121.2334	2.0
      3	TP53	CDK4	EGFR	PIK3R1	CDKN2A	MDM2	STAT1	ERBB2	VAV1	RB1	121.2334	2.0
      4	EGFR	TP53	CDK4	STAT1	PIK3R1	MDM2	RB1	ERBB2	VAV1	CDKN2A	121.2334	3.0
      5	EGFR	CDK4	TP53	PIK3R1	CDKN2A	RB1	ERBB2	STAT1	MDM2	VAV1	121.2334	2.0
      6	TP53	CDKN2A	PIK3R1	MDM2	EGFR	CDK4	ERBB2	RB1	VAV1	STAT1	121.2334	2.0
      7	TP53	CDKN2A	EGFR	CDK4	RB1	ERBB2	PIK3R1	STAT1	MDM2	VAV1	121.2334	2.0
      8	CDK4	TP53	RB1	CDKN2A	EGFR	MDM2	MDM4	STAT1	PIK3R1	ERBB2	120.7249	2.0
      9	TP53	MDM2	EGFR	STAT1	ERBB2	CDKN2A	RB1	CDK4	PIK3R1	MDM4	120.7249	3.0
      10	CDKN2A	TP53	MDM2	MDM4	CDK4	RB1	EGFR	ERBB2	PIK3R1	STAT1	120.7249	2.0

      The gen is: 32
      The accuracy of the algorithm is: 99.88%
   ```
9. The results output in the file are shown as follows:
   ```
      dataFile//A_cut_off-2-sample_cut_off-0_GBM-geneCount-440.txt   networkFile//mergeNet.txt   modelName = model_NCM   K = 10   Ng = 110   pmChromosome = 0.01   Np = 14   w = 0.4   C1=C2 = 1.0   numberAlgorithm = 10   standard = 121.2334：
            1：1.3441s           CDK4      TP53      CDKN2A    MDM2      RB1       EGFR      ERBB2     PIK3R1    STAT1     VAV1      121.2334 
            2：0.6130s           TP53      CDKN2A    CDK4      MDM2      RB1       MDM4      EGFR      STAT1     ERBB2     PIK3R1    120.7249 
            3：1.1960s           CDK4      TP53      EGFR      PIK3R1    ERBB2     STAT1     CDKN2A    RB1       VAV1      MDM2      121.2334 
            4：0.6801s           EGFR      TP53      RB1       PIK3R1    MDM2      STAT1     CDKN2A    VAV1      CDK4      ERBB2     121.2334 
            5：0.6150s           EGFR      CDK4      TP53      ERBB2     STAT1     RB1       PIK3R1    MDM2      VAV1      CDKN2A    121.2334 
            6：0.7591s           TP53      CDK4      RB1       CDKN2A    EGFR      MDM2      STAT1     PIK3R1    ERBB2     VAV1      121.2334 
            7：0.6491s           TP53      CDKN2A    RB1       CDK4      EGFR      STAT1     ERBB2     PIK3R1    MDM2      VAV1      121.2334 
            8：0.9520s           CDKN2A    TP53      EGFR      MDM2      STAT1     ERBB2     PIK3R1    CDK4      VAV1      RB1       121.2334 
            9：0.6861s           CDKN2A    TP53      EGFR      CDK4      MDM2      RB1       PIK3R1    ERBB2     VAV1      STAT1     121.2334 
            10：0.6050s          TP53      CDKN2A    CDK4      ERBB2     EGFR      PIK3R1    VAV1      RB1       MDM2      STAT1     121.2334 
            The average time of 10 executions is: 0.81s
            The accuracy of the algorithm is: 99.96
   ```

  
## Some supplementary notes

* If you input other custom file, please check the file format and adjust the parameters.
* In addition, the `main` method in the `PSO_Main.java` file in package `PSO_Resources` is the entry of the particle swarm optimization algorithm, which can be used to perform the particle swarm optimization method separately for the adapted maximum weight submatrix problem.
* The entire project can be imported with reference to the following steps (compiler `MyEclipse`):
   ![image](https://user-images.githubusercontent.com/23414623/166691738-15ce6a1e-ea6a-428b-a857-2683b8573626.png)
   ![image](https://user-images.githubusercontent.com/23414623/166691819-93ac1ec9-6c5e-4739-9148-204c7dc7e92e.png)
   ![image](https://user-images.githubusercontent.com/23414623/166691855-7ea1747d-acbc-4fa2-9661-5da5fe4cc7d5.png)
   ![image](https://user-images.githubusercontent.com/23414623/166692130-144eb3dc-51ad-4421-bb2b-37fb25645046.png)
   ![image](https://user-images.githubusercontent.com/23414623/166692225-4263a8f0-4df0-4381-a451-a8a85d112ade.png)
   ![image](https://user-images.githubusercontent.com/23414623/166692466-0bffda51-c148-4bd1-93be-f9f5dd9e2247.png)






   
  
