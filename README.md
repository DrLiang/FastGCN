FastGCN
=======

Features of the Software
First genetic entropy was exploited to filter out some genes with no mutation or small mutation in the raw data preprocessing step. Then Pearson Correlation Coefficients were calculated. At last we normalized these coefficients and used False Discovery Rate to control the multiple tests. All these calculating procedures were implemented on GPU. We also compressed the coefficient matrix to save the space. 
The R code was developed by Meimei Liang.The GPU code, OpenMP code for Multi-core computing and Single Thread C/C++ code were developed by Futao Zhang.

What are in FastGCN repository

	Source Code
		1. Source code of GPU implementation. 
		2. Source code of Multi-core CPU implementation. 
		3. Source code of Single-thread C/C++ implementation. 
		4. Source code of Single-thread R implementation.

		Version v1.1 (10/26/2014) Missing data Imputation and Modules Identification were added!
	
		Version v1.0 (07/01/2014)

	Simulation Datasets
		Data for test: test_data.txt
		Datasets for testing performance: http://ibi.zju.edu.cn/software/FastGCN/eData.zip

	Executable Binary Code Downloads
		1. GPU Implementation for Windows: x64 x86
		2. Multi-core Implementation for Windows: x64 x86
		3. C/C++ Implementation for Windows: x64 x86
		4. GPU Implementation for RHEL6.5: FastGCN GPU Linux
		5. Multi-core Implementation for RHEL6.5: FastGCN OpenMP Linux
		6. C/C++ Implementation for RHEL6.5: FastGCN C/C++ Linu

		Version v1.0 (07/01/2014)

System Requirements
Operation System: Windows, Linux and Mac OS.

Citation
Meimei Liang, Futao Zhang, Gulei Jin and Jun Zhu* (2014) FastGCN: A GPU Accelerated Tool for Fast Gene Co-expression Networks. (Submitted)