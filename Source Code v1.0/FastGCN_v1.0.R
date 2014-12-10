 #############################################################################
 #FastGCN (a GPU accelerated tool for Fast Gene Co-expression Networks) 1.0   #
 #Copyright (c) 2010-2014 by Meimei Liang, Futao Zhang, Gulei Jin and Jun Zhu#
 #Institute of Bioinformatics, Zhejiang University, China                    #
 #First use Information Entropy to pre-process and cut off some genes.       # 
 #Then calculate the Pairwise Correlation Coefficients.                      #
 #Z-normalization of the coefficients                                        #
 #############################################################################

entropy_pre=function(eData,remain_para)
{
	#the first column is the gene names
	eData_tmp = eData[,-1]
	gNum=dim(eData_tmp)[1]
	sNum=dim(eData_tmp)[2]
	eData_tmp=abs(eData_tmp)
	#option 1
	rowSum=rowSums(eData_tmp)
	p=eData_tmp/rowSum
	logp=log(p)
	plogp=p*logp
	entropy=-rowSums(plogp)
	
	#option 2: use the code below, but it is too slow.
	#entropy=c()
	#for( i in 1:gNum)
	#{
  	#    expr_sum=sum(abs(eData_tmp[i,]))
  	#    pij=abs(eData_tmp[i,])/expr_sum
  	#    entropy=c(entropy,-sum(pij*log(pij)))
	#}
	
	sorted_entropy=sort(entropy)
	if(remain_para>1.0) pos=round(remain_para)
	else pos=ceiling(remain_para*gNum)
	
	threshold=sorted_entropy[pos]
	idx=which(entropy <= threshold)
	eData=eData[idx,]
	return(eData)	
}

cal_pcc=function(eData,USE_WGCNA)
{
	eData_tmp = eData[,-1]
	eData_tmp = as.matrix(eData_tmp)
	gNum=dim(eData_tmp)[1]
	sNum=dim(eData_tmp)[2]
	if(USE_WGCNA){
		print("using WGCNA package.")
		coeff=corFast(t(eData_tmp))
	}else{
		coeff=cor(t(eData_tmp),method="pearson")
	}
	# or use the code below, but it is too slow.
	#coeff=c()
	#for(i in 1:gNum)
	#{
	#	j=i+1
	#	while(j<=gNum)
	#	{
			
	#		tmp=cor(eData_tmp[i,],eData_tmp[j,], method="pearson")
		
	#		coeff=c(coeff,tmp)		
			
	#		j=j+1
	#	}
	#}	
	return(coeff)
}
getRow=function(vecid)
{
	rownum=round(sqrt(vecid))
	while(TRUE){
		if((rownum*(rownum-1)/2)>=vecid) break
		rownum=rownum+1
	}
	return(rownum)
}

USE_WGCNA=FALSE
if(USE_WGCNA){
 library(WGCNA)
 enableWGCNAThreads()
 #print(WGCNAnThreads())
}

remain_para=1
Z_threshold=3.29
FDR=0.05
FDR_switcher=TRUE
eData = read.table("eData_g2000.txt",head=TRUE)#the first row is the individual names

timestamp=Sys.time()

#step1: Filter with Genetic Entropy	
if(abs(remain_para-1.0)>1e-6) eData=entropy_pre(eData,remain_para)

#step2: calculate the PCC and get the coefficient matrix
coeff=cal_pcc(eData,USE_WGCNA)

#step3: Normalization and FDR control
coeff=coeff[upper.tri(coeff)]
abs_coeff=abs(coeff)
coeff_sd = sd(abs_coeff)
coeff_mean = mean(abs_coeff)
coeff_z=(abs_coeff-coeff_mean)/coeff_sd
pval=1-pnorm(coeff_z)
qval= p.adjust(pval,method="fdr",length(pval))

Sys.time()-timestamp

#step4 output
if(FDR_switcher)
{
	idx = which(qval<=FDR)
}else 
{
	idx=which(coeff_z>=Z_threshold)
}

rlt<-data.frame() 
if(length(idx)>0)
{
	for(i in 1:length(idx))
	{
		vec_id=idx[i]
		row_id=getRow(vec_id)
		col_id= vec_id-(row_id-1)*(row_id-2)/2
		gName1=as.character(eData[row_id,1])
		gName2=as.character(eData[col_id,1])
		coeff_pcc=coeff[vec_id]	
		z_value=coeff_z[vec_id]
		p_value=pval[vec_id]
		q_value=qval[vec_id]
		rlt[i,"gene1"]=gName1
		rlt[i,"gene2"]=gName2
		rlt[i,"coeff_pcc"]=coeff_pcc
		rlt[i,"z-score"]=z_value
		rlt[i,"p_value"]=p_value
		rlt[i,"q_value"]=q_value
	}
}
write.csv(rlt, "output_2k.csv")
