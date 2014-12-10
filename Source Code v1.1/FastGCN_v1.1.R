 #############################################################################
 #FastGCN (a GPU accelerated tool for Fast Gene Co-expression Networks) 1.1   #
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

# Imputation with Overall Mean of the whole values
OM_impute=function(eData)
{
	miss_idx=is.na(eData)
	if(length(which(miss_idx==TRUE))>0)
	{
		etmp=eData[,-1]
		etmp=as.vector(as.matrix(etmp))
		idx=which(is.na(etmp))
		etmp=etmp[-idx]
		iValue=mean(etmp)
		eData[miss_idx]=iValue	
	}
	return(eData)
}
# Imputation with Overall Mean of each Gene expression values
OMG_impute=function(eData)
{
	miss_idx=which(is.na(eData))
	if(length(miss_idx)>0)
	{
		geneNum=dim(eData)[1]
		for( i in 1:geneNum)
		{
			idx=which(is.na(eData[i,]))
			if(length(idx)>0)
			{
				vec=eData[i,]
				vec[idx]=0
				rMean=sum(vec[-1])/(length(vec)-length(idx)-1)
				eData[i,idx]=rMean
			}
		}
	}
	
	return(eData)
}

#Modules Identification. Input is the Pearson Correlation Coefficient Symmetric Matrix
#z-score as the evaluation criteria.
Modls_Idnt=function(coeff,gname)
{
	rlt=list()
	adja_mat=(abs(coeff)-coeff_mean)/coeff_sd
	diag(adja_mat)=0 # set the diagonal zero 
	idx=adja_mat>=Z_threshold
	adja_mat[idx]=1
	adja_mat[!idx]=0 #got the adjacency matrix
	K=rowSums(adja_mat) # got the node connectivity 
	idx=which(K==0)
	adja_mat=adja_mat[-idx,-idx] #remove the genes whose K is zero
	K=K[-idx]	
	rlt$gname=gname[-idx] #adjust the gene names
	H=adja_mat%*%adja_mat
	diag(H)=0
	H_A=H+adja_mat
	slt_num=length(K)
	K_mat=matrix(rep(K,slt_num),slt_num,slt_num)
	KI_KJ=K_mat+t(K_mat)
	Denominator=KI_KJ-H_A
	rlt$S=H_A/Denominator # got the similarity matrix
	return(rlt)
}

Modls_Idnt_=function(coeff)
{
	adja_mat=(abs(coeff)-coeff_mean)/coeff_sd # currently adja_mat is the z-score matrix
	diag(adja_mat)=0 # set the diagonal zero 
	idx=adja_mat>=Z_threshold
	adja_mat[idx]=1
	adja_mat[!idx]=0 #got the adjacency matrix
	K=rowSums(adja_mat) # got the node connectivity 
	H=adja_mat%*%adja_mat
	diag(H)=0
	H_A=H+adja_mat
	slt_num=length(K)
	K_mat=matrix(rep(K,slt_num),slt_num,slt_num)
	KI_KJ=K_mat+t(K_mat)
	Denominator=KI_KJ-H_A
	S=H_A/Denominator # got the similarity matrix
	return(S)
}

USE_WGCNA=FALSE
if(USE_WGCNA){
 library(WGCNA)
 enableWGCNAThreads()
 #print(WGCNAnThreads())
}

remain_para=1
Z_threshold=6.8
FDR=0.05
FDR_switcher=FALSE
eData = read.table("eData_g2000.txt",head=TRUE)#the first row is the individual names

timestamp=Sys.time()

#imputation with overall mean
eData=OMG_impute(eData)
#imputation completed


#step1: Filter with Genetic Entropy	
if(abs(remain_para-1.0)>1e-6) eData=entropy_pre(eData,remain_para)

#step2: calculate the PCC and get the coefficient matrix
coeff=cal_pcc(eData,USE_WGCNA)

#step3: Normalization and FDR control
coeff_vec=coeff[upper.tri(coeff)]
abs_coeff=abs(coeff_vec)
coeff_sd = sd(abs_coeff)
coeff_mean = mean(abs_coeff)
coeff_z=(abs_coeff-coeff_mean)/coeff_sd
pval=1-pnorm(coeff_z)
qval= p.adjust(pval,method="fdr",length(pval))



#step4: Modules Identification
gname=as.character(eData[,1])
rlt=Modls_Idnt(coeff,gname)
gname=rlt$gname
S=rlt$S

Sys.time()-timestamp

#step5 output the co-expression network
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
		coeff_pcc=coeff_vec[vec_id]	
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
