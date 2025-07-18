
setwd("C:/Users/ et etc")

pheno=read.table("pheno.txt", header=T)       #read phenotypes
dim(pheno)                                    #check how many
summary(pheno)                                #check the values
head(pheno)

geno = as.matrix(read.table("geno50k.txt"))    #read genotypes
dim(geno)                                     #check how many
head(geno[1:10,1:10])

map=read.table("map.txt", header=T)         #read the map file, with poistion of each SNP
dim(map)
nmarkers=nrow(map)                         #check that the number fits the number in genotype file       

#code is similar to prcatical 6, just realize that genotypes per SNP are now in columns, not in rows

p=colMeans(geno)/2		 # calculate allele frequency (p) for every SNP, p=sum/2
geno=t(geno)-2*p 		 # adjust genotypes for column mean (-2p), column centres at mean
geno=t(geno)			        # the t() function takes a transpose, used here to make arrays fit
geno=as.matrix(geno)


#Linear regression  GWAS 
snp_effects=numeric(ncol(geno))
p_values=numeric(ncol(geno))

for (i in 1:ncol(geno))     #note that this loop can take a minute
{
  LM1=summary(lm(pheno[,1]~geno[,i]))
  snp_effects[i]=LM1$coeff[2,1]
  p_values[i]=LM1$coeff[2,4]
}

# plot for the most significant SNP the different genotypes vs their phenotypes 
plot(geno[,which(p_values==min(p_values))],pheno$phenotype,
     xlab="genotypes", ylab="phenotypes",
     main=paste("effect size:",round(snp_effects[which(p_values==min(p_values))],2)))
mod=lm(pheno$phenotype~geno[,which(p_values==min(p_values))])
abline(mod,lwd=2,col="blue") # add regression line


# plot for the least significant SNP the different genotypes vs their phenotypes 
plot(geno[,which(p_values==max(p_values))],pheno$phenotype,
     xlab="genotypes",ylab="phenotypes",
     main=paste("effect size:",round(snp_effects[which(p_values==max(p_values))],2)))
mod=lm(pheno$phenotype~geno[,which(p_values==max(p_values))])
abline(mod,lwd=2,col="blue")

plot(density(-log10(p_values),na.rm=T),main="-log10(p-values)")
plot(density(snp_effects,na.rm=T),main="allelic effects")
plot(abs(snp_effects),-log10(p_values))         #plot estimated effect vs p-value of SNPs

# manhattan plot
#This uses the packages qqman,, you need to install the package an define the library
install.packages("qqman")   #you get a couple of messages, ignore them 
library(qqman)

map=read.table("map.txt", header=T, sep="\t")       # map file should have three columns with SNP CHR  BP as colnames. 
gwas = (cbind (map,p_values))                 # a file, with Chrom, Positions (bp) and p value
colnames(gwas)=c("SNP","CHR","BP","P")


jpeg('Manhattan.jpg',width=1100)      #name of output file for Manhattan Plot
manhattan(gwas,chr="CHR",ylim=c(0,10),suggestiveline=4,genomewideline=F,cex=1.2,cex.axis=1.1,col=c("blue","green"))
dev.off()

jpeg('QQ.jpg',width=1100) #name of output file for QQ  Plot (this plot sugegsts inflated p-values)
qq(p_values)    
dev.off()


# The next is to plot all markers in a manhattan plot, 
#   or markers of a particular chomosomes, or even parts of a single chromosome
# it doesn't have the nice color scheme as the manhattan() function

mbcol=as.factor(map$Chr)
chroms=unique(map$Chr)        

chromseps=numeric(length(chroms))

xdist=length(map[,3])
cum=0
chrpos=numeric(length(chroms))    #in the next loop we find the starting positions (cumulative genomewide) of the first SNP at each chrom, 
for (i in 1:length(chroms))
   {
     index=which(map[,2]==chroms[i])
     chromseps[i]=index[length(index)]
     xdist[index]=map[,3][index]+cum
     cum=cum+map[,3][index[length(index)]]
     chrpos[i]=cum
}

chromseps[0]=1
chrpos[2:length(chroms)]=
   chrpos[2:length(chroms)]-
   ((chrpos[2:length(chroms)]
       -chrpos[1:(length(chroms)-1)])/2)
chrpos[1]=chrpos[1]/2

par(mfrow = c(1,1))

# plot all markers
plot(xdist,-log10(p_values),pch=20,
     xlab="chromosome",
     ylab=expression(paste("-", log[10],"p-value", sep="")),
     axes=F,main="single SNP regressions")

abline(h=-log10(0.05/nrow(map)), lty=2)
axis(1, at=chrpos, labels=chroms,las=1)
axis(2,)

#make a plot for chrom i
i=3
if (i==1) 
{  index=1:chromseps[i]  } else
{ index=chromseps[i-1]:chromseps[i]}


plot(xdist[index],-log10(p_values[index]),pch=20,col="blue",xlab="map position",
     ylab=expression(paste("-", log[10],"p-value", sep="")),
     main=paste("chromosome",i), ylim=c(0, 7))
abline(h=-log10(0.05/nrow(map)), lty=2)
axis(1, at=chrpos, labels=chroms,las=1)




