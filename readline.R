library(stringr)
MSP_Path ="C:/Users/Public/專題研究/HMDB.txt"

tmp = readLines(MSP_Path)
### 哪裡出現關鍵字，例如Name
sIndex_Name = grep("NAME:", tmp) #正規表示式，應該像是擷取文字信息的某個特徵並標示在第幾行
sIndex_PrecursorMZ = grep("PRECURSORMZ: ", tmp)
sIndex_PrecursorTYPE = grep("PRECURSORTYPE: ", tmp)
sIndex_NumPeaks = grep("Num Peaks: ", tmp)

#把冒號前的字串丟掉
Name = gsub('.*:' , "", x = tmp[sIndex_Name])
print(Name)
PrecursorMZ = gsub(pattern = ".*: ", "", x = tmp[sIndex_PrecursorMZ])

NumPeaks = as.numeric(gsub("Num Peaks: " ,"" , tmp[sIndex_NumPeaks]))


S = sIndex_NumPeaks + 1#第幾行開始有Numpeaks數值
E = sIndex_NumPeaks + NumPeaks#Numpeaks數值在第幾行結束?
pIndex = cbind(S, E)#透過column合併資料

i = 1
Spec = list()
for(i in 1:nrow(pIndex))
{
  sIndex = pIndex[i, "S"] : pIndex[i, "E"]#pIndex為column資料pIndex[i, "S"]表示第i層的S
  
  #stringr::str_split(tmp[sIndex], "\t")
  mz = as.numeric(gsub("\t.*", "", tmp[sIndex]))#tmp = readLines(MSP_Path),
  intensity = as.numeric(gsub(".*\t", "", tmp[sIndex]))
  Spec[[i]] = cbind(mz, intensity)
}

  Spec[[i]] = cbind(mz, intensity)


Dot <- function(I1, I2)
{
  m=sum(I1*I2)
 return(m)
}

VecLen <- function(I1,I2)
{
  return(Dot(I1, I2) ^ 0.5)
}



#利用mapple()做內積?該如何把內積加總
DotSpec=list()
SumDotSpec=list()
square=list()
for(i in 1:nrow(pIndex)){
  DotSpec[[i]]=mapply(Dot,Spec[[i]][,2],Spec[[i]][,2])
  SumDotSpec[[i]]=mapply(sum,DotSpec[i])#?
  
}

Dot(Spec[[i]][,2], Spec[[i]][,2])

x = Spec[[i]][,2]
x = matrix(Spec[[i]][,2], ncol = 1)
x
t(x) %*% (x)


