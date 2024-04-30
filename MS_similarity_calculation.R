install.packages("msentropy")
install.packages('tidyverse')
library(tidyverse)
library(stringr)
library(matrixsimilarity)
library(msentropy)
library(purrr)
library(ggplot2)

#建立資料庫
MSP_Path ="C:/Users/Public/MSMSdotproductsimilarity/HMDB.txt"

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

#未知的compound
ms_unknown <- matrix(c(108.889,113.407,114.285,114.536,136.499,136.876,154.447,155.451,155.577, 15.431,53.522,100.0,47.916,22.35,17.753,11.029,32.825,24.099), nrow = 9, ncol = 2, byrow = TRUE,
               dimnames = list(NULL,c("mz", "intensity")))
print(ms_unknown)

#計算entropy_similarity

## seq_along(Spec)創建一個序列，其長度與 Spec 列表相同。此序列將用於循環遍歷 Spec
##~ {...} 部分定義了一個匿名函數，該函數將應用於 Spec 列表中的每個元素
##[.x] 是 R 中一種特殊的運算符，用於根據索引訪問嵌套列表或數據框中的元素。在這種情況下，它與 Spec[[.x]] 相同，但更簡潔地用於匿名函數
similarity_Score <- map(seq_along(Spec), ~ { 
  calculate_entropy_similarity(
    peaks_a = ms_unknown,
    peaks_b = Spec[[.x]],
    ms2_tolerance_in_da = 0.01,
    ms2_tolerance_in_ppm = 10,
    clean_spectra = FALSE,
    min_mz = 0,
    max_mz = 1000,
    noise_threshold = 0,
    max_peak_num = 10
  )
})

#畫質譜比對圖
for(i in 1:length(Spec)){
    png(paste0("~/Spec_similarity/", i,".png"), width  = 800, height = 800, units = "px" )  
    plot(ms_unknown[,1], ms_unknown[,2], type = "h",##h=histogram 
         xlab ="mz", ylab = "Intensity", ylim = c(-100, 100), lwd = 2,
         main = "Spectral Entropy Similarity")
    abline(h = 0)
    lines(Spec[[i]][,1],-Spec[[i]][,2], type ="h",
          col = 2, lwd = 2)
    mtext(paste0("Similarity Score:", similarity_Score[[i]]), side = 3, line = 1,adj =0)
    mtext(paste('black line: unknown','red line: reference',sep =' '),side = 1, line =1,adj =1,padj =1)
    dev.off()
  
}
