# Walidacja kandydatow JADROWYCH (z 29_lokalizacja.py) na danych SCA7.
# Uzycie: Rscript 30_kandydaci_jadrowi.R "<tag>" <plik.rds.gz> <UMImin> "<gen1,gen2,...>"
suppressMessages(library(Matrix))
a <- commandArgs(TRUE); tag <- a[1]; f <- a[2]; UMIMIN <- as.numeric(a[3])
x <- readRDS(gzcon(gzfile(f,"rb"))); rna <- attributes(x)$assays[[1]]
M <- attributes(rna)$counts; md <- attributes(x)$meta.data; rm(rna,x); invisible(gc())
cat(sprintf("\n================= %s =================\n", tag))
umi <- Matrix::colSums(M); cl <- as.character(md$seurat_clusters)
typ <- as.character(md$Type); an <- as.character(md$MULTIseq_group)
cnt <- function(g) if (g %in% rownames(M)) as.numeric(M[g,]) else rep(0,ncol(M))
PC <- c("Calb1","Car8","Pcp2","Ppp1r17","Itpr1","Slc1a6")
sig <- numeric(ncol(M)); for (g in PC) sig <- sig + cnt(g)/pmax(umi,1)*1e4
sb <- sort(tapply(sig,cl,median),decreasing=TRUE)
keep <- which(cl %in% names(sb)[sb>=sb[[1]]*0.35] & umi >= UMIMIN)
cat(sprintf("komorki Purkinjego powyzej %d UMI: %d  (medUMI WT %.0f / SCA7 %.0f)\n",
    UMIMIN, length(keep), median(umi[keep][typ[keep]=="WT"]), median(umi[keep][typ[keep]=="SCA7"])))

cp <- function(g) cnt(g)[keep]/umi[keep]*1e4
idx <- log2((cp("Aldoc")+1)/(cp("Plcb4")+1))
ty <- typ[keep]; am <- an[keep]
THR <- median(idx[ty=="WT"])            # prog uczony WYLACZNIE na WT -> 50% Aldoc+ u WT
cat(sprintf("prog osi uczony na WT (mediana indeksu) = %.4f  -> u WT z definicji 50%%\n", THR))
ap <- idx > THR
cat(sprintf("udzial Aldoc+ : WT %.1f%%   SCA7 %.1f%%\n", 100*mean(ap[ty=="WT"]), 100*mean(ap[ty=="SCA7"])))

cat("\n--- rozdzielczosc kandydatow (d = f_pos|Aldoc+ minus f_pos|Aldoc-), przy pelnej glebokosci ---\n")
kand <- strsplit(a[4], ",")[[1]]
cat(sprintf("%-8s %10s %10s %8s | %10s %10s %8s   %s\n","gen","WT f|A+","WT f|A-","WT d","SC f|A+","SC f|A-","SC d","d_SCA7/d_WT"))
for (g in kand) {
  v <- cp(g) > 0; r <- c()
  for (t in c("WT","SCA7")) { m <- ty==t
    r <- c(r, mean(v[m & ap]), mean(v[m & !ap])) }
  dW <- r[1]-r[2]; dS <- r[3]-r[4]
  cat(sprintf("%-8s %10.3f %10.3f %8.3f | %10.3f %10.3f %8.3f   %.2f\n", g, r[1],r[2],dW, r[3],r[4],dS, dS/dW))
}

cat("\n--- per zwierze: udzial Aldoc+ i d dla Cux2 ---\n")
cat(sprintf("%-14s %5s %6s %9s %8s\n","zwierze","typ","n","%Aldoc+","d_Cux2"))
cx <- cp("Cux2") > 0; ans <- sort(unique(am)); D <- AP <- setNames(numeric(length(ans)),ans)
for (A in ans) { w <- am==A
  AP[A] <- mean(ap[w]); D[A] <- mean(cx[w & ap]) - mean(cx[w & !ap])
  cat(sprintf("%-14s %5s %6d %8.1f%% %8.3f\n", A, unique(ty[w])[1], sum(w), 100*AP[A], D[A])) }
gt <- sapply(ans, function(A) unique(ty[am==A])[1])
for (nm in c("AP","D")) { V <- get(nm); W <- V[gt=="WT"]; S <- V[gt=="SCA7"]
  cat(sprintf("\n  %s\n    WT   : %s (sr %.3f)\n    SCA7 : %s (sr %.3f)\n    Wilcoxon p = %.4f   SCA7/WT = %.2f\n",
      ifelse(nm=="AP","udzial Aldoc+","rozdzielczosc Cux2 (d)"),
      paste(sprintf("%.3f",W),collapse=", "), mean(W),
      paste(sprintf("%.3f",S),collapse=", "), mean(S),
      wilcox.test(W,S)$p.value, mean(S)/mean(W))) }
