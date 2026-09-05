# Czy Cux2 nadal nadaje sie na target sortowania jader u SCA7?
# (1) dwumodalnosc sygnalu Cux2 - pytanie o bramke, niezalezne od etykiet Aldoc
# (2) zgodnosc z osia Aldoc/Plcb4 - obciazone kolowoscia, patrz komentarz w raporcie
suppressMessages(library(Matrix))
a <- commandArgs(TRUE); tag <- a[1]; f <- a[2]
x <- readRDS(gzcon(gzfile(f,"rb"))); rna <- attributes(x)$assays[[1]]
M <- attributes(rna)$counts; md <- attributes(x)$meta.data; rm(rna,x); invisible(gc())
cat(sprintf("\n================= %s =================\n", tag))
umi <- Matrix::colSums(M); cl <- as.character(md$seurat_clusters)
typ <- as.character(md$Type); an <- as.character(md$MULTIseq_group)
has <- function(g) g %in% rownames(M)
cnt <- function(g) if (has(g)) as.numeric(M[g,]) else rep(0, ncol(M))

PC <- c("Calb1","Car8","Pcp2","Ppp1r17","Itpr1","Slc1a6"); PC <- PC[sapply(PC,has)]
sig <- numeric(ncol(M)); for (g in PC) sig <- sig + cnt(g)/pmax(umi,1)*1e4
sb <- sort(tapply(sig, cl, median), decreasing=TRUE)
sel <- cl %in% names(sb)[sb >= sb[[1]]*0.35]

T <- floor(quantile(umi[sel], 0.10)); keep <- which(sel & umi >= T); p <- T/umi[keep]
cat(sprintf("komorki Purkinjego %d, po odsianiu ponizej T=%d UMI: %d\n", sum(sel), T, length(keep)))

kand <- c("Cux2","Aldoc","Plcb4","Cpne9","Wscd2","Camk2d","Sorcs2","Mpped2","Kctd12","Corin","Lhx5","Ebf2")
cat("obecne w macierzy:", paste(kand[sapply(kand,has)], collapse=", "), "\n")
if (any(!sapply(kand,has))) cat("BRAK:", paste(kand[!sapply(kand,has)], collapse=", "), "\n")
RAW <- sapply(kand, function(g) cnt(g)[keep])

set.seed(1)
TH <- sapply(kand, function(g) rbinom(length(keep), RAW[,g], p))   # wyrownana glebokosc
colnames(TH) <- kand
ty <- typ[keep]; am <- an[keep]

cat("\n--- (1) rozklad sygnalu Cux2 po wyrownaniu glebokosci (T UMI na komorke) ---\n")
cat(sprintf("%-6s %6s %8s %8s %8s %8s %8s\n","typ","n","%zero","%1","%2","%>=3","med|>0"))
for (t in c("WT","SCA7")) {
  v <- TH[ty==t,"Cux2"]
  cat(sprintf("%-6s %6d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %8.1f\n", t, length(v),
      100*mean(v==0), 100*mean(v==1), 100*mean(v==2), 100*mean(v>=3),
      if (any(v>0)) median(v[v>0]) else NA))
}

cat("\n--- (2) os Aldoc/Plcb4 i zgodnosc Cux2, per zwierze ---\n")
idx <- log2((TH[,"Aldoc"]+1)/(TH[,"Plcb4"]+1))
cat(sprintf("%-14s %5s %6s %9s %9s %9s %9s\n","zwierze","typ","n","%Aldoc+","f_Cux2|A+","f_Cux2|A-","d"))
res <- list()
for (A in sort(unique(am))) {
  w <- am==A; ii <- which(w)
  ap <- idx[ii] > 0                                   # Aldoc przewaza nad Plcb4
  cp <- TH[ii,"Cux2"] > 0
  fp <- mean(cp[ap]); fn <- mean(cp[!ap]); d <- fp-fn
  res[[A]] <- c(d=d, fp=fp, fn=fn, ap=mean(ap), t=NA)
  cat(sprintf("%-14s %5s %6d %8.1f%% %9.3f %9.3f %9.3f\n", A, unique(ty[w])[1], sum(w),
      100*mean(ap), fp, fn, d))
}
gt <- sapply(sort(unique(am)), function(A) unique(ty[am==A])[1])
D  <- sapply(res, function(r) r["d"]); AP <- sapply(res, function(r) r["ap"])
for (nm in list(c("d","rozdzielczosc Cux2 (d)"), c("ap","udzial Aldoc+"))) {
  V <- if (nm[1]=="d") D else AP
  W <- V[gt=="WT"]; S <- V[gt=="SCA7"]
  cat(sprintf("\n  %s\n    WT   : %s  (srednia %.3f)\n    SCA7 : %s  (srednia %.3f)\n",
      nm[2], paste(sprintf("%.3f",W),collapse=", "), mean(W),
      paste(sprintf("%.3f",S),collapse=", "), mean(S)))
  cat(sprintf("    Wilcoxon p = %.4f   [min dla 4v4 = 0.0286]   spadek SCA7/WT = %.2f\n",
      wilcox.test(W,S)$p.value, mean(S)/mean(W)))
}
