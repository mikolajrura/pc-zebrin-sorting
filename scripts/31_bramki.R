# Czystosc i wydajnosc bramek cytometrycznych dla sortowania jader Aldoc+/-.
# czystosc = jaki % jader w bramce to naprawde Aldoc+ ; wydajnosc = jaki % Aldoc+ zlapiemy
suppressMessages(library(Matrix))
a <- commandArgs(TRUE); tag <- a[1]; f <- a[2]; UMIMIN <- as.numeric(a[3])
x <- readRDS(gzcon(gzfile(f,"rb"))); rna <- attributes(x)$assays[[1]]
M <- attributes(rna)$counts; md <- attributes(x)$meta.data; rm(rna,x); invisible(gc())
cat(sprintf("\n================= %s =================\n", tag))
umi <- Matrix::colSums(M); cl <- as.character(md$seurat_clusters)
typ <- as.character(md$Type); an <- as.character(md$MULTIseq_group)
cnt <- function(g) if (g %in% rownames(M)) as.numeric(M[g,]) else rep(0,ncol(M))
PC <- c("Calb1","Car8","Pcp2","Ppp1r17","Itpr1","Slc1a6")
sg <- numeric(ncol(M)); for (g in PC) sg <- sg + cnt(g)/pmax(umi,1)*1e4
sb <- sort(tapply(sg,cl,median),decreasing=TRUE)
keep <- which(cl %in% names(sb)[sb>=sb[[1]]*0.35] & umi >= UMIMIN)
cp <- function(g) cnt(g)[keep]/umi[keep]*1e4
ty <- typ[keep]; am <- an[keep]
idx <- log2((cp("Aldoc")+1)/(cp("Plcb4")+1))
ap  <- idx > median(idx[ty=="WT"])            # prawda odniesienia, prog uczony na WT
cx  <- cp("Cux2") > 0; eb <- cp("Ebf2") > 0

bramki <- list(
  "Cux2+"            = cx,
  "Ebf2-"            = !eb,
  "Cux2+ i Ebf2-"    = cx & !eb,
  "Cux2+ lub Ebf2-"  = cx | !eb)

cat(sprintf("\n%-18s %6s %10s %10s %10s %10s\n","bramka","typ","%wbramce","czystosc","wydajnosc","wzbogac."))
for (nm in names(bramki)) {
  g <- bramki[[nm]]
  for (t in c("WT","SCA7")) {
    m <- ty==t; G <- g[m]; A <- ap[m]
    czyst <- mean(A[G]); wyd <- mean(G[A]); baza <- mean(A)
    cat(sprintf("%-18s %6s %9.1f%% %9.1f%% %9.1f%% %9.2fx\n",
        nm, t, 100*mean(G), 100*czyst, 100*wyd, czyst/baza))
  }
}

cat("\n--- czystosc bramki 'Cux2+ i Ebf2-' per zwierze ---\n")
g <- bramki[["Cux2+ i Ebf2-"]]; ans <- sort(unique(am))
CZ <- setNames(numeric(length(ans)), ans)
cat(sprintf("%-14s %5s %8s %10s %10s\n","zwierze","typ","%bramka","czystosc","wydajnosc"))
for (A in ans) { w <- am==A
  CZ[A] <- mean(ap[w & g])
  cat(sprintf("%-14s %5s %7.1f%% %9.1f%% %9.1f%%\n", A, unique(ty[w])[1],
      100*mean(g[w]), 100*CZ[A], 100*mean(g[w & ap]))) }
gt <- sapply(ans, function(A) unique(ty[am==A])[1])
W <- CZ[gt=="WT"]; S <- CZ[gt=="SCA7"]
cat(sprintf("\n  czystosc WT   : %s (sr %.3f)\n  czystosc SCA7 : %s (sr %.3f)\n  Wilcoxon p = %.4f   SCA7/WT = %.2f\n",
    paste(sprintf("%.3f",W),collapse=", "), mean(W),
    paste(sprintf("%.3f",S),collapse=", "), mean(S),
    wilcox.test(W,S)$p.value, mean(S)/mean(W)))
