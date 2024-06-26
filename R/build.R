#' Open a SNP GDS file and extract information.
#'
#' @param gds_data File name
#'
#' @param caco An object of class caco. Output of \code{case_control} function.
#' @param gwas_generator logical; if \code{TRUE} an object of class result_snps
#' will be saved inside output list.
#'
#' @return a list of seven objects; including smp.id, snp.id, snp.pos, smp.indx,
#' smp.snp (a matrix with samples in rows and snp in columns), caco,
#' snps(column sum of smp.snp for each case control)
#'
#'
#' @export
build_gwas <- function(gds_data = "name.gds", caco = "name.Rda", gwas_generator = TRUE) {
  if (missing(caco)) {
    stop(" provide case_control list object or case_control rda (contains list of case_control) file name ")
  }
  if (!is.list(caco)) {
    caco <- case_control(caco)
  }
  genoRA <- SNPRelate::snpgdsOpen(gds_data)
  smp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "sample.id"))
  snp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.rs.id"))
  snp.pos <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "snp.position"))
  smp.indx <- which(smp.id %in% unique(unlist(caco)))
  smp.snp <- list()
  for (j in 1:length(caco)) {
    smp.snp[[j]] <- (gdsfmt::read.gdsn(gdsfmt::index.gdsn(genoRA, "genotype"))[which(smp.id[smp.indx] %in% caco[[j]]), ])
    rownames(smp.snp[[j]]) <- smp.id[which(smp.id[smp.indx] %in% caco[[j]])]
    colnames(smp.snp[[j]]) <- snp.id
    smp.snp[[j]][smp.snp[[j]] == 3] <- NA
  }
  names(smp.snp) <- names(caco)
  SNPRelate::snpgdsClose(genoRA)
  output <- list(smp.id = smp.id, snp.id = snp.id, snp.pos = snp.pos, smp.indx = smp.indx, smp.snp = smp.snp, caco = caco)
  class(output) <- "gwas"
  if (gwas_generator) {
    output$snps <- extract(output)
  }
  return(output)
}


#' Read .vcf structured text format files and reduce the size of file.
#'
#' @param phased_vcf A file name for a variant call format (vcf) file.
#'
#' @param caco An object of class caco. Output of \code{case_control} function.
#'
#' @return the output will be a a list of class phase contains two sparse matrix
#' for each haplotype.
#'
#' @export
build_phase <- function(phased_vcf = "name.vcf", caco) {
  if (missing(caco) || is.null(phased_vcf)) stop("case_control and 'phased vcf' are needed")
  phased <- vector(mode = "list", length = 2)
  tmp <- data.table::fread(phased_vcf, nrows = 0)
  tmp2 <- data.table::fread(phased_vcf, select = which(colnames(tmp) %in% unique(unlist(caco))))
  ind1 <- which(tmp2[, lapply(.SD, substr, 1, 1)] == "1", arr.ind = TRUE)
  ind2 <- which(tmp2[, lapply(.SD, substr, 3, 3)] == "1", arr.ind = TRUE)
  phased[[1]] <- Matrix::sparseMatrix(i = ind1[, 1], j = ind1[, 2], x = 1, dims = c(nrow(tmp2), ncol(tmp2)))
  phased[[2]] <- Matrix::sparseMatrix(i = ind2[, 1], j = ind2[, 2], x = 1, dims = c(nrow(tmp2), ncol(tmp2)))
  colnames(phased[[1]]) <- colnames(phased[[2]]) <- colnames(tmp2)
  names(phased) <- c("Hap.1", "Hap.2")
  class(phased) <- "phase"
  return(phased)
}

#' Open a ibd file and extract information.
#'
#' @param ibd_data a file name for output of \href{http://faculty.washington.edu/browning/refined-ibd.html}{Refined IBD}
#'
#' @param gwas object of class gwas
#' @param gwid_generator logical; if \code{TRUE} an object of class result_snps
#' will be saved inside output list.
#'
#' @return the output will be a object(list) of class gwid contains
#' profile object, IBD object and result_snps object.
#'
#' @export
build_gwid <- function(ibd_data = "name.ibd", gwas = "object of class gwas", gwid_generator = TRUE) {
  ibd <- data.table::fread(ibd_data)
  V1 <- V2 <- V3 <- V4 <- V5 <- V6 <- V7 <- V8 <- V9 <- NULL
  ibd <- ibd[V1 %in% unlist(unique(gwas[["caco"]])) & V3 %in% unlist(unique(gwas[["caco"]]))]
  class(ibd) <- append("IBD", class(ibd))
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  profile <- ind <- vector(mode = "list", length = length(gwas[["caco"]])) # list length 6
  for (j in seq_along(gwas[["caco"]])) {
    ind[[j]] <- which(ibd$V1 %in% gwas[["caco"]][[j]] & ibd$V3 %in% gwas[["caco"]][[j]])
    a1 <- ibd[V1 %in% gwas[["caco"]][[j]] & V3 %in% gwas[["caco"]][[j]]]
    a2 <- seq2(match(a1$V6, gwas$snp.pos), match(a1$V7, gwas$snp.pos))
    Un1 <- unlist(a2)
    profile[[j]] <- Matrix::sparseMatrix(
      i = rep(seq_along(a2), lengths(a2)),
      j = Un1,
      x = 1
    )
    # if (ncol(profile[[j]]) < length(gwas$snp.pos) & sum(diff(sort(unique(Un1))))!=0 ){
    #   ind_not1 <- which(diff(sort(unique(Un1)))!=1)
    #   for (k in seq_along(ind_not1)) {
    #     mytemp <- Matrix(0, nrow = nrow(profile[[j]]), ncol = ( sort(unique(Un1))[ind_not1[k]+1] - sort(unique(Un1))[ind_not1[k]] - 1))
    #     profile[[j]] <- cbind(profile[[j]][,1:(ind_not1[[k]])] , mytemp , profile[[j]][,(ind_not1[[k]]+1):ncol(profile[[j]])])
    #   }
    # }

    if (ncol(profile[[j]]) < length(gwas$snp.pos)) {
      if (min(Un1) > 1) {
        mytemp <- Matrix(0, nrow = nrow(profile[[j]]), ncol = (min(Un1) - 1))
        profile[[j]] <- cbind(mytemp, profile[[j]][,min(Un1):max(Un1)])
      }

      if (max(Un1) < length(gwas$snp.pos)) {
        mytemp <- Matrix(0, nrow = nrow(profile[[j]]), ncol = (length(gwas$snp.pos) - max(Un1)))
        profile[[j]] <- cbind(profile[[j]], mytemp)
      }

      profile[[j]] <- methods::as(profile[[j]], "sparseMatrix")
    }
  }

  names(profile) <- names(ind) <- names(gwas[["caco"]])

  class(profile) <- "profile"
  output <- list(profile = profile, IND = ind)
  output$snp_pos <- gwas[["snp.pos"]]
  output$ibd <- ibd
  class(output) <- "gwid"
  if (gwid_generator) {
    output$res <- extract(output)
  }
  return(output)
}
