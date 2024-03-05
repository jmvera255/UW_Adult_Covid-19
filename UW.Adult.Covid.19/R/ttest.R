#' Title
#'
#' @param data_in input data.frame
#' @param col1 vector of column names for 1st condition
#' @param col2 vector of column names for 2nd condition
#' @param col3 not used
#' @param label.col which column to use for rows
#' @param var.equal var.equal parameter for t.test call
#' @param alternative alternative parameter for t.test call
#' @param debug print debugging information
#'
#' @return data.frame of results
#' @export
#'
doTTest<-function(data_in, col1, col2, col3 = NA, label.col = "row.names", var.equal = TRUE, 
                  alternative="two.sided",debug=FALSE, is_paired=FALSE) {
  if (is_paired) {
    message("Running a paired t-test as requested. If unpaired test desired set is_paired = FALSE (default).")
  }

  if (!("ID" %in% colnames(data_in))) {
    if (label.col == "row.names") {
      data_in$ID = rownames(data_in);
    }   else if (label.col %in% colnames(data_in)) {
      data_in$ID = data_in[,label.col];
    } else {
      data_in$ID = 1:nrow(data_in);
    }
  }
  if(length(unique(data_in$ID)) != nrow(data_in)) {
    stop("label.col ",label.col, " or ID is not a unique col!");
  }

  group2=data.frame(data_in[,col2]);
  group1=data.frame(data_in[,col1]);
  base_mean        = rep(NA, nrow(data_in));
  mean1         = rep(NA, nrow(data_in));
  mean2         = rep(NA, nrow(data_in));
  fold_changes  = rep(NA, nrow(data_in));
  stderr1       = rep(NA, nrow(data_in));
  stderr2       = rep(NA, nrow(data_in));
  pvalues       = rep(NA, nrow(data_in));
  npvalues      = rep(NA, nrow(data_in));
  ttest_values  = rep(NA, nrow(data_in));
  base_stderr      = rep(NA, nrow(data_in));

  df = rep(NA, nrow(data_in))
  if (debug) {cat("group1:",dim(group1),"\n");}
  if (debug) {cat("group2:",dim(group2),"\n");}

  #group1/group2
  if (debug) {cat("iterating\n");}
  for (idx in 1:nrow(data_in)) {

    group1_gene=as.numeric(stats::na.omit(t(group1[idx,])));
    group2_gene=as.numeric(stats::na.omit(t(group2[idx,])));

    n1=length(group1_gene);
    n2=length(group2_gene);
    if(debug){cat("n1:",n1," n2:",n2,"\n");}
    if(debug){cat(dim(group1_gene),"\n");}
    if(debug){cat(dim(group2_gene),"\n");}
    current_base_mean = NA;
    current_mean1 = NA;
    current_mean2 = NA;

    current_foldchange = NA;
    current_pvalue = NA;
    current_ttest_value = NA;

    if (n1 > 0 || n2 > 0)
    {
      current_base_mean = mean(c(group1_gene, group2_gene));
    }
    if (n1 > 0) {
      current_mean1 = mean(group1_gene);
    }
    if (n2 > 0) {
      current_mean2 = mean(group2_gene);
    }
    if (n1 >= 1 && n2 >= 1 && (n1+n2) >= 3) {
      #print(group1_gene);
      #print(group2_gene);

      temp=try(stats::t.test(x=group1_gene, y=group2_gene, var.equal=var.equal, 
                             alternative=alternative, paired = is_paired));
      if ("try-error" %in% class(temp)) {
        current_foldchange = 0;
        current_pvalue = 1;
        current_ttest_value = 0;
      } else {
        # make compatible with paired test
        current_mean1 = mean(group1_gene);
        current_mean2 = mean(group2_gene);
        current_ttest_value = -temp$statistic;
        if (is_paired) {
          current_foldchange <- -1*temp$estimate[1] # mean of the differences
        } else {
          current_foldchange = as.numeric(temp$estimate[2]-temp$estimate[1]);
        }
        current_pvalue = as.numeric(temp$p.value);
        if (current_pvalue == "NaN") {
            current_pvalue = 1.0;
        }
        stderr1[idx] = sd(group1_gene);
        stderr2[idx] = sd(group2_gene);
        base_stderr[idx] = sd(c(group1_gene, group2_gene))
        df[idx] = temp$parameter;
      }
    } else if (n1 >= 1 && n2 >= 1) {
      if (is_paired) {
        current_foldchange <- -1*temp$estimate[1] # mean of the differences
      } else {
        current_foldchange= (current_mean2 - current_mean1);
      }
    }
    base_mean[idx] = current_base_mean;
    mean1[idx] = current_mean1;
    mean2[idx] = current_mean2;
    ttest_values[idx] = current_ttest_value;
    fold_changes[idx] = current_foldchange;
    pvalues[idx] = current_pvalue;
    npvalues[idx] = 1 - current_pvalue;

    if (idx %% 100000 == 0) {
      cat("processed ", idx, " out of ",nrow(data_in), "\n");
    }
  }

  #Calculate the adjusted p-values using the b-h algorithm
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("multtest")
  #
  library(multtest)
  res<-multtest::mt.rawp2adjp(pvalues,c("BH"),na.rm=TRUE);
  adjp=res$adjp[order(res$index),];

  res<-mt.rawp2adjp(npvalues,c("BH"),na.rm=TRUE);
  adjnp=res$adjp[order(res$index),];
  
  if (is_paired) {
    ans=data.frame(
    ID=data_in$ID,
    base_mean = base_mean,
    mean1 = mean1,
    mean2 = mean2,
    statistic = ttest_values,
    stderr1 = stderr1,
    stderr2 = stderr2,
    base_stderr = base_stderr,
    df = df,
    mean_of_diffs=2^fold_changes,
    log2_mean_of_diffs = fold_changes,
    pvalue=pvalues,
    npvalue = npvalues,
    padj=adjp[,2],
    npadj=adjnp[,2],
    stringsAsFactors=FALSE)
    ans$foldChange[is.na(ans$mean_of_diffs)] = 1.0;
  } else {
    ans=data.frame(
    ID=data_in$ID,
    base_mean = base_mean,
    mean1 = mean1,
    mean2 = mean2,
    statistic = ttest_values,
    stderr1 = stderr1,
    stderr2 = stderr2,
    base_stderr = base_stderr,
    df = df,
    foldChange=2^fold_changes,
    log2FC = fold_changes,
    pvalue=pvalues,
    npvalue = npvalues,
    padj=adjp[,2],
    npadj=adjnp[,2],
    stringsAsFactors=FALSE)
    ans$foldChange[is.na(ans$foldChange)] = 1.0;
  }
  rownames(ans) = ans$ID;
  #ans = ans[!is.na(ans$foldChange),];
  ans$pvalue[is.na(ans$pvalue)] = 1.0;
  ans$padj[is.na(ans$padj)] = 1.0;
  ans = ans[order(ans$padj,ans$pvalue),];
  cat("ans:",dim(ans),"\n");

  return(ans);
}

