
getMetricStats<-function(blocks, probes, metric, label, proteins=getProteinLabel(probes), pos = getProteinStart(probes)) {

  cat("getMetricStats - start\n");  
  min_value = rep(NA, nrow(blocks))
  max_value = rep(NA, nrow(blocks))
  mean_value = rep(NA, nrow(blocks))
  
  
  
  for (idx in 1:nrow(blocks)) {
    metric_values = metric[proteins == blocks$Protein[idx] & pos >= blocks$Start[idx] & pos <= blocks$Stop[idx]];
    min_value[idx] = min(metric_values, na.rm=TRUE);
    max_value[idx] = max(metric_values, na.rm=TRUE);
    mean_value[idx] = mean(metric_values, na.rm=TRUE);    
  }
  
  ans = data.frame(
    min_value = min_value,
    max_value = max_value,
    mean_value = mean_value
  );
  colnames(ans) = c(
    paste0("Min.",label),
    paste0("Max.",label),
    paste0("Mean.",label)
  )
  return(ans);
}


addSequenceAnnotations<-function(blocks_df, probe_meta) {
  
  ans_df = blocks_df
  umeta = unique(probe_meta[,c("PROBE_ID","PROBE_SEQUENCE")]);
  rownames(umeta) = umeta$PROBE_ID;
  
    
  
  first_probe = paste0(ans_df$Protein,";",ans_df$Start);
  last_probe = paste0(ans_df$Protein,";",ans_df$Stop);
  ans_df$First.Sequence = umeta[first_probe,"PROBE_SEQUENCE"];
  ans_df$Last.Sequence = umeta[last_probe,"PROBE_SEQUENCE"];
  ans_df$Overlap.Sequence = NA;
  ans_df$Full.Sequence = NA;
  nprobes = ans_df$NPeptides;
  for (idx in 1:nrow(ans_df)) {
    start = nprobes[idx];
    stop = nchar(ans_df$First.Sequence[idx])
    ans_df$Overlap.Sequence[idx] = substr(ans_df$First.Sequence[idx], start, stop);
    probes = paste0(ans_df$Protein[idx], ";",ans_df$Start[idx]:ans_df$Stop[idx]);
    
      if (nprobes[idx] == 1) { #If just one probe, full sequence is just the probe.
        ans_df$Full.Sequence[idx] = ans_df$First.Sequence[idx];
      } else {
        #If first and last sequence overlap, then just cat those two.  Else, stitch the whole thing together.
        if (start <= stop) {
          ans_df$Full.Sequence[idx] = catSequences(c(ans_df$Start[idx], ans_df$Stop[idx]),c(ans_df$First.Sequence[idx],ans_df$Last.Sequence[idx]))
        } else {      
          ans_df$Full.Sequence[idx] = catSequences(ans_df$Start[idx]:ans_df$Stop[idx],umeta[probes,"PROBE_SEQUENCE"]);
        }
      }
#      if (idx %% 1000 == 0) {
#        cat(idx, " of ",nrow(ans_df),"\n");
#      }
    }
  return(ans_df);
}


findBlocks<-function(protein.df) {
  
  if (nrow(protein.df)  == 1) {
    return(
      data.frame(
        Protein = protein.df$Protein,
        Start = protein.df$Pos,
        Stop = protein.df$Pos,
        Min.Pvalue = protein.df$Pvalue,
        Max.Pvalue = protein.df$Pvalue,
        Mean.Pvalue = protein.df$Pvalue,
        Min.logFC = protein.df$logFC,
        Max.logFC = protein.df$logFC,
        Mean.logFC = protein.df$logFC,
        Min.Signal = protein.df$Mean.Signal,
        Max.Signal = protein.df$Mean.Signal,
        Mean.Signal = protein.df$Mean.Signal,
        stringsAsFactors = FALSE
      )
    );
  }
  
  protein.df = protein.df[order(protein.df$Pos, decreasing=FALSE),];
  ans.df = NULL;
  start_idx=1;
  start_pos = protein.df$Pos[start_idx];
  for (idx in 2:nrow(protein.df)) {
    current_idx = idx;
    current_pos = protein.df$Pos[current_idx];
    prev_idx = idx-1;
    prev_pos = protein.df$Pos[prev_idx];
    d = current_pos - prev_pos;
    
    if (d > 1) {
      #Block is start_pos to prev_pos
      #cat("Add ",start_idx," ",prev_idx,"\n");
      pvalues = protein.df$Pvalue[start_idx:prev_idx];
      folds = protein.df$logFC[start_idx:prev_idx];
      signals = protein.df$Mean.Signal[start_idx:prev_idx];
      #print(pvalues)
      ans.df = rbind(
        ans.df,
        data.frame(
          Protein=protein.df$Protein[1],
          Start = start_pos,
          Stop = prev_pos,
          Min.Pvalue = min(pvalues),
          Max.Pvalue = max(pvalues),
          Mean.Pvalue = mean(pvalues),
          Min.logFC = min(folds),
          Max.logFC = max(folds),
          Mean.logFC = mean(folds),
          Min.Signal = min(signals),
          Max.Signal = max(signals),
          Mean.Signal = mean(signals),
          stringsAsFactors = FALSE
        )
      )
      #Update start_pos
      start_idx = current_idx;
      start_pos = current_pos;
    }
  }
  if (start_idx != nrow(protein.df)) {
    #Add in last block
    #cat("Add last ",start_idx, " ", nrow(protein.df),"\n");
    pvalues = protein.df$Pvalue[start_idx:nrow(protein.df)]
    folds = protein.df$logFC[start_idx:prev_idx];
    signals = protein.df$Mean.Signal[start_idx:prev_idx];
    #print(pvalues)
    ans.df = rbind(
      ans.df,
      data.frame(
        Protein=protein.df$Protein[1],
        Start = start_pos,
        Stop = protein.df$Pos[nrow(protein.df)],
        Min.Pvalue = min(pvalues),
        Max.Pvalue = max(pvalues),
        Mean.Pvalue = mean(pvalues),
        Min.logFC = min(folds),
        Max.logFC = max(folds),
        Mean.logFC = mean(folds),
        Min.Signal = min(signals),
        Max.Signal = max(signals),
        Mean.Signal = mean(signals),
        stringsAsFactors = FALSE
      )
    )
  }
  
  
  return(ans.df);
  
}

catSequences<-function(positions, sequences,debug=FALSE) {

  if (length(sequences) == 1) {
    return(sequences);
  }
  positions = positions-positions[1]+1;
  
  seq = rep("",1000);
  for (idx in 1:length(positions)) {
    pos = positions[idx];
    aa_vec = strsplit(sequences[idx],"")[[1]]
    startp = pos;
    endp = pos + length(aa_vec) - 1
    seq[startp:endp] = aa_vec;    
  }
  
  seq = seq[seq!=""]
  seqc = paste0(seq, collapse="",sep="");
  return(seqc);
}

findEpitopes<-function(probes, logfc, pvalues, mean.signal, probe_meta, lfc.threshold = 0, pvalue.threshold = 0.05, two.sided = FALSE) {
  
  
  sig = pvalues < pvalue.threshold & logfc > lfc.threshold
  
  
  protein.df = data.frame(
    Probe = probes[sig],
    Protein = getProteinLabel(probes[sig]),
    Pos = getProteinStart(probes[sig]),
    logFC = logfc[sig],
    Pvalue = pvalues[sig],
    Mean.Signal = mean.signal[sig],
    stringsAsFactors=FALSE
  );
  
  protein_list = split(protein.df, protein.df$Protein);
  #print(names(protein_list))
  ans_list = lapply(protein_list, findBlocks)
  
  ans_dt = data.table::rbindlist(ans_list);
  ans_df = as.data.frame(ans_dt, stringsAsFactors=FALSE);
  
  if (!missing(probe_meta)) {
  
    umeta = unique(probe_meta[,c("PROBE_ID","PROBE_SEQUENCE")]);
    rownames(umeta) = umeta$PROBE_ID;
  
    
  
    first_probe = paste0(ans_df$Protein,";",ans_df$Start);
    last_probe = paste0(ans_df$Protein,";",ans_df$Stop);
    ans_df$First.Sequence = umeta[first_probe,"PROBE_SEQUENCE"];
    ans_df$Last.Sequence = umeta[last_probe,"PROBE_SEQUENCE"];
    ans_df$Overlap.Sequence = NA;
    ans_df$Full.Sequence = NA;
    nprobes = ans_df$Stop - ans_df$Start + 1;
    
    for (idx in 1:nrow(ans_df)) {
      start = nprobes[idx];
      stop = nchar(ans_df$First.Sequence[idx])
      ans_df$Overlap.Sequence[idx] = substr(ans_df$First.Sequence[idx], start, stop);
      probes = paste0(ans_df$Protein[idx], ";",ans_df$Start[idx]:ans_df$Stop[idx]);
      
      if (nprobes[idx] == 1) { #If just one probe, full sequence is just the probe.
        ans_df$Full.Sequence[idx] = ans_df$First.Sequence[idx];
      } else {
        #If first and last sequence overlap, then just cat those two.  Else, stitch the whole thing together.
        if (start <= stop) {
          ans_df$Full.Sequence[idx] = catSequences(c(ans_df$Start[idx], ans_df$Stop[idx]),c(ans_df$First.Sequence[idx],ans_df$Last.Sequence[idx]))
        } else {      
          ans_df$Full.Sequence[idx] = catSequences(ans_df$Start[idx]:ans_df$Stop[idx],umeta[probes,"PROBE_SEQUENCE"]);
        }
      }
#      if (idx %% 1000 == 0) {
#        cat(idx, " of ",nrow(ans_df),"\n");
#      }
    }

  }
  
  if (two.sided) {
    cat("Getting other side\n");
    ans_df_2 = findEpitopes(
      probes = probes,
      logfc = -logfc,
      pvalues = pvalues,
      mean.signal = mean.signal,
      probe_meta = probe_meta,
      lfc.threshold = lfc.threshold,
      pvalue.threshold = pvalue.threshold,
      two.sided = FALSE
    );
    if (nrow(ans_df_2) > 0) {
      ans_df_2$Max.logFC=-ans_df_2$Max.logFC;
      ans_df_2$Min.logFC=-ans_df_2$Min.logFC
      ans_df_2$Mean.logFC = -ans_df_2$Mean.logFC;
      ans_df = rbind(ans_df, ans_df_2);
    }
  }
  if (nrow(ans_df) > 0) {
    ans_df = ans_df[order(ans_df$Min.Pvalue),]
  }
  attr(ans_df,"protein.df") = protein.df;
  
  
  return(ans_df);
}

getEpitopeSequence<-function(protein, start, stop, probe_meta, umeta = unique(probe_meta[,c("PROBE_ID","PROBE_SEQUENCE")])) {

  rownames(umeta) = umeta$PROBE_ID;

  ans = rep(NA, length(protein))
  for (idx in 1:length(protein)) {
    probes = paste0(protein[idx],";",start[idx]:stop[idx])
    if (length(probes) == 1) {
      ans[idx] = umeta[probes,"PROBE_SEQUENCE"]
    } else {
      ans[idx] = catSequences(start[idx]:stop[idx],umeta[probes,"PROBE_SEQUENCE"])
    }
  }
  return(ans);
}


findEpitopesTTestProbe<-function(ttest_result_probe, probe_meta, lfc.threshold, pvalue.threshold) {
  ans = findEpitopes(
    rownames(ttest_result_probe),
    ttest_result_probe$log2FC,
    ttest_result_probe$padj,
    ttest_result_probe$mean12,
    lfc.threshold = lfc.threshold,
    pvalue.threshold = pvalue.threshold,
    probe_meta = probe_meta
  );
  
}


getEpitopeMat<-function(probe_mat, epitopes_df) {

  Epitope_Mat = NULL;

  Epitope_Mat_list = list();

  for (epitope_idx in 1:nrow(epitopes_df)) {
    probes = paste0(
      epitopes_df$Protein[epitope_idx],";",
      epitopes_df$Start[epitope_idx]:epitopes_df$Stop[epitope_idx]
    );

    epitope_id = epitopes_df$Epitope_ID[epitope_idx];
    values = probe_mat[probes,]
    if (length(probes) > 1) {
      values = colMeans(values);
    }
    Epitope_Mat_list[[epitope_id]] = data.frame(t(values), check.names=FALSE);

  }
  library(data.table)
  Epitope_Mat_dt = rbindlist(Epitope_Mat_list);
  Epitope_Mat = as.data.frame(Epitope_Mat_dt, stringsAsFactors=FALSE);
  rownames(Epitope_Mat) = names(Epitope_Mat_list);
  
  return (Epitope_Mat);
  
}

