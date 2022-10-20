
#' Load the probe metadata
#'
#' @param file_name input filename
#'
#' @return data.frame of meta data from probes
#' @export
#'
loadProbeMeta<-function(file_name="all_sequences_except_wi.tsv.gz") {

  probe_meta = read.table(file_name, sep="\t",header=TRUE);

  probe_meta$PROBE_ID = paste0(probe_meta$SEQ_ID,";",probe_meta$POSITION);
  return(probe_meta);
}


#' Load the sequence matrix
#'
#' @param file_name input filename
#' @param col_name signal column to use (Default INTENSITY)
#'
#' @return
#' @export
#'
loadSeqMat<-function(file_name="df_stacked.tsv", col_name="INTENSITY") {
  IgG_stacked = read.table(file_name,sep="\t",header=TRUE);
  #cat("Generating ustack\n");
  ustack = unique(IgG_stacked[,c("SAMPLE_NAME","PROBE_SEQUENCE",col_name)]);

  colnames(ustack)[3] = "lSignal";

  sample_ids = unique(ustack$SAMPLE_NAME);
  sample_df = data.frame();

  for (sample_id in sample_ids) {
  #cat("sample_id:",sample_id,"\n");
    idx = ustack$SAMPLE_NAME == sample_id;
    sub.df = data.frame(
      seq = ustack$PROBE_SEQUENCE[idx],
      lSignal = ustack[idx,3]
    );
    if (nrow(sample_df) == 0) {
      sample_df = sub.df
    }   else {
      sample_df = merge(sample_df, sub.df, by.x="seq", by.y="seq");
    }
    colnames(sample_df)[ncol(sample_df)] = paste0("Sample_",sample_id);
  }
  rownames(sample_df) = sample_df$seq;

  sample_meta = unique(IgG_stacked[,c("SAMPLE_NAME","COVID_POSITIVE")]);
  sample_meta$SAMPLE_NAME = paste0("Sample_",sample_meta$SAMPLE_NAME);

  attr(sample_df, "sample_meta") = sample_meta


  return(sample_df);
}


getProteinLabel<-function(probe.ids) {

  probe.list = strsplit(probe.ids,";");
  protein.list = lapply(
    probe.list, function(x){
      n = length(x)-1;
      return(paste(x[1:n],sep=";",collapse=";"))
    }
  )
  protein.labels = unlist(protein.list);
  return(protein.labels);

}

getProteinStart<-function(probe.ids) {
  probe.list = strsplit(probe.ids, ";");
  protein.start.list = lapply(
    probe.list,
    function(x) {

      return(x[length(x)]);
    }
  );
  return(as.integer(unlist(protein.start.list)));
}


getEpitopeProtein<-function(epitope_id) {
  ans = unlist(
    lapply(
      strsplit(epitope_id,"_"),
      function(l) {
        ans = l[1:(length(l)-2)]
        return(paste(ans,collapse="_"))
      }
    )
  )
  return(ans);
}

getEpitopeStart<-function(epitope_id) {
  ans = unlist(
    lapply(
      strsplit(epitope_id,"_"),
      function(l) {
        return(l[length(l)-1])
      }
    )
  )
  return(as.integer(ans));
}

getEpitopeStop<-function(epitope_id) {
  ans = unlist(
    lapply(
      strsplit(epitope_id,"_"),
      function(l) {
        return(l[length(l)])
      }
    )
  )
  return(as.integer(ans));
}

getEpitopeProbeIDs<-function(epitope_id) {
  protein = getEpitopeProtein(epitope_id);
  start = getEpitopeStart(epitope_id)
  stop = getEpitopeStop(epitope_id);

  ans = paste0(protein,";",start:stop);
  return(ans);
}

#' Make Epitope Id
#'
#' @param protein protein of the epitope
#' @param start first probe start position
#' @param stop last probe start position
#'
#' @return epitope identifer string
#' @export
#'
#' @examples
#' getEpitopeID("A", 1, 3)
getEpitopeID<-function(protein, start, stop) {
  return(paste(protein,start,stop, sep="_"))
}


#' Title
#'
#' @param in_mat
#'
#' @return
#' @export
#'
#' @examples
normalizeQuantile<-function(in_mat) {
  norm_mat = preprocessCore::normalize.quantiles(as.matrix(in_mat));
  norm_mat = as.data.frame(norm_mat);
  colnames(norm_mat) = colnames(in_mat);
  rownames(norm_mat) = rownames(in_mat);
  return(norm_mat);
}

getSequenceMatToProbeMat<-function(probe_meta, file_path) {
  if (!missing(file_path) && file.exists(file_path)) {
    cat("Loading seq_to_probe from ",file_path,"\n")
    load(file_path);
  } else {
    cat("Generating seq_to_probe\n")
    library(Matrix);

    umeta = unique(probe_meta[,c("PROBE_ID","PROBE_SEQUENCE")])


    uniq_seq = unique(umeta$PROBE_SEQUENCE);
    probe_idx = 1:nrow(umeta);
    seq_idx = 1:length(uniq_seq)
    names(seq_idx) = uniq_seq;

    seq_idx2 = seq_idx[umeta$PROBE_SEQUENCE]

    seq_to_probe = sparseMatrix(probe_idx, seq_idx2, x =1);
    rownames(seq_to_probe) = umeta$PROBE_ID;
    colnames(seq_to_probe) = uniq_seq;

    if (!missing(file_path)) {
      cat("Saving to ",file_path,"\n");
      save(list=c("seq_to_probe"), file=file_path)
    }
  }
  return(seq_to_probe);
}



#' Title
#'
#' @param seq_mat
#' @param probe_meta
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
convertSequenceMatToProbeMat<-function(seq_mat, probe_meta, file_path) {
  library(Matrix);
  seq_to_probe = getSequenceMatToProbeMat(probe_meta, file_path);
  probe_mat = as.matrix(seq_to_probe[,rownames(seq_mat)] %*% as.matrix(seq_mat));
  return(probe_mat)
}
