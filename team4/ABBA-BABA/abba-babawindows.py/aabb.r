#This code is based on the excellent description here:
#https://www.math.wustl.edu/~sawyer/handouts/Jackknife.pdf
running the R comands and correct the conmands of jacknife.R to the following comands.

get.block.indices <- function(block_size, positions, chromosomes = NULL){
    if (is.null(chromosomes) == TRUE) {
        block_starts <- seq(min(positions), max(positions), block_size)
        
        block_ends <- block_starts + block_size - 1
        
        block_indices <- lapply(1:length(block_starts), function(x) which(positions >= block_starts[x] &
                                                                          positions <= block_ends[x]))
        }
    else {
        chrom_names <- unique(chromosomes)
        
        block_starts <- lapply(chrom_names, function(chrom_name) seq(min(positions[chromosomes==chrom_name]),
                                                                     max(positions[chromosomes==chrom_name]), block_size)) 
        
        block_chroms <- unlist(lapply(1:length(block_starts), function(x) rep(chrom_names[x], length(block_starts[[x]]))))
        
        block_starts <- unlist(block_starts)
        
        block_ends <- block_starts + block_size - 1
        
        block_indices <- lapply(1:length(block_starts), function(x) which(chromosomes == block_chroms[x] &
                                                                          positions >= block_starts[x] &
                                                                          positions <= block_ends[x]))
        }
    
    #identify blocks that are not empty
    not_empty <- which(sapply(block_indices, length) > 0)
    
    block_indices[not_empty]
    
    }

#this function runs the jackknife procedure by calculating pseudovalues by removing one block at a time
#if the arguments specified by "..." are vectors, they will be indexed as they are.
#if they have two dimensions, they will be indexed along the first dimension
block.jackknife <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))
        }
    else{
        pseudovalues <- sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))
        }
    
    mean <- mean(pseudovalues)
    
    v <- var(pseudovalues)
    
    std_dev <- sd(pseudovalues)
    
    se <- std_dev/sqrt(n_blocks)
    
    list(mean=mean, variance=v, standard_deviation=std_dev, standard_error=se)
}


get_genome_blocks <- function(block_size, chrom_lengths) {
  block_starts <- sapply(chrom_lengths, function(l) seq(1, l, block_size))
  data.frame(start = unlist(block_starts),
             end = unlist(block_starts) + block_size - 1,
             chrom = rep(names(block_starts), sapply(block_starts, length)))
}

get_genome_jackknife_indices <- function(chromosome, position, block_info){
  lapply(1:nrow(block_info), function(x) !(chromosome == block_info$chrom[x] &
                                             position >= block_info$start[x] &
                                             position <= block_info$end[x]))
}


get_jackknife_sd <- function(FUN, input_dataframe, jackknife_indices){
  n_blocks <- length(jackknife_indices)
  overall_mean <- FUN(input_dataframe)
  sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - FUN(input_dataframe[jackknife_indices[[i]],])*(n_blocks-1)))
}
