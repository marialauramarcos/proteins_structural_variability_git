WindowsRMSDcontacts <- function(kij,
                                r.p.1,
                                r.p.2,
                                aligned.p.1.index.3N,
                                aligned.p.2.index.3N) {
  # select aligned sites          
  r.p.1 = r.p.1[aligned.p.1.index.3N]
  r.p.2 = r.p.2[aligned.p.2.index.3N] 
  
  # calculate the number of aligned sites
  n.aligned = length(r.p.1)/3
  
  # change kij
  diag(kij) = 1

  # build a vector to save the data
  score.windows = c()
  
  # sart a loop to aligne each window and calculate the score
  for (i in (1:n.aligned)) {
    
    indexes.window.i = which(kij[i, ] == 1)
    indexes.window.i.3N = sort(c(indexes.window.i * 3, indexes.window.i * 3 - 2,  indexes.window.i * 3 - 1))

    r.p.2.window = fit.xyz(fixed = r.p.1[indexes.window.i.3N],
                          mobile = r.p.2[indexes.window.i.3N],
                      fixed.inds = (1:length(indexes.window.i.3N)),
                     mobile.inds = (1:length(indexes.window.i.3N)))
    
    score.windows.i = mean(colSums((matrix((r.p.2.window - r.p.1[indexes.window.i.3N]), nrow = 3)) ^ 2)) 
    score.windows[i] = score.windows.i
  }
  as.vector(score.windows)
}
