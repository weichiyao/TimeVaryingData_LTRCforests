###=========================================================
### Complete data case, i.e. time-varying covariates
###=========================================================

library(survival)
make_realset <- function(x, naremove = TRUE){
  first <- with(pbcseq, c(TRUE, diff(id) !=0)) #first id for each subject
  last  <- c(first[-1], TRUE)  #last id
  time1 <- with(pbcseq, ifelse(first, 0, day))
  time2 <- with(pbcseq, ifelse(last,  futime, c(day[-1], 0)))
  event <- with(pbcseq, ifelse(last,  status, 0))
  event <- 1*(event==2)
  
  pbcseq$start <- time1
  pbcseq$stop <- time2
  pbcseq$event <-  event
  
  if (naremove == TRUE) {
    pbcseq$albumin<-log(pbcseq$albumin)
    pbcseq$bili<-log(pbcseq$bili)
    pbcseq$protime<-log(pbcseq$protime)
    
    ## only leave with the necessary covariates
    pbcseq <- pbcseq[,c("id","start","stop","event","age","edema","alk.phos","chol","ast","platelet","spiders","hepato","ascites","albumin","bili","protime")]
    names(pbcseq)[1:4] <- c("ID","Start","Stop","Event")
    
  
    naidx = which(rowSums(is.na(pbcseq))>0)
    naid_uniq = unique(pbcseq[naidx,]$ID)
    
    
    pbcseq0 = pbcseq
    for (ii in naid_uniq){
      ni = sum(pbcseq$ID == ii)
      ## find out which covariate has na values
      nacolidx_ii = which(colSums(is.na(pbcseq[pbcseq$ID == ii,]))>0)
      ## for each of those covariates that have na values
      for (kk in nacolidx_ii){
        ## find out in which line na values lie
        naidx_kk = which(is.na(pbcseq[pbcseq$ID == ii,kk])==1)
        ## find out in which line non-na values lie
        Nnaidx_kk = which(is.na(pbcseq[pbcseq$ID == ii,kk])==0)
        ## If there are at least one observed value
        if (length(Nnaidx_kk)>0){
          ## for each of those lines where na values lie
          for (ll in naidx_kk){
            ## for those na appear before the first non-na then go to the first non-na value that comes next
            if (ll < Nnaidx_kk[1]){
              pbcseq[pbcseq$ID == ii,kk][1] = pbcseq[pbcseq$ID == ii,kk][Nnaidx_kk[1]]
            } else {
              ## if not, fill in the blank with the closest one that comes before
              fillidx_ll = max(Nnaidx_kk[Nnaidx_kk<ll])
              # print(c(ii,kk,ll,fillidx_ll))
              pbcseq[pbcseq$ID == ii,kk][ll] = pbcseq[pbcseq$ID == ii,kk][fillidx_ll]
            }
          }
        }
      }
    }
    
    ## Get rid of those with all values in one column absent 
    pbcseq = pbcseq[rowSums(is.na(pbcseq)) == 0, ]
  } else {
    pbcseq <- subset(pbcseq, select = -c(futime, status, day))
    pbcseq$sex = ifelse(pbcseq$sex == "f", 1, 0)
    pbcseq$tdo = pbcseq$stop - pbcseq$start
  }
  
  return(pbcseq)
}  

# k=0
# for (i in 1:312){
#   if (sum(pbc$Event[pbc$ID == i])== 1){
#     k = k+1
#   }
# }
