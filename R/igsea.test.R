igsea.test <-
  function(gel, pheno, ssize, gind, gsind, B = 500, vtype = "binary", method = "AT", alpha1 = 0.0253, pihat = 1) {
    if (dim(gel)[2] != sum(ssize) | dim(gel)[1] != dim(gind)[1] | length(ssize) != dim(gind)[2] | length(pheno) != sum(ssize) | dim(gel)[1] != dim(gsind)[1] | is.element(method, c("AT", "FE", "RE")) == "FALSE" | is.element(vtype, c("binary", "continuous")) == "FALSE") {
      stop ("Wrong settings!")
    } else { 
      gcount = dim(gind)[1]
      nstudy = dim(gind)[2]
      ngs = dim(gsind)[2]
      KS = function(Z, set_id){
        ind_s = which(gsind[, set_id] == 1)
        M = length(ind_s)
        hs = ks.test(Z[ind_s], Z[-ind_s], alternative = "greater")$statistic
        es = hs/sqrt(1/M + 1/(gcount - M))
        return(es)
      }
      if (vtype == "binary") {
        uvcal = function(X, Y){
          itcpt = glm(Y ~ 1, family = "binomial")$coefficients[1]
          ez = exp(itcpt)
          utemp = sum((Y - ez/(1 + ez)) * X)
          vtemp = sum(ez/(1 + ez)^2 * X^2)
          cb = complex(real = utemp, imaginary = vtemp)
          return(cb)
        }
        uv = matrix(0, gcount, nstudy)
        for (g in 1:gcount) {
          for (k in 1:nstudy) {
            if(gind[g, k] == 1) {
              uv[g, k] = uvcal(gel[g, (1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])], pheno[(1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])])
            }
          }
        }
        v = Im(uv)
        u = Re(uv)
        if (method == "FE") {
          Cfe = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cfe[g] = sum(utemp2)^2/sum(vtemp2)
          }
          ESfe = numeric(ngs)
          for (s in 1:ngs) ESfe[s] = KS(-Cfe, s)
          ESfe_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESfe_permu[s, ] = replicate(B, KS(sample(-Cfe, gcount), s))
          qfe = numeric(ngs)
          for (s in 1:ngs) qfe[s] = pihat * sum(ESfe_permu >= ESfe[s])/B/rank(-ESfe, ties.method = "max")[s]
          result = qfe
        } else if (method == "RE") {
          Cre = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cre[g] = sum(utemp2)^2/sum(vtemp2) + 0.5*(sum(utemp2^2) - sum(vtemp2))^2/sum(vtemp2^2)
          }
          ESre = numeric(ngs)
          for (s in 1:ngs) ESre[s] = KS(-Cre, s)
          ESre_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESre_permu[s, ] = replicate(B, KS(sample(-Cre, gcount), s))
          qre = numeric(ngs)
          for (s in 1:ngs) qre[s] = pihat * sum(ESre_permu >= ESre[s])/B/rank(-ESre, ties.method = "max")[s]
          result = qre
        } else if (method == "AT") {
          Cfe = numeric(gcount)
          heter = function(X, Y){
            glmtemp = glm(Y ~ X, family = "binomial")
            recordglm = complex(real = glmtemp$coefficients[2], imaginary = summary(glmtemp)$coefficients[4])
            return(recordglm)
          }
          betasigma = matrix(0, gcount, nstudy)
          for (g in 1:gcount) {
            for (k in 1:nstudy) {
              if(gind[g, k] == 1) {
                betasigma[g, k] = heter(gel[gcount, (1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])], pheno[(1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])])
              }
            }
          }
          beta = Re(betasigma)
          sd_beta = Im(betasigma)
          weight = 1/sd_beta^2
          mmmean = numeric(gcount)
          cochrQ = numeric(gcount)
          p_stage1 = numeric(gcount)
          p_stage2 = numeric(gcount)
          p_cb = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cfe[g] = sum(utemp2)^2/sum(vtemp2)
            p_stage2[g] = 1 - pchisq(Cfe[g], df = 1)
            if (length(indtemp) == 1) {
              p_cb[g] = p_stage2[g]
            } else {
              betatemp = beta[g, ][indtemp]
              weighttemp = weight[g, ][indtemp]
              mmmean[g] = sum(betatemp*weighttemp)/sum(weighttemp)
              cochrQ[g] = sum(weighttemp*(betatemp - mmmean[g])^2)
              p_stage1[g] = 1 - pchisq(cochrQ[g], df = length(indtemp) - 1)
              if (p_stage1[g] < alpha1) p_cb[g] = p_stage1[g]
              else p_cb[g] = alpha1 + p_stage2[g]*(1 - alpha1)
            } 
          }
          ESat = numeric(ngs)
          for (s in 1:ngs) ESat[s] = KS(p_cb, s)
          ESat_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESat_permu[s, ] = replicate(B, KS(sample(p_cb, gcount), s))
          qat = numeric(ngs)
          for (s in 1:ngs) qat[s] = pihat * sum(ESat_permu >= ESat[s])/B/rank(-ESat, ties.method = "max")[s]
          result = qat
        } else {
          result = "method can only be AT, FE or RE."
        }
      } else if (vtype == "continuous") {
        u = matrix(0, gcount, nstudy)
        v = matrix(0, gcount, nstudy)
        sigma = matrix(0, gcount, nstudy)
        beta = matrix(0, gcount, nstudy)
        sd_beta = matrix(0, gcount, nstudy)
        for (k in 1:nstudy) {                                            
          for (g in 1:gcount){
            Y = pheno[(1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])]
            X = gel[g, (1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])] 
            sigma[g, k] = sd(pheno[(1 + sum(ssize[0:(k - 1)])) : sum(ssize[0:k])])                                    
            temp = summary(lm(Y~X))                          
            beta[g, k] = temp$coefficients[2]                           
            sd_beta[g, k] = temp$coefficients[4]                       
            u[g, k] = (1 / sigma[g, k])^2 * sum((Y - mean(Y)) * X)     
            v[g, k] = (1 / sigma[g, k])^2 * sum(X * X) 
          }
        }
        if (method == "FE") {
          Cfe = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cfe[g] = sum(utemp2)^2/sum(vtemp2)
          }
          ESfe = numeric(ngs)
          for (s in 1:ngs) ESfe[s] = KS(-Cfe, s)
          ESfe_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESfe_permu[s, ] = replicate(B, KS(sample(-Cfe, gcount), s))
          qfe = numeric(ngs)
          for (s in 1:ngs) qfe[s] = pihat * sum(ESfe_permu >= ESfe[s])/B/rank(-ESfe, ties.method = "max")[s]
          result = qfe
        } else if (method == "RE") {
          Cre = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cre[g] = sum(utemp2)^2/sum(vtemp2) + 0.5*(sum(utemp2^2) - sum(vtemp2))^2/sum(vtemp2^2)
          }
          ESre = numeric(ngs)
          for (s in 1:ngs) ESre[s] = KS(-Cre, s)
          ESre_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESre_permu[s, ] = replicate(B, KS(sample(-Cre, gcount), s))
          qre = numeric(ngs)
          for (s in 1:ngs) qre[s] = pihat * sum(ESre_permu >= ESre[s])/B/rank(-ESre, ties.method = "max")[s]
          result = qre
        } else if (method == "AT") {
          Cfe = numeric(gcount)
          weight = 1/sd_beta^2
          mmmean = numeric(gcount)
          cochrQ = numeric(gcount)
          p_stage1 = numeric(gcount)
          p_stage2 = numeric(gcount)
          p_cb = numeric(gcount)
          for (g in 1:gcount) {
            indtemp = which(gind[g, ] > 0)
            utemp2 = u[g, ][indtemp]
            vtemp2 = v[g, ][indtemp]
            Cfe[g] = sum(utemp2)^2/sum(vtemp2)
            p_stage2[g] = 1 - pchisq(Cfe[g], df = 1)
            if (length(indtemp) == 1) {
              p_cb[g] = p_stage2[g]
            } else {
              betatemp = beta[g, ][indtemp]
              weighttemp = weight[g, ][indtemp]
              mmmean[g] = sum(betatemp*weighttemp)/sum(weighttemp)
              cochrQ[g] = sum(weighttemp*(betatemp - mmmean[g])^2)
              p_stage1[g] = 1 - pchisq(cochrQ[g], df = length(indtemp) - 1)
              if (p_stage1[g] < alpha1) p_cb[g] = p_stage1[g]
              else p_cb[g] = alpha1 + p_stage2[g]*(1 - alpha1)
            } 
          }
          ESat = numeric(ngs)
          for (s in 1:ngs) ESat[s] = KS(p_cb, s)
          ESat_permu = matrix(NA, ncol = B, nrow = ngs)
          for (s in 1:ngs) ESat_permu[s, ] = replicate(B, KS(sample(p_cb, gcount), s))
          qat = numeric(ngs)
          for (s in 1:ngs) qat[s] = pihat * sum(ESat_permu >= ESat[s])/B/rank(-ESat, ties.method = "max")[s]
          result = qat
        } else {
          result = "method can only be AT, FE or RE."
        }
      } else {
        result = "vtype can only be binary or continuous."
      }  
      return(result)
    }
  }
