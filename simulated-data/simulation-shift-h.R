source("methods/methods.R")
library(tidyverse)

out_rate <- .05
p <- 25
n <- 100
cl <- 0.25 #proportion of samples of each cluster (between 0 and 1)
t <- seq(0, 1, len = p)
d <- as.matrix(dist(t, upper = T, diag = T))

# Generate model
t.cov <- 0.3 * exp(-1 / 0.3 * d)  # covariance function in time
L <- chol(t.cov)  # Cholesky decomposition
mu <- 30 * t * (1 - t) ^ 2

methods = c(
            # "FOADA-TH",
            # "FOADA-RMAH",
            # "FADALARA",
            # "RMAH",
            # "ISFE",
            # "OUG",
            # "LRT",
            # "TRIM",
            # "POND",
            # "FOM",
            # "FB",
            # "HDR",
            # "AAKNN",
            "TRIMKMEANS"
            )
positives = c("TP", "FP")

nrep = 10
res <- array(0, dim = c(length(methods), nrep, 2))
dimnames(res)[[1]] <- methods
dimnames(res)[[3]] <- positives

clusters = c("single", "multiple")

for (cluster in clusters) {
    for (rep in 1:nrep) {
        cat(cluster, "-", rep, "\n")
        set.seed(rep)
        
        ### Generate Data
        if (cluster == "single") {
            e = matrix(rnorm(n * p), p, n)
            y = mu + t(L) %*% e
        } else {
            # Generate model
            e = matrix(rnorm(n * p * cl), p, n * cl)
            y1 = mu + t(L) %*% e
            
            e = matrix(rnorm(n * p * (1 - cl)), p, n * (1 - cl))
            t2 <- t + 0.3
            mu2 <- 30 * t2 * (1 - t2) ^ 2
            
            y2 = mu2 + t(L) %*% e
            
            y = cbind(y1, y2)
        }
        
        ### Generate Outliers
        out <- sort(sample(1:n, n * out_rate))
        for (i in out) {
            e = rnorm(p)
            t3 <- t + 0.15
            mu3 <- 30 * t3 * (1 - t3) ^ 2
            y[, i] = mu3 + t(L) %*% e
        }
        
        data = t(y)
        
        # out_partial <- OTRESH(data)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["FOADA-TH", rep, ] = res_partial
        # 
        # out_partial <- OMULTI(data)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["FOADA-RMAH", rep, ] = res_partial
        # 
        # out_partial <- FADALARA(data, k=3)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["FADALARA", rep, ] = res_partial
        # 
        # out_partial <- RMAH(data, t)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["RMAH", rep, ] = res_partial
        # 
        # out_partial <- ISE(data, t)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["ISFE", rep, ] = res_partial
        # 
        # 
        # out_partial <- OGH(data, t)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["OUG", rep, ] = res_partial
        #  
        # out_partial <- LRT(data, t, out_rate = out_rate)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["LRT", rep, ] = res_partial
        # 
        # out_partial <- TRIM(data, t, out_rate)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["TRIM", rep, ] = res_partial
        # 
        # out_partial <- POND(data, t, out_rate)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["POND", rep, ] = res_partial
        # 
        # out_partial <- FOM(data)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["FOM", rep, ] = res_partial
        # 
        # out_partial <- FB(data)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["FB", rep, ] = res_partial
        # 
        # out_partial <- HDR(data, t, out_rate)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["HDR", rep, ] = res_partial
        # 
        # out_partial <- AAKNN(data, vi = 5, vf = 15, narch = 6)
        # res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        # res["AAKNN", rep, ] = res_partial
        
        k <- if(cluster == "single") 1 else 2
        out_partial <- TRIMKMEANS(data, k, out_rate)
        res_partial = c(length(setdiff(out, out_partial)), length(setdiff(out_partial, out)))
        res["TRIMKMEANS", rep, ] = res_partial
    }
    
    tp = res[, , "TP", drop = FALSE]
    tpr = 1 - tp / length(out)
    
    tpr_mean = apply(tpr, 1, mean)
    tpr_sd = apply(tpr, 1, sd)
    
    fp = res[, , "FP", drop = FALSE]
    fpr = fp / (n - length(out))
    
    fpr_mean = apply(fpr, 1, mean)
    fpr_sd = apply(fpr, 1, sd)
    
    fpr_tibble <-
        tibble(
            method = row.names(fpr),
            mean = fpr_mean,
            sd = fpr_sd,
            rate = "false positive"
        )
    
    tpr_tibble <-
        tibble(
            method = row.names(tpr),
            mean = tpr_mean,
            sd = tpr_sd,
            rate = "true positive"
        )
    
    pr <-
        bind_rows(tpr_tibble, fpr_tibble) %>%
        add_column(cluster = cluster) %>%
        transform(out_rate=method %in% c("TRIM", "HDR", "LRT", "POND", "FADALARA", "TRIMKMEANS"))
    
    write_csv(pr, paste("simulated-data/results/csv/shift-h-", cluster, "-cluster-v2.csv", sep = ""))
}