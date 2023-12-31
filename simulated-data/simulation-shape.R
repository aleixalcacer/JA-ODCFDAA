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
mu <- 30 * t * (1 - t) ^ (3 / 2)

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
        e = matrix(rnorm(n * p), p, n)
        y = mu + t(L) %*% e
        
        iso_interval = seq(-2, 2, .5)
        
        if (cluster == "multiple") {
            n_2 <- ceiling(n * cl)
            mu_2 <- 10 * t^(1/3) * (1 - t)^3
            e_2 = matrix(rnorm(n_2 * p), p, n_2)
            y_2 = mu_2 + t(L) %*% e_2
            y[, 1:n_2] = y_2
        }
        
        ### Generate Outliers
        out<-sample((n/2):n, n * out_rate)
        for (i in out){
            mu_3 = 15 * t^2 * (1-t)^2
            e_3 = rnorm(p)
            y[,i] =  mu_3 + t(L) %*% e_3
        }
        
        data = t(y)
        
        # matplot(y, type = "l")
        
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
        # out_partial <- AAKNN(data, vi = 5, vf = 10, narch = 6)
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
        transform(out_rate=method %in% c("TRIM", "HDR", "LRT", "POND", "TRIMKMEANS"))
    
    write_csv(pr, paste("simulated-data/results/csv/shape-", cluster, "-cluster-v2.csv", sep = ""))
}