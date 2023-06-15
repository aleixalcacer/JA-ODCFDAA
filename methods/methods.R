



source('methods/aux/LFunctions_to_calculate_real_archetypes_with_swap2.R')
source('methods/aux/LstepArchetypesMod.R')
source('methods/aux/LstepArchetypoids.R')
source('methods/aux/robout.R')
source('methods/aux/fhdr.R')

library(archetypes)
library(Anthropometry)
library(roahd)
library(rainbow)
library(mrfDepth)
library(fda)
library(ks)
library(cluster)
library(fda.usc)
library(DDoutlier)


OTRESH <- function(data, th = 0.5) {
    y = t(data)
    X = t(y)
    set.seed(1234)
    lass10d <- stepLArchetypoids3(data = X, k = 3, 10)
    
    a3 = t(y[, lass10d[[1]][[1]]])
    data = t(y)
    huge = 200
    k = 3
    
    n <- ncol(t(data))
    x_gvv <- rbind(t(data), rep(huge, n))
    zs <- rbind(t(a3), rep(huge, k))
    zs <- as.matrix(zs)
    alphascce <- matrix(0, nrow = k, ncol = n)
    
    for (j in 1:n) {
        alphascce[, j] = coef(nnls(zs, x_gvv[, j]))
    }
    
    wm = which.min(c(
        sum(alphascce[1, ] > th[1]),
        sum(alphascce[2, ] > th[1]),
        sum(alphascce[3, ] > th[1])
    ))
    out = which(alphascce[wm, ] > th[1])
    return(out)
}


OMULTI <- function(data) {
    y = t(data)
    X = t(y)
    set.seed(1234)
    lass10d <- stepLArchetypoids3(data = X, k = 3, 10)
    
    a3 = t(y[, lass10d[[1]][[1]]])
    data = t(y)
    huge = 200
    k = 3
    
    n <- ncol(t(data))
    x_gvv <- rbind(t(data), rep(huge, n))
    zs <- rbind(t(a3), rep(huge, k))
    zs <- as.matrix(zs)
    alphascce <- matrix(0, nrow = k, ncol = n)
    
    for (j in 1:n) {
        alphascce[, j] = coef(nnls(zs, x_gvv[, j]))
    }
    
    sco = t(alphascce)
    rownames(sco) = as.numeric(1:n)
    s = cbind(sco, rep(1, n))
    
    out = robout(s, 1, "mcd")$outliers
    return(out)
}


OGH <- function(data, t) {
    fD = fData(t, data)
    out = outliergram(fD, display = F)$ID_outliers
    return(out)
}

RMAH <- function(data, t) {
    yf = fds(t, t(data))
    n <- ncol(t(data))
    colnames(yf$y) = 1:n
    out = foutliers(yf, method = "robMah")$outliers
    return(out)
}

ISE <- function(data, t) {
    yf = fds(t, t(data))
    n <- ncol(t(data))
    colnames(yf$y) = 1:n
    out = foutliers(
        data = yf,
        method = "HUoutliers",
        order = 2,
        lambda = 3.29
    )$outliers
    return(out)
}


LRT <- function(data, t, out_rate = 0.05) {
    Datosf <- fdata(data, t)
    out <-
        outliers.lrt(
            Datosf,
            nb = 200,
            smo = 0.05,
            dfunc = depth.mode,
            trim = out_rate
        )$outliers
    return(out)
}

TRIM <- function(data, t, out_rate = 0.05) {
    # Depth Based Trimming (Febrero, Galeano, Gonzalez-Manteiga, 2008)
    
    Datosf <- fdata(data, t)
    out <-
        outliers.depth.trim(
            Datosf,
            nb = 200,
            smo = 0.05,
            trim = out_rate,
            dfunc = depth.mode
        )$outliers
    return(out)
}

POND <- function(data, t) {
    # Depth Based Weighting (Febrero, Galeano, Gonzalez-Manteiga, 2008)
    Datosf <- fdata(data, t)
    out <-
        outliers.depth.pond(Datosf,
                            nb = 200,
                            smo = 0.05,
                            dfunc = depth.mode)$outliers
}


FOM <- function(data) {
    y = t(data)
    
    yt = array(y, dim = c(dim(y), 1))
    Result <- fOutl(yt,  type = "fAO", diagnostic = TRUE)
    
    
    # The user may opt to draw a cut off line separating the outliers.
    # which will be plotted in red
    ffom = fom(Result, cutoff = TRUE)
    out = which(ffom$data[, 4] == 'red')
    return(out)
}

FB <- function(data) {
    y = t(data)
    #Functional Boxplots (Sun and Genton 2011)
    out = fda::fbplot(y, plot = F)$outpoint
    return(out)
}


HDR <- function(data, t, out_rate) {
    # High density regions (Hydmann and Shang 2010)
    y = t(data)
    Datos = fds(t, y)
    out <-
        fhdr(
            Datos,
            alpha = c(out_rate, 0.5),
            projmethod = "PCAproj",
            xlab = '',
            ylab = '',
            plotlegend = F
        )
    return(out)
}

AAKNN <- function(data,
                  e = 2,
                  narch = 6,
                  vi = 5,
                  vf = 10, ai = NULL) {
    set.seed(1234)
    X = data
    nada = narch
    if (is.null(ai)) {
        ai = LstepArchetypesMod(X, 1:nada, nrep = 10, verbose = FALSE)
    }
    #Draw the RSS graph in order to find the elbow (number of initial archetypes)
    #screeplot(ai)
    
    h = 1
    pii = matrix(0, nrow = (vf - vi + 1), ncol = dim(X)[1])
    for (i in vi:vf) {
        pei = rep(0, dim(X)[1])
        for (j in e:nada) {
            a = bestModel(ai[[j]])
            pei = pei + KNN_SUM(a$alphas, k = i)
            
        }
    
        b = boxplot(pei, plot=FALSE)
        pii[h, ] = as.integer(is.element(pei, b$out))
        h = h + 1
    }
    
    oaaknn = as.integer(apply(pii, 2, sum) > (((vf - vi + 1)) * 0.5))
    out = which(oaaknn == 1)
    
    return(out)
}


FADALARA <- function(data, k, prob=0.75) {
    
    rownames(data) = 1:(dim(data)[1])
    X = t(data)
    num_pts = dim(X)[1]
    range = 1.5
    
    # ADA robust:
    # Cleaning:
    outl_ada_clean <- adamethods::do_clean(X, num_pts, range, 80)
    outl_ada_clean1 <- adamethods::do_clean(X, num_pts, 3)
    outl_ada_clean= union( outl_ada_clean, outl_ada_clean1)
    
    if (length(outl_ada_clean) == 0) { # For cases where cleaning didn't detect any outliers.
        X1 <- data
    } else{
        X1 <- data[-outl_ada_clean,]
    }
    ada <- adamethods::do_ada_robust(subset = X1, numArchoid = k, numRep = 20, huge = 200,
                                     prob = prob, compare = FALSE, method = "adjbox")
    # parse_number(rownames(X1)[ada$outliers]) is needed to identify the right outliers in
    # the data frame where some points where removed from the cleaning step.
    outl_ada <- c(outl_ada_clean, parse_number(rownames(X1)[ada$outliers]))
}