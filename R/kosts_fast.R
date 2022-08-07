
#' Calculate Kost's approximation of covariance
#'
#' @param cor_in correlation between two n x n data vectors.
#'
#' @return Kost's approximation of the covariance between the -log cumulative distributions. This is calculated with a cubic polynomial fit.
#' @export
#'
#' @examples
kostPolyFitFast <- function(cor_in) {
    a1 <- 3.263
    a2 <- 0.710
    a3 <- 0.027 # Kost's cubic coefficients
    (a1*cor_in + a2*cor_in^2 + a3*cor_in^3)
}

#' Calculate Kost Covariance
#'
#' @param data_matrix An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numeric matrix.
#' Note: Method does not deal with missing values.
#' @return An m x m matrix of pairwise covariances between the data vectors calculated using Kost's polynomial fit and the pearson correlation function.
#' @export
#'
#' @examples
calculateKostCovarianceFast <- function(data_matrix) {
    #cat("faster implementation\n");
    m = nrow(data_matrix)
    covar_matrix = mat.or.vec(m, m)

    covar_matrix = stats::cor(t(data_matrix));
    covar_matrix = kostPolyFitFast(covar_matrix);
    return(covar_matrix);
}

#calculateKostCovarianceOld <- function(data_matrix) {
#    m = nrow(data_matrix)
#    covar_matrix = mat.or.vec(m, m)
#    for (i in 1:m) {
#        for (j in i:m) {
#            res0 <- cor.test(data_matrix[i,], data_matrix[j,])
#            cor <- res0$estimate
#            p_val <- res0$p.value
#            covar = kostPolyFit(cor)
#            covar_matrix[i, j] = covar
#            covar_matrix[j, i] = covar
#          }
#        }
#    return(covar_matrix)
#}


#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numeric matrix
#       A numeric vector of m P-values to combine.
#Output: A combined P-value using Kost's Method.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
kostsMethodFast <- function(data_matrix, p_values, extra_info = FALSE) {
    covar_matrix <- calculateKostCovarianceFast(data_matrix)
    combinePValues(covar_matrix, p_values, extra_info = extra_info)
}


#Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
        df_brown = df_fisher
        c = 1.0
    }
    x = 2.0*sum( -log(p_values) )

    p_brown = stats::pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
    p_fisher = stats::pchisq(df=df_fisher, q=x, lower.tail=FALSE)

    if (extra_info) {
        return(list(P_test=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF=df_brown))
    }
    else {
        return(p_brown)
    }
}

