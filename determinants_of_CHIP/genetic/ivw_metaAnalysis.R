
### Perform inverse-variance weighted fixed effect meta-analysis in R
# P values are from two-sided Z-tests
################### Inverse-variance weighted meta-analysis ####
ivw_meta <- function(b, bSE) {
  # Calculate weights
  weights <- 1 / (bSE^2)
  
  # Calculate overall standard error
  overall_se <- sqrt(1 / sum(weights))
  
  # Calculate overall effect size
  overall_beta <- sum(b * weights) / sum(weights)
  
  # Calculate Z-score
  Z <- overall_beta / overall_se
  
  # Calculate p-value
    # P values are from two-sided Z-tests
  pval <- 2 * (1 - pnorm(abs(Z)))
  # pval <- 2 * pnorm(abs(Z), lower.tail = FALSE)
  return(list(beta = overall_beta, 
              se = overall_se, 
              Z = Z, 
              pval = pval))
  
}


#ivw_meta2 <- function(b1, b2, bSE1, bSE2){
#	
	# inverse-variance weight
#	w1 <- 1/(bSE1^2)
	
#	w2 <- 1/(bSE2^2)
      	
	# SE meta
 #       se <- sqrt(1/sum(w1, w2) )
	
	# beta meta
#	beta <- sum(b1*w1, b2*w2)/sum(w1, w2)
	
	# Z meta    
#	Z <- beta/se
	        
#	pval <- 2 * pnorm(abs(Z), lower.tail = FALSE)
		    
	# 
#	return(list(beta, se, Z, pval))
# }


# ivw_meta <- function(b, bSE) {
# 	# Calculate weights
# 	weights <- 1 / (bSE^2)
#       
#         # Calculate overall standard error
# 	overall_se <- sqrt(1 / sum(weights))
# 	    
# 	# Calculate overall effect size
# 	overall_beta <- sum(b * weights) / sum(weights)
# 	
# 	# Calculate heterogeneity statistics
# 	        # Cochran's Q statistic is a statistical test used to assess the heterogeneity of effect sizes across studies in a meta-analysis. It quantifies the extent of inconsistency or variability among the effect sizes.
# 	Q <- sum(weights * (b - overall_beta)^2)
# 	df <- length(b) - 1
# 	I2 <- ((Q - df) / Q) * 100
# 	# heterogeneity_pval <- 1 - pchisq(Q, df)
# 	cochran_pval <- 1 - pchisq(Q, df)
# 
# 	# Calculate Z-score
# 	Z <- overall_beta / overall_se
# 		        
# 	# Calculate p-value
# 	pval <- 2 * (1 - pnorm(abs(Z)))
# 
# 	### 
# 	return(list(beta = overall_beta, se = overall_se, Z = Z, pval = pval, 
# 	q_statistic = Q, q_p_value = cochran_pval, i2 = I2))
# 
# }



######################################
