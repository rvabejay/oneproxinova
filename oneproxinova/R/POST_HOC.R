#Data Pre-processing

#Function checks if inputs are valid
input.prelim <- function(data, factor, obs, alpha){

  #Check if data exists
  if(exists(deparse(substitute(data)))){

    #Check if data is a dataset
    if(!(is.data.frame(data))){
      cat("Convert the data as a dataframe.\n")
      return(FALSE)
    }else{

      #For tibbles (type of dataframe, done to avoid some errors)
      data <- as.data.frame(data)

      is.factor.valid <- (is.character(factor)) || (is.numeric(factor))
      is.obs.valid <- (is.character(obs)) || (is.numeric(obs))

      #Check if factor and obs are valid (characters or integers)
      if(!(is.factor.valid && is.obs.valid)){
        cat("Inputs for factor and/or obs are not valid. Please enter the column name or the column index.\n")
        return(FALSE)
      }else{

        #Checks if factor and obs are correct (correct indexes or correct column names)

        is.factor.correct <- NULL
        is.obs.correct <- NULL

        if(is.numeric(factor)){
          if((factor <= ncol(data)) && (factor %% 1 == 0)){
            is.factor.correct <- TRUE
          }else{
            is.factor.correct <- FALSE
          }
        }else{
          if(factor %in% colnames(data)){
            is.factor.correct <- TRUE
          }else{
            is.factor.correct <- FALSE
          }
        }


        if(is.numeric(obs)){
          if((obs <= ncol(data)) && (obs %% 1 == 0)){
            is.obs.correct <- TRUE
          }else{
            is.obs.correct <- FALSE
          }
        }else{
          if(obs %in% colnames(data)){
            is.obs.correct <- TRUE
          }else{
            is.obs.correct <- FALSE
          }
        }

        if(is.factor.correct && is.obs.correct){

          #Check if alpha input is valid
          if((alpha < 1) && (alpha > 0)){
            return(TRUE)
          }else{
            cat("Input for alpha is invalid. Please enter an appopriate significance level.")
            return(FALSE)
          }

        }else{
          cat("Input for factor/obs invalid. Please enter the correct column index/name/s.")
          return(FALSE)
        }
      }
    }
  }else{
    cat("The data does not exist.\n")
    return(FALSE)
  }
}


#Creates a subset of the data with only the observations of the variables of interest
base <- function(data, factor, obs){
  factor.vector <- data[, factor]
  obs.vector <- data[, obs]
  newdata <- cbind(factor.vector, obs.vector)
  newdata <- as.data.frame(newdata)
  colnames(newdata) <- c(factor, obs)
  return (newdata)
}


#Creates a dataframe with columns as the factor levels
subdata <- function(data, factor, obs){
  column <- as.vector(data[,factor])

  #Get the unique values of the column vector to get the factor levels
  factor.levels <- unique(column)

  #Create a data frame with column length as the number of factor levels
  grouped.data <- data.frame(matrix(ncol = length(factor.levels)))

  #Fill the dataframe
  for(i in 1:length(factor.levels)){
    level.data <- c()

    #Check if the value in the dataframe corresponds to a certain factor level
    for(j in 1:nrow(data)){
      if(data[j,factor] == factor.levels[i]){
        level.data <- append(level.data, data[j,obs])
      }
    }

    #Bind to the dataframe
    grouped.data <- cbind(grouped.data, level.data)

    i <- i + 1
  }

  #Remove empty columns
  grouped.data <- subset(grouped.data, select = -c(1:length(factor.levels)))

  #Set the factor levels as the column names of the data frame
  colnames(grouped.data) <- factor.levels

  return(grouped.data)
}

#Function that computes for the sum of squares and mean sum of squares
ANOVA.SS_MS <- function(main, factor.data, obs){
  #Correction Factor
  CF <- ((sum(main[,obs]))^2)/(nrow(main))

  #Total SS
  TSS <- (sum(main[,obs]^2)) - CF

  #Treatment SS (Factor SS)
  total.level.sum <- 0

  for(i in 1:ncol(factor.data)){
    level.sum <- sum(factor.data[,i])
    level.sum.squared <- level.sum^2
    level.weight <- level.sum.squared/nrow(factor.data)

    total.level.sum <- total.level.sum + level.weight
  }

  TrSS <- total.level.sum - CF

  #Error SS
  ESS <- TSS - TrSS

  #MSTr
  MSTr <- TrSS/(ncol(factor.data)-1)

  #MSE
  MSE <- ESS/(nrow(main)-ncol(factor.data))

  return(c(TrSS, ESS, TSS, MSTr, MSE))
}

#computes for the point estimate of the mean difference between the two factor levels
point.estimate <- function(data, x, y){

  #create separate vectors
  x.data <- data[,x]
  y.data <- data[,y]

  #compute for the individual means
  x.mean <- mean(x.data)
  y.mean <- mean(y.data)

  #compute for the difference
  p.estimate <- x.mean - y.mean
  return(p.estimate)
}

#computes for the standard error
standard.error <- function(MSE, reps){
  std.error <- sqrt(MSE*(2*(1/reps)))
  return(std.error)
}


#computes for the critical value
crit.val <- function(alpha, noOfLevels, EDF, MSE, reps){
  #compute for the q-value
  q.val <- qtukey(alpha, noOfLevels, EDF, lower.tail = FALSE)
  constant <- sqrt(MSE/reps)

  crit <- q.val*constant
  return(crit)
}



#computes for the t-approximate value using q quantile value
t.approx <- function(mean, se){
  tapprox <- abs(mean)/se
  return(tapprox)
}


#returns a matrix of the possible level combinations
combinations <- function(grouped.data){
  levels <- colnames(grouped.data)
  combs <- combn(levels, 2)
  return(combs)
}


#Summary function for post-hoc analysis
#Function only works for datasets with equal observations for every factor level
POST.HOC <- function(data, factor, obs, alpha = 0.05){
  #Check if function inputs are valid
  if(input.prelim(data, factor, obs, alpha)){

    #Convert numeric inputs to the column/variable name corresponding it
    if(is.numeric(factor)){
      factor <- colnames(data)[factor]
    }

    if(is.numeric(obs)){
      obs <- colnames(data)[obs]
    }

    #Create a dataframe containing the variables of interest
    main <- base(data, factor, obs)

    #Create dataframe with columns as the factors
    factor.data <- subdata(main, factor, obs)

    #Identify possible combinations of factors
    factor.combinations <- combinations(factor.data)

    #No. of Levels
    no.of.levels <- ncol(factor.data)

    #Error DF
    EDF <- nrow(main) - ncol(factor.data)

    #MSE
    MSE <- ANOVA.SS_MS(main, factor.data, obs)[5]

    #Number of repetitions per factor level (assumed to be equal)
    reps <- nrow(factor.data)

    #Empty dataframe for the output
    POSTHOC_SUMMARY = data.frame()

    #Vectors to store outputs at each iteration
    p.estimates <- numeric()
    s.errors <- numeric()
    t.values <- numeric()
    significant <- c()

    #Critical value
    critical <- crit.val(alpha, no.of.levels, EDF, MSE, reps)

    #Repeat the process for every possible combination
    for(i in 1:ncol(factor.combinations)){
      #Point Estimate
      x <- factor.combinations[1,i]
      y <- factor.combinations[2,i]
      point.est <- point.estimate(factor.data, x, y)

      #Standard Error
      se <- standard.error(MSE, reps)

      #T-approximate value
      t.val <- t.approx(point.est, se)

      #Append to Vectors
      p.estimates <- append(p.estimates, abs(round(point.est,4)))
      s.errors <- append(s.errors, signif(se,4))
      t.values <- append(t.values, signif(t.val,4))

      if(signif(t.val, 4) > critical){
        significant <- append(significant, "*")
      }else{
        significant <- append(significant, " ")
      }
    }

    #Bind vectors vertically
    POSTHOC_SUMMARY <- cbind(p.estimates, s.errors, t.values, significant)

    #Create vector for the factor comparison combinations
    comparisons <- c()
    for(i in 1:ncol(factor.combinations)){
      text1 <- as.character(factor.combinations[1, i])
      text2 <- as.character(factor.combinations[2, i])
      text3 <- paste(text1, "-", text2, "== 0")
      comparisons <- append(comparisons, text3)
    }

    POSTHOC_SUMMARY <- as.data.frame(POSTHOC_SUMMARY)

    #Set column names and row names
    colnames(POSTHOC_SUMMARY) <- c("Pt. Estimate", "Std. Error", "t-value", "")
    rownames(POSTHOC_SUMMARY) <- comparisons

    #Print results
    cat("\n\t Results of Pairwise Mean Comparisons at", round((100*alpha), 2), "% level of signficance\n\n")
    cat("Factor:", factor,"\tResponse:",obs, "\n\n")
    print(POSTHOC_SUMMARY)

    cat("\nNull hypothesis: The mean of the two factor levels are equal.\n")
    cat("Reject Ho if t-value is greater than the critical value =", critical)
    cat("\n\nNote: Significance of t-value is approximate.")
    cat("\nUse the R function TukeyHSD to determine the p-values at each comparison.\n")
  }
}



