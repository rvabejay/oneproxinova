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



#Data processing


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


#Function that computes for the test statistic
F.test <- function(MSTr, MSE){
  fstat <- MSTr/MSE
  return(fstat)
}


#Function that computes for the critical value
#Using Paulson's Approximate transformation
F.tab <- function(df1, df2, alpha){
  #To obtain the z-value
  z = qnorm(1-alpha)

  a = ((1-(2/(9*df2)))^2)-((2*(z^2))/(9*(df2)))
  b = -(((1-(2/(9*df1))))*(1-(2/(9*(df2)))))
  c = ((1-(2/(9*df1)))^2)-((2*(z^2))/(9*(df1)))

  #By manipulating Paulson's approximation
  fcrit = ((-b + (sqrt(b^2-(a*c))))/a)^3
  return(fcrit)
}


#Computing for p-value
#Using Closed Form Approximation by Smillie and Anstey
P.F <- function(fstat, df1, df2){
  #defining constants
  a1 <- 0.278393
  a2 <- 0.230389
  a3 <- 0.000972
  a4 <- 0.078108

  u = ((1-(2/(9*df2)))*(fstat^(1/3))-(1-(2/(9*df1))))/sqrt(2/(9*df2)*(fstat^(2/3)) + 2/(9*df1))
  y = u/sqrt(2)

  f = 1-(0.5*(1+(a1*y)+(a2*y^2)+(a3*y^3)+(a4*y^4))^-4)
  return(1-f)
}


#Summary function for the Test on More than Two Populations
#Data should be a dataframe
#Factor and obs can be an integer (column index) or a character (column name)
#Function does not work for factor levels with unequal observations
SPM <- function(data, factor, obs, alpha = 0.05){
  #Check if function inputs are valid
  if(input.prelim(data, factor, obs, alpha)){

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

    #Degrees of Freedom
    #Treatment DF
    TrDF <- ncol(factor.data)-1

    #Exponential Error DF
    EDF <- nrow(main) - ncol(factor.data)

    #Total DF
    TDF <- nrow(main)-1

    #Vector for SS and MS
    SS_MS <- c()

    SS_MS <- ANOVA.SS_MS(main, factor.data, obs)

    TrSS <- SS_MS[1]
    ESS <- SS_MS[2]
    TSS <- SS_MS[3]
    MSTr <- SS_MS[4]
    MSE <- SS_MS[5]


    #Test Statistic
    test.stat <- F.test(MSTr, MSE)

    #Critical Value
    crit.val <- F.tab(TrDF, EDF, alpha)

    #P-value
    p.val <- P.F(test.stat, TrDF, EDF)

    #Vectors for Data Presentation
    source.variation <- c(factor, "Residuals")
    dfs <- c(TrDF, EDF)
    SS <- c(round(TrSS, 4), round(ESS, 4))
    MS <- c(round(MSTr, 4), round(MSE, 4))
    Fval <- c(round(test.stat, 4), NA)
    Pval <- c(signif(p.val, 4), NA)

    if(p.val < alpha){
      code <- c("*", NA)
    }

    #Bind vectors to a dataframe
    if(p.val < alpha){
      SUMMARY.ANOVA <- cbind(source.variation, dfs, SS, MS, Fval, Pval, code)
    }else{
      SUMMARY.ANOVA <- cbind(source.variation, dfs, SS, MS, Fval, Pval)
    }

    SUMMARY.ANOVA <- as.data.frame(SUMMARY.ANOVA)

    #Change column names
    if(p.val < alpha){
      colnames(SUMMARY.ANOVA) <- c("Source", "DF", "Sum Squares", "Mean Sum Squares", "F value", "P-value", " ")
    }else{
      colnames(SUMMARY.ANOVA) <- c("Souce", "DF", "Sum Squares", "Mean Sum Squares", "F value", "P-value")
    }

    #Print results
    cat("\n\t Results of One-Way ANOVA at", round((100*alpha), 2), "% level of signficance\n\n")
    cat("Factor:", factor,"\tResponse:",obs, "\n\n")
    print(SUMMARY.ANOVA, na.print = "", row.names = FALSE)

    if(p.val < alpha){
      cat("\nTest statistic", round(test.stat, 4), "is greater than the critical value", round(crit.val, 4), "at", TrDF, "and", EDF, "DF at alpha = ", alpha, ".\n")
      cat("P-value", signif(p.val, 4), "is less than alpha = ", alpha, ".\n")
      cat("\nAlternative Hypothesis:  At least one factor level has a different mean response. \n")
      cat("Post hoc analysis is recommended.\n\n")
    }else{
      cat("\nTest statistic", round(test.stat, 4), "is less than the critical value", round(crit.val, 4), "at", TrDF, "and", EDF, "DF at alpha = ", alpha, ".\n")
      cat("P-value", signif(p.val, 4), "is greater than alpha = ", alpha, ".\n")
      cat("\nNull Hypothesis: The means among the factor levels are all equal.\n\n")
    }

  }
}
