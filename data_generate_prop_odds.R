## Note that the dataset assumes the reverse order of the ordinal endpoint (i.e., 6 is death 
## instead of 1, 1 is back to normal activities instead of 6
## This changes the subtractive constants to additive constants (e.g., add 0.411 instead of subtracting
## 0.411

## Read in data
file.dir <- "C:/Users/pet00180/Downloads/"
flu = read.csv(paste(file.dir, "OrdScore.csv", sep =""))
flu2 = subset(flu, flutype <= 3 & symdur <= 7 & news >= 2 & !is.na(score1) & !is.na(score2)
              & !is.na(score3) & !is.na(score4) & !is.na(score5) & !is.na(score6) 
              & !is.na(score7))

## One patient was discharged from the hospital on the same day as randomization. Altering this patient's
## status to be in hospital, not in ICU, not on O2
flu2$score0_alt = factor(ifelse(flu2$score0 == 2, 3, flu2$score0))

## FLU-IVIG distributions of the placebo group across the follow-up period
day0p = prop.table(table(flu2$score0_alt))
day1p = prop.table(table(flu2$score1))
day2p = prop.table(table(flu2$score2))
day3p = prop.table(table(flu2$score3))
day4p = prop.table(table(flu2$score4))
day5p = prop.table(table(flu2$score5))
day6p = prop.table(table(flu2$score6))
day7p = prop.table(table(flu2$score7))

## FLU-IVIG placebo group transition matrices across the follow-up period
day0_to_1 = prop.table(table(flu2$score0_alt,flu2$score1),1)
day1_to_2 = prop.table(table(flu2$score1,flu2$score2),1)
day2_to_3 = prop.table(table(flu2$score2,flu2$score3),1)
day3_to_4 = prop.table(table(flu2$score3,flu2$score4),1)
day4_to_5 = prop.table(table(flu2$score4,flu2$score5),1)
day5_to_6 = prop.table(table(flu2$score5,flu2$score6),1)
day6_to_7 = prop.table(table(flu2$score6,flu2$score7),1)

## Derives distribution of time to event endpoint
time_event <- function(x) {
  time_to_event <- which(x < 3)
  if(length(time_to_event) > 0) {
    event <- 1
    time <- time_to_event[1] - 1
  } else {
    event <- 0
    time <- 7
  }
  return(c(time, event))
}

## Calculates empirical power for a given column of squared test statistic values
powfun = function(col){
  return(mean(col > qchisq(0.95,df=1)))
}

## The expit function
expit = function(x){exp(x)/(1+exp(x))}

## Converts cumulative log odds to probabilities
invlog = function(x){
  if(length(x) == 1){
    prob = expit(x)
    return(c(prob,1-prob))
  } else{
    len = length(x) + 1
    p = rep(NA,len)
    p[1] = expit(x[1])
    for(i in 2:(len-1)){
      p[i] = expit(x[i]) - expit(x[i - 1])
    }
    p[len] = 1 - sum(p[1:(len-1)])
    return(p)
  }
}

## Add constants to cumulative log odds
invcumlog = function(p,c=0){
  nzs = which(p>0)
  start = min(nzs)
  finish = max(nzs)
  n_log = log(cumsum(p[start:finish])/(1-cumsum(p[start:finish]))) + c
  f_log = n_log[!is.infinite(n_log)]
  probs = invlog(f_log)
  f_probs = c(rep(0,start-1),probs,rep(0,6 - finish))
  return(f_probs)
}

## Given a transition matrix, the sampling distribution on the current day, and a column number 
## to indicate the current day, finds the sampling distribution distribution of the next day
next_day = function(mat,t_day,col){
	next_vec = rep(0,160)
      lvs = as.numeric(levels(as.factor(mat[,col])))
	for(i in lvs){
		inds = which(mat[,col] == i)
		t_row = t_day[as.character(i),]
		t_row2 = t_row[t_row>0]
		if(length(t_row2) > 1){
			next_day = sample(x=as.numeric(names(t_row2)),size=length(inds),replace = TRUE,prob=t_row2)
		}else{
			next_day = rep(as.numeric(names(t_row2)),length(inds))
		}
		next_vec[inds] = next_day
	}
	return(next_vec)
}


## The add_ functions add constants to the cumulative log odds of the rows of the given
## transition matrix under certain conditions

add_c1 = function(mat,constant){
	new_mat = apply(mat,1,invcumlog,c = constant)
	z = t(new_mat)
	colnames(z) = c(1:6)
	rownames(z) = c(3:5)
	return(z)
}

add_c1111 = function(mat,constant){
	new_mat = t(apply(mat[2:3,],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1,],new_mat)
	rownames(merge_mat) = c(3:5)
	return(merge_mat)
}

add_c2 = function(mat,constant){
	ones = ceiling(which(mat==1)/ncol(mat))
	new_mat = t(apply(mat[(ones[1]+1):(ones[2]-1),],1,invcumlog, c = constant))
	ones_mat = mat[ones,]
	merge_mat = rbind(ones_mat[1,],new_mat,ones_mat[2:nrow(ones_mat),])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c22 = function(mat,constant){
	ones = ceiling(which(mat==1)/ncol(mat))
	new_mat = t(apply(mat[(ones[1]+2):(ones[2]-1),],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1:2,],new_mat,mat[ones[-1],])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c222 = function(mat,constant){
	ones = ceiling(which(mat==1)/ncol(mat))
	half_mat = invcumlog(mat[2,],c = constant * 3/4)
	new_mat = t(apply(mat[(ones[1]+2):(ones[2]-1),],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1,],half_mat,new_mat,mat[ones[-1],])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c2222 = function(mat,constant){
	ones = ceiling(which(mat==1)/ncol(mat))
	new_mat = t(apply(mat[(ones[1]+3):(ones[2]-1),],1,invcumlog, c = constant))
	ones_mat = mat[ones,]
	merge_mat = rbind(ones_mat[1,],mat[2:3,],new_mat,ones_mat[2:nrow(ones_mat),])
      rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c22222 = function(mat,constant){
	ones = ceiling(which(mat==1)/ncol(mat))
	new_mat = invcumlog(mat[(ones[1]+3):(ones[2]-1),],1)
	ones_mat = mat[ones,]
	merge_mat = rbind(ones_mat[1,],mat[2:3,],new_mat,ones_mat[2:nrow(ones_mat),])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c3 = function(mat,constant){
	new_mat = t(apply(mat[1:5,],1,invcumlog, c = constant))
	merge_mat = rbind(new_mat,mat[6,])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c33 = function(mat,constant){
	hosp_mat = t(apply(mat[3:5,],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1:2,],hosp_mat,mat[6,])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c333 = function(mat,constant){
	half_mat = invcumlog(mat[2,],c = constant * 3/4)
	hosp_mat = t(apply(mat[3:5,],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1,],half_mat,hosp_mat,mat[6,])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}

add_c3333 = function(mat,constant){
	hosp_mat = t(apply(mat[4:5,],1,invcumlog, c = constant))
	merge_mat = rbind(mat[1:3,],hosp_mat,mat[6,])
	rownames(merge_mat) = c(1:6)
	return(merge_mat)
}


## Add (subtract) a constant (e.g., 0.5) to add fewer (more) severe cases of influenza
## Add 0 to use the FLU-IVIG placebo group transition matrices
p = 0
day0_to_1 = add_c1(day0_to_1,p)
day1_to_2 = add_c2(day1_to_2,p)
day2_to_3 = add_c2(day2_to_3,p)
day3_to_4 = add_c2(day3_to_4,p)
day4_to_5 = add_c3(day4_to_5,p)
day5_to_6 = add_c2(day5_to_6,p)
day6_to_7 = add_c2(day6_to_7,p)

## T1: Treatment effect remains constant for all categories across all seven days of follow-up.
g = 0.411
e0_to_1 = add_c1(day0_to_1,g)
e1_to_2 = add_c2(day1_to_2,g)
e2_to_3 = add_c2(day2_to_3,g)
e3_to_4 = add_c2(day3_to_4,g)
e4_to_5 = add_c3(day4_to_5,g)
e5_to_6 = add_c2(day5_to_6,g)
e6_to_7 = add_c2(day6_to_7,g)

## T2: Treatment effect only benefits patients for the first three days after randomization 
## for all categories.
#g = 0.891
#e0_to_1 = add_c1(day0_to_1,g)
#e1_to_2 = add_c2(day1_to_2,g)
#e2_to_3 = add_c2(day2_to_3,g)
#e3_to_4 = add_c2(day3_to_4,0)
#e4_to_5 = add_c3(day4_to_5,0)
#e5_to_6 = add_c2(day5_to_6,0)
#e6_to_7 = add_c2(day6_to_7,0)

## T3: Treatment effect constantly decreases with each successive day with no additional 
## benefit on day 7 for all categories
#g = 0.787
#sl = -g/6
#e0_to_1 = add_c1(day0_to_1,g)
#e1_to_2 = add_c2(day1_to_2,g + sl)
#e2_to_3 = add_c2(day2_to_3,g + 2*sl)
#e3_to_4 = add_c2(day3_to_4,g + 3*sl)
#e4_to_5 = add_c3(day4_to_5,g + 4*sl)
#e5_to_6 = add_c2(day5_to_6,g + 5*sl)
#e6_to_7 = add_c2(day6_to_7,0)

## T4: Treatment effect remains constant across all seven days of follow-up but is 75% as 
## effective for discharged patients.
#g = 0.424
#e0_to_1 = add_c1(day0_to_1,g)
#e1_to_2 = add_c222(day1_to_2,g)
#e2_to_3 = add_c222(day2_to_3,g)
#e3_to_4 = add_c222(day3_to_4,g)
#e4_to_5 = add_c333(day4_to_5,g)
#e5_to_6 = add_c222(day5_to_6,g)
#e6_to_7 = add_c222(day6_to_7,g)

## T5: Treatment effect only benefits patients in the In ICU and non-ICU, on oxygen categories 
## across all seven days of follow-up.
#g = 0.834
#e0_to_1 = add_c1111(day0_to_1,g)
#e1_to_2 = add_c2222(day1_to_2,g)
#e2_to_3 = add_c2222(day2_to_3,g)
#e3_to_4 = add_c2222(day3_to_4,g)
#e4_to_5 = add_c3333(day4_to_5,g)
#e5_to_6 = add_c2222(day5_to_6,g)
#e6_to_7 = add_c22222(day6_to_7,g)

## T6: Treatment effect only benefits patients for the last three days of follow-up for all 
## categories.
#g = 0.716
#e0_to_1 = add_c1(day0_to_1,0)
#e1_to_2 = add_c2(day1_to_2,0)
#e2_to_3 = add_c2(day2_to_3,0)
#e3_to_4 = add_c2(day3_to_4,0)
#e4_to_5 = add_c3(day4_to_5,g)
#e5_to_6 = add_c2(day5_to_6,g)
#e6_to_7 = add_c2(day6_to_7,g)

## T7: Treatment effect constantly increases with each successive day, starting on day 2 for 
## all categories. No treatment effect on day 1.
#sl = 0.666/6
#e0_to_1 = add_c1(day0_to_1,0)
#e1_to_2 = add_c2(day1_to_2,sl)
#e2_to_3 = add_c2(day2_to_3,2*sl)
#e3_to_4 = add_c2(day3_to_4,3*sl)
#e4_to_5 = add_c3(day4_to_5,4*sl)
#e5_to_6 = add_c2(day5_to_6,5*sl)
#e6_to_7 = add_c2(day6_to_7,6*sl)

## More and less severe day 0 placebo group distributions
day0p_severe = day0p
day0p_better = day0p
day0p_severe[1:3] = invcumlog(day0p, -0.5)[1:3]
day0p_better[1:3] = invcumlog(day0p, 0.5)[1:3]

## Data generating function
t_gen = function(){

	day0vec = sample(3:5, 160, day0p, replace = TRUE)
	s_mat = cbind(day0vec)
	s_mat2 = cbind(s_mat,next_day(s_mat,day0_to_1,1))
	s_mat3 = cbind(s_mat2,next_day(s_mat2,day1_to_2,2))
	s_mat4 = cbind(s_mat3,next_day(s_mat3,day2_to_3,3))
	s_mat5 = cbind(s_mat4,next_day(s_mat4,day3_to_4,4))
	s_mat6 = cbind(s_mat5,next_day(s_mat5,day4_to_5,5))
	s_mat7 = cbind(s_mat6,next_day(s_mat6,day5_to_6,6))
	s_mat8 = cbind(s_mat7,next_day(s_mat7,day6_to_7,7))
	s_mat_plac = cbind(s_mat8,rep(0,160))

	day0vec = sample(3:5, 160, day0p, replace = TRUE)
	e_mat = cbind(day0vec)
	e_mat2 = cbind(e_mat,next_day(e_mat,e0_to_1,1))
	e_mat3 = cbind(e_mat2,next_day(e_mat2,e1_to_2,2))
	e_mat4 = cbind(e_mat3,next_day(e_mat3,e2_to_3,3))
	e_mat5 = cbind(e_mat4,next_day(e_mat4,e3_to_4,4))
	e_mat6 = cbind(e_mat5,next_day(e_mat5,e4_to_5,5))
	e_mat7 = cbind(e_mat6,next_day(e_mat6,e5_to_6,6))
	e_mat8 = cbind(e_mat7,next_day(e_mat7,e6_to_7,7))
	e_mat_trt = cbind(e_mat8,rep(1,160))

	c_mat = rbind(e_mat_trt,s_mat_plac)
	colnames(c_mat) = c("Day_0", "Day_1", "Day_2", "Day_3", "Day_4", "Day_5", "Day_6", "Day_7", "trt")
	d_mat = as.data.frame(c_mat)

	return(d_mat)
}

## Calculates empirical power across results from parallel simulation
mp = function(results,l){
	full = NULL
	for(i in 1:l){
	full = rbind(full,results[[i]]$test_stats)
	}
	return(full)
}

## Calculates average coefficient values across results from parallel simulation
mt = function(results,l){
	full = NULL
	for(i in 1:l){
	full = rbind(full,results[[i]]$trt_stats)
	}
	return(full)
}

ms = function(results,l){
	full = NULL
	for(i in 1:l){
	full = rbind(full,results[[i]]$day7_bin)
	}
	return(full)
}

## Vector of starting values
guesses = c(-4.4107760,-2.6661592,-1.2253290,-0.5279292,1.0201407,-0.5709795)

## Given a vector of guesses for estimating the parameter values of the proportional odds model (gs), 
## a placebo group (plac_prob), and an experimental group (exp_prob), returns the values of the 
## partial derivatives of the expected score functions.
## Used to find the analytic odds ratio for the proportional odds model.
## Method devised by Peterson et al. (2017)
paramfnGen = function(gs,plac_prob,exp_prob){
  lvs = length(plac_prob)
  store = NULL
  dgpen = 0
  dgind2 = 0
  dgult = 0
  
  dg1 = 
    1/(exp(gs[1]) + 1) * plac_prob[1] + 
    (exp(gs[1])*(exp(gs[2]) + 1))/((exp(gs[1]) + 1)*(exp(gs[1]) - exp(gs[2]))) * plac_prob[2] + 
    1/(exp(gs[1] + gs[lvs]) + 1) * exp_prob[1] +
    (exp(gs[1])*(exp(gs[2] + gs[lvs]) + 1))/((exp(gs[1] + gs[lvs]) + 1)*(exp(gs[1]) - exp(gs[2]))) * exp_prob[2]
  
  if(lvs > 3){
    for(i in 1:(lvs - 3)){
      dgx = 
        (exp(gs[i+1])*(exp(gs[i]) + 1))/((exp(gs[i+1]) + 1)*(exp(gs[i+1]) - exp(gs[i]))) * plac_prob[i+1] + 
        (exp(gs[i+1])*(exp(gs[i+2]) + 1))/((exp(gs[i+1]) + 1)*(exp(gs[i+1]) - exp(gs[i+2]))) * plac_prob[i+2] +
        (exp(gs[i+1])*(exp(gs[i] + gs[lvs]) + 1))/((exp(gs[i+1] + gs[lvs]) + 1)*(exp(gs[i+1]) - exp(gs[i]))) * exp_prob[i+1] +
        (exp(gs[i+1])*(exp(gs[i+2] + gs[lvs]) + 1))/((exp(gs[i+1] + gs[lvs]) + 1)*(exp(gs[i+1]) - exp(gs[i+2]))) * exp_prob[i+2]
      
      store = c(store, dgx)
      }
    }
    
  for(i in 1:(lvs - 2)){
  dgind = 
    exp(gs[lvs-1])/(exp(gs[lvs-1]) + 1) * plac_prob[i] +
    exp(gs[lvs-1] + gs[lvs])/(exp(gs[lvs-1] + gs[lvs]) + 1) * exp_prob[i]
    dgpen = dgpen + dgind
    }
  dgpen = dgpen + 
    exp(gs[lvs-1])/(exp(gs[lvs-1]) - exp(gs[lvs-2])) * plac_prob[lvs-1] - exp(gs[lvs-1])/(exp(gs[lvs-1]) + 1) + 
    exp(gs[lvs-1])/(exp(gs[lvs-1]) - exp(gs[lvs-2])) * exp_prob[lvs-1] - exp(gs[lvs-1] + gs[lvs])/(exp(gs[lvs-1] + gs[lvs]) + 1)
  
  dgult =
    (1/(exp(gs[1] + gs[lvs]) + 1) - 1/(exp(gs[lvs-1] + gs[lvs]) + 1) + 1) * exp_prob[1]

  if(lvs > 3){
    for(i in 1:(lvs - 3)){
        dgind2 = dgind2 + (1/(exp(gs[i] + gs[lvs]) + 1) + 1/(exp(gs[i+1] + gs[lvs]) + 1) - 1/(exp(gs[lvs - 1] + gs[lvs]) + 1)) * exp_prob[i + 1]
      }
  }
  dgult = dgult + dgind2 + 1/(exp(gs[lvs - 2] + gs[lvs]) + 1) * exp_prob[lvs - 1] - exp(gs[lvs-1] + gs[lvs])/(exp(gs[lvs-1] + gs[lvs]) + 1)
  
  return(c(dg1, store, dgpen, dgult))
}


## Given a dataset, calculates the win ratio and corresponding test statistic
WR_fun1 = function(comb){
  n = 320
  ntrt = n/2
  nw = 0
  nl = 0
  ut = NULL
  for (j in 1:n){ # For each patient
    uw = 0
    ul = 0
    for(k in 1:(n-1)){
      if(comb[j] < comb[-j][k]){ # Compare each patient to all others
        uw = uw + 1
        if(j < (ntrt + 1) & k > (ntrt - 1)){ # Only if comparing experimental with placebo
          nw = nw + 1                                 
        }
      } else if (comb[j] > comb[-j][k]){
        ul = ul + 1
        if(j < (ntrt + 1) & k > (ntrt - 1)){
          nl = nl + 1
        }
      }
    }
    ut = rbind(ut, uw - ul)
  }
  variance = (ntrt^2 / (n * (n-1)) * t(ut)%*%ut)/nl^2
  win_ratio = nw/nl
  z = (win_ratio - 1)/sqrt(variance)
  return(list(z = z, win_ratio = win_ratio))
}

## Given a dataset, calculates the win ratio and corresponding test statistic with adjustment for
## baseline status
WR_fun5 = function(dset){
  freq = table(dset$Day_0)
  numstrat = length(freq)
  mincat = as.integer(names(freq)[1])
  maxcat = as.integer(names(freq)[length(freq)])
  nws = rep(0,numstrat)
  nls = rep(0,numstrat)
  vars = rep(0,numstrat)
  for(i in mincat:maxcat){ # For each strata
    trt = dset[dset$trt == 1 & dset$Day_0 == i,]
    plac = dset[dset$trt == 0 & dset$Day_0 == i,]
    comb = rbind(trt,plac)$Day_7
    ntrt = nrow(trt)
    nplac = nrow(plac)
    tot = nplac + ntrt
    nw = 0
    nl = 0
    ut = NULL
    if(tot > 1){
      for (j in 1:tot){ # For each patient
        uw = 0
        ul = 0
        for(k in 1:(tot-1)){
          if(comb[j] < comb[-j][k]){ # Compare each patient to all others
            uw = uw + 1
            if(j < (ntrt + 1) & k > (ntrt - 1)){ # Only if comparing experimental
              nw = nw + 1                                  # patient with placebo patient
            }
          } else if (comb[j] > comb[-j][k]){
            ul = ul + 1
            if(j < (ntrt + 1) & k > (ntrt - 1)){
              nl = nl + 1
            }
          } 
        }
        ut = rbind(ut, uw - ul)
      }
      nws[i] = nw # Store stratum wins
      nls[i] = nl
      vars[i] = ntrt * nplac / (tot * (tot-1)) * t(ut)%*%ut # Store stratum variance
    }
  }
  win_ratio = sum(nws)/sum(nls)
  vd = sum(vars)/sum(nls)^2
  z = (sum(nws)/sum(nls) - 1)/sqrt(vd)
  return(list(z = z, win_ratio = win_ratio, nws = sum(nws), nls = sum(nls)))
}
