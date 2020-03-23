




onestep <- function (x,t, params) {  #function to calculate one step of stochastic SIR
  
  with(
    as.list(c(x, parms)), #lets us access variables and parameters stored in y and pars by name
    {
   
    #compute various additional populations  
    I.detected <-I1+I2+I3+I4
    I.undetected <- Iu1+Iu2+Iu3+Iu4
    I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infectious
    N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
    gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
    sigmai <- 6*sigma  # multiplier 6 for pseudo stages
    etat <- eta(t,w)     # case notification rate
    betat <- beta(t,w)   # time dependent transmissibility, presymptomatic=1 causes this transmissibility to apply to late stage latent cases as well

    rates <- as.numeric(c(betat*I.detected/N+betat*c*I.undetected/N+presymptomatic*betat*c*E6/N,                           # movements out of S
                          sigmai, # E1 to E2 flow
                          sigmai, # E2 to E3 flow
                          sigmai, 
                          sigmai, 
                          sigmai, 
                          sigmai, # E6 to either I1 or Iu1 based on detection rate  
                          gammai, 
                          gammai, 
                          gammai, 
                          gammai,                   # movements out of I (detected)
                    b, b, b, b,                                       # movements out of I (undetected)
                    etat))

    # transition probabilities
    p <- matrix(0, nrow=length(rates),ncol=length(x))   # matrix to hold transitions probs
                
         p[1,]  <- c(exp(-rates[1]*dt),1-exp(-rates[1]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # S-> E
         p[2,]  <- c(0, exp(-rates[2]*dt), 1-exp(-rates[2]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[3,]  <- c(0, 0, exp(-rates[3]*dt), 1-exp(-rates[3]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[4,]  <- c(0, 0, 0, exp(-rates[4]*dt), 1-exp(-rates[4]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[5,]  <- c(0, 0, 0, 0, exp(-rates[5]*dt), 1-exp(-rates[5]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[6,]  <- c(0, 0, 0, 0, 0, exp(-rates[6]*dt), 1-exp(-rates[6]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[7,]  <- c(0, 0, 0, 0, 0, 0, exp(-rates[7]*dt), (1-exp(-rates[7]*dt))*q(w), 0, 0, 0, (1-exp(-rates[7]*dt))*(1-q(w)), 0, 0, 0, 0, 0, 0)    # Transitions out of E
         p[8,]  <- c(0, 0, 0, 0, 0, 0, 0, exp(-rates[8]*dt), 1-exp(-rates[8]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I (detected), 1st dummy to 2nd dummy
         p[9,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[9]*dt), 1-exp(-rates[9]*dt), 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I, 2nd dummy to 3rd dummy
         p[10,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[10]*dt), 1-exp(-rates[10]*dt), 0, 0, 0, 0, 0, 0, 0)  # Transitions out of I, 3rd dummy to 4th dummy
         p[11,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[11]*dt), 0, 0, 0, 0, 1-exp(-rates[11]*dt), 0, 0)  # Transitions out of I -> H, 4th dummy to 
         p[12,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[12]*dt), 1-exp(-rates[12]*dt), 0, 0, 0, 0, 0)    # Transitions out of I (undetected)
         p[13,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[13]*dt), 1-exp(-rates[13]*dt), 0, 0, 0, 0)    # Transitions out of I
         p[14,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[14]*dt), 1-exp(-rates[14]*dt), 0, 0, 0)  # Transitions out of I
         p[15,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[15]*dt), 0, 1-exp(-rates[15]*dt), 0)  # Transitions out of I -> R_u
         p[16,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[16]*dt), 0, 1-exp(-rates[16]*dt))  # Transitions R_d -> to C (notification)

         # update states
         #matrix to hold future states
         xnew <- matrix(0, nrow=length(rates),ncol=length(x))                                #
         
         for(i in 1:length(rates)){
           states1[i,] <- t(rmultinom(1, x[i], p[i,])) 
         }
         
         states1 <- colSums(states1)
         states1[17] <- states1[17]+Ru  #add formerly Recovered undetected cases
         states1[18] <- states1[18]+C  #add formerly notified cases
         
         return(x <- c(dt, states1, tail(x,1)+dt))
       
         } #close as-list command
       )
}

model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,length(x)))         #set up array to store results
  colnames(output) <- c("time","S",
                        "E1", "E2", "E3", "E4", "E5", "E6",
                        "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4",
                        "H", "Ru", "C", "cum.time") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- as.numeric(onestep(x,params))
  }
  output                                    #return output
}

# moving from symptomatic (detected) to hospitalized
# takes 2 values, one before some intervention, one after
gamma <- function(z = 12, b=0.143, a0=1/1.5, t){
  # piecewise function
  # default parameters z = 12, b=1/7, a0=1/1.5
  #    z: time at start of intervention (notionally March 12)
  #    b: intercept (positive)
  #    a0: post intervention isolation ratae
  #    t: time in the model
  
  gamma <- ifelse(t<=z, gamma <- b, gamma <- a0)
  return(gamma)
}

eta <- function(t, w=12) ifelse(t<=w,1/3,1/3) #fraction of cases that are reported/notified. Different from those detected. Maybe not relevant for GA.

q <- function(t, w=12, q0=1, q1=1) ifelse(t<=w,q0,q1) #case detection rate, affects movement from E to either I-detected or I-undetected

beta <- function(t, w=12, beta0=0.6584, beta.factor=2) ifelse(t<=w,beta0,beta0/beta.factor)

evaluate.model <- function(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, c=1, presymptomatic=1, dt=0.05),
                                       init = list(S=10600000, E1=0, E2=0, E3=0, E4=0, E5=6, E6=0,
                                                   I1 = 1, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=2, nstep=NULL, start=as.Date("2020-03-01"),today=Sys.Date()){
 
  if(is.null(nstep)) nstep <- (as.numeric(today-start)+1+28)/params$dt #run simulation from start to current time plus four weeks

  xstart <- c(time=0, unlist(init), cum.time = 0) #initial conditions

  data <- vector(mode='list',length=nsims) #initialize list to store the output

  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    
  }

  return(data)
}


scenarios <- read.csv('Georgia scenarios - Sheet1.csv')

scenarios.output <- list()
for(i in 1:dim(scenarios)[1]){
  parms <- as.list(scenarios[i,2:10])
  init <- as.list(scenarios[i,11:28])
  scenarios.output[[i]] <- evaluate.model(params=parms, init = init, nsims=25, nstep=NULL, start=as.Date("2020-03-01"))
}

get.range <- function(simulations){
  # function to get the min and max total number of cases
  min <- min(unlist(lapply(simulations, function(x) tail(x$C,1))))
  max <- max(unlist(lapply(simulations, function(x) tail(x$C,1))))
  range <- c(min, max)
}

out <- lapply(scenarios.output, get.range)


start=as.Date("2020-03-01")
set.seed(1292020)                #set seed
out64.natural <- evaluate.model(params=list(beta0=0.6584, sigma=1/6.4, z=1200, b=0.143, a0=1/1.5, w=100, presymptomatic=1, c=1, dt=0.05),
                                       init = list(S=10600000, E1=6, E2=6, E3=6, E4=6, E5=6, E6=6,
                                                   I1 = 7, I2= 7, I3=7, I4=7, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=15, nstep=NULL, start=as.Date("2020-03-01"))

out64.baseline <- evaluate.model(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=100, c=1, presymptomatic=1, dt=0.05),
                                       init = list(S=10600000, E1=6, E2=6, E3=6, E4=6, E5=6, E6=6,
                                                   I1 = 7, I2= 7, I3=7, I4=7, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=15, nstep=NULL, start=as.Date("2020-03-01"))
plot.model(out64.baseline, log='y', title='Baseline: Rapid case identification but no social distancing')

out8 <- evaluate.model(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, c=1, presymptomatic=1, dt=0.05),
                                       init = list(S=10600000, E1=1, E2=1, E3=1, E4=1, E5=1, E6=1,
                                                   I1 = 1, I2= 1, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=15, nstep=NULL, start=as.Date("2020-03-01"))

out64 <- evaluate.model(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, c=1, presymptomatic=1, dt=0.05),
                                       init = list(S=10600000, E1=6, E2=6, E3=6, E4=6, E5=6, E6=6,
                                                   I1 = 7, I2= 7, I3=7, I4=7, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=15, nstep=NULL, start=as.Date("2020-03-01"))

out128 <- evaluate.model(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, c=1, presymptomatic=1, dt=0.05),
                                       init = list(S=10600000, E1=13, E2=13, E3=13, E4=13, E5=13, E6=13,
                                                   I1 = 13, I2= 13, I3=12, I4=12, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=15, nstep=NULL, start=as.Date("2020-03-01"))
plot.model(out128, log='y', title='Outbreak started with 128 cases on 1 March and interventions on 12 March')
plot.cases(out128, georgia)
