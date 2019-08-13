library(ggplot2)
library(nlmixr)
library(RxODE)
library(dplyr)
library(tidyr)
str(theo_sd)

#derivative of concentration with respect to time (change in concentration over time) = dC/Dt = (-Cl/V))*C = -K * C
#Cl = clearance = dose/Area under the curve
#V = Distribution volume = Dose/Initial concentration
#Kel = elimination rate constant = ln(1st concentration)-ln(last concentration)/change in time
#ka = constant of absorption
#lag = lag time

#one compartmental iv bolus dosing model

#Cl = clearance
#V = Distribution volume 
#ka = constant of absorption

one_cmpt_iv_bolus_dosing_model <- function() { #initialize the model
  ini({ #initialize the parameters of the model
    tka <- .5
    tcl <- -3.2
    tv <- -1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    add.err <- 0.1
  })
  model({ #model it
    ka <- exp(tka +eta.ka)
    Cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.err)
  })
}





#fit function runs the nlmixr model 
fit <- nlmixr(nikhil_one_cmpt_iv_bolus_dosing_model, theo_sd, est="nlme") #one.cmt = above initialized function with parameters defined, data, estimation method


nlmixr::vpc(fit, n=600, show=list(obs_dv=T), log_y=T, xlab="Time (h)", ylab="Concentration") #generate visual predictive check
plot(fit) #make plots, including residuals and individual plots 


ggplot(fit,aes(TIME, DV))+geom_point()+facet_wrap(~ID)+geom_line(aes(TIME,IPRED),col="red",size=1.2)+geom_line(aes(TIME,PRED),col="blue",size=1.2);

#plots that predict points outside of observed data
plot(augPred(fit));

#simulate a new regimen
et <-eventTable() %>%
  add.dosing(4.7, nbr.doses=2, dosing.interval=12) %>%
  add.sampling(seq(0,24,length.out=51))

n <- 300
bid <-simulate(fit,events=et,nSub=n*n)

p <-c(0.05, 0.5, 0.95);
s <-bid %>% mutate(id=sim.id%%n) %>% group_by(id,time) %>%
  do(data.frame(p1=p, eff=quantile(.$sim, probs=p))) %>%
  group_by(pl,time) %>%
  do(data.frame(p2=p, eff=quantile(.$eff, probs=p))) %>%
  ungroup()  %>% mutate(p2=sprintf("p%02d",(p2*100))) %>%
  spread(p2,eff) %>% mutate(Percentile=factor(sprintf("%02d%%",round(p1*100))))

ggplot(s,aes(time,p50,col=Percentile,fill=Percentile) )+ geom_ribbon(aes(ymin=p05,ymax=p95),alpha=0.5) + geom_line(size=1.2) + theme_light(base_size=18) + xlab("Time (hr)") + ylab("Theophylline (\u03BCg/mL)")


#one compartmental oral dosing model

nlmixr::Oral_1CPT #dataset from nlmixr documentation that can be used to validate the accuracy of the model


one_cmpt_oral_dosing_model <- function() { 
  ini({
    tka <- 0.45 # Log Ka
    tcl <- 1 # Log Cl
    tv <- 3.45    # Log V
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.err <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.err)
  })
}

#run the model
fit <- nlmixr(one_cmpt_oral_dosing_model, Oral_1CPT, est="nlme") 


#test

oral_1_compartmental_dosing_model <- function(){
  ini({
    lCl <- 1.6      #log Cl (L/hr) #clearance
    lVc <- log(90)  #log Vc (L) # volume of distribution
    lKA <- 0.1      #log Ka (1/hr) #constant of absorption
    prop.err <- c(0, 0.2, 1)
    eta.Cl ~ 0.1 ## BSV Cl
    eta.Vc ~ 0.1 ## BSV Vc
    eta.KA ~ 0.1 ## BSV Ka
  })
  model({
    ## First parameters are defined in terms of the initial estimates
    ## parameter names.
    Cl <- exp(lCl + eta.Cl)
    Vc <- exp(lVc + eta.Vc)
    KA <- exp(lKA + eta.KA)
    ## After the differential equations are defined
    kel <- Cl / Vc; #constant of elimination 
    
    #the code below is a way of doing 
    d/dt(depot)    = -KA*depot;
    d/dt(centr)  =  KA*depot-kel*centr;
    ## And the concentration is then calculated
    cp = centr / Vc;
    ## Last, nlmixr is told that the plasma concentration follows
    ## a proportional error (estimated by the parameter prop.err)
    cp ~ prop(prop.err)
  })
}

fit.s <- nlmixr(oral_1_compartmental_dosing_model,Oral_1CPT,est="saem",control=saemControl);
nlmixr::vpc(fit.s, n=600, show=list(obs_dv=T), log_y=T, xlab="Time (h)", ylab="Concentration") #generate visual predictive check
plot(fit.s) #make plots, including residuals and individual plots 
