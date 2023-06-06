options(warn=2) # warnings are treated as errors

# source: https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector
number2binary <- function(number, noBits) 
{
	binary_vector = rev(as.numeric(intToBits(number)))
	if(missing(noBits)) 
	{
		return(binary_vector)
	} else 
	{
		binary_vector[-(1:(length(binary_vector) - noBits))]
	}
}

# call parameters
source("pars.R")

tick_max=399
Cdh1_array <- c()
SBF_array <- c()
Cln2_array <- c()
Clb5_array <- c()
Clb2G_array <- c()
Clb2M_array <- c()
Cdc20_array <- c()
Size_array <- c()

##################################
## PART 1 = Follow mother cells ##
##################################

initials.D <- c() # to hold initial conditions for new daughter cells
	
j=64 # j=64 is for state '1000000' (G1); j=65 is for state '1000001' (cell at birth)	
	
for(rep in 1:20) # repeats (20 cells)
{
	number2binary(j, 7) -> initial
	
	# Initial condition
	Cdh1=initial[1] 
	SBF=initial[2]
	Cln2=initial[3]
	Clb5=initial[4]
	Clb2G=initial[5]
	Clb2M=initial[6]
	Cdc20=initial[7]
	Size=0.65
	tr=0;
	
	tick <- 0
	tick_timecourse <- c()
	tr_timecourse <- c()
	Cdh1_timecourse <- c()
	SBF_timecourse <- c()
	Cln2_timecourse <- c()
	Clb5_timecourse <- c()
	Clb2G_timecourse <- c()
	Clb2M_timecourse <- c()
	Cdc20_timecourse <- c()
	Size_timecourse <- c()
	
	state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
	colnames(state)[10] <- 'Phase'
			
	S0=rlnorm(1,log(S0.mean), S0.CV)
	
	while(1) # will break out when tick > tick_max
	{
	    	# Calculate new values by Boolean functions
	    	Cdh1new = !( Cln2 || Clb5 || Clb2G || Clb2M ) || ( Cdc20 && !( Cln2 || Clb5 || Clb2M ) )
		SBFnew = ( SBF && !( Clb2G || Clb2M ) ) || 
			( 
				( Cdh1 && !SBF && !Cln2 && !Clb5 && !Clb2G && !Clb2M && !Cdc20 ) &&
			  	( Size > S0 ) && ( runif(1) < (Size-S0)^2 )
			)
		Cln2new = SBF
		Clb5new = ( Clb5 || SBF ) && !( Cdh1 || Cdc20 )
		Clb2Gnew = Clb2M || ( ( Clb5 || Clb2G ) && !Cdh1 )
		Clb2Mnew = ( Clb2G || Clb2M ) && !Cdh1 && ( !Cdc20 || Clb5 ) && !Cln2
		Cdc20new = Clb2M || ( Clb2G && Cdc20 )

	    	# Which variables change? del...=0 if no change.
		delCdh1 = Cdh1new - Cdh1
		delSBF  = SBFnew - SBF
		delCln2 = Cln2new - Cln2
		delClb5 = Clb5new - Clb5
		delClb2G = Clb2Gnew - Clb2G
		delClb2M = Clb2Mnew - Clb2M
		delCdc20 = Cdc20new - Cdc20
	
	    	# Propensity
	    	x1 = abs(delCdh1)*pCdh1
		x2 = x1 + abs(delSBF)*pSBF
		x3 = x2 + abs(delCln2)*pCln2
		x4 = x3 + abs(delClb5)*pClb5
		x5 = x4 + abs(delClb2G)*pClb2G
		x6 = x5 + abs(delClb2M)*pClb2M
		x7 = x6 + abs(delCdc20)*pCdc20
	
		# If no variables change, x7=0 so set propensity = pG1
	    	if(x7>0) { x8=x7 } else 
	    	{ x8=pG1 }
	
		# Determine which variable actually changes?
	    	s2 = runif(1)*x8
		y1 = (s2<x1)
		y2 = (s2>=x1)&&(s2<x2)
		y3 = (s2>=x2)&&(s2<x3)
		y4 = (s2>=x3)&&(s2<x4)
		y5 = (s2>=x4)&&(s2<x5)
		y6 = (s2>=x5)&&(s2<x6)
		y7 = (s2>=x6)&&(s2<x7)
		y8 = (s2>=x7)
		
		# Update real time given the probability per unit time of the possible changes
		delt = -log(runif(1))/x8
		
		# delt during Clb2M or Cdc20 activation are drawn from a lognormal distribution
		# This will overwrite above delt
		if(y6==1 && Clb2M==0)
		{
			delt <- rlnorm(1,log(tM.mean), tM.CV)
		}
		if(y7==1 && Cdc20==0)
		{
			delt <- rlnorm(1,log(tM.mean), tM.CV)
		}
		
		tr = tr + delt
		
		# Before updating vairables, check if tr hits the tick    
	    	Size -> Size.tick
	    	while(tr>tick) 
	    	{
			tick_timecourse <- c(tick_timecourse, tick)
	    		tr_timecourse <- c(tr_timecourse, tr)
			Cdh1_timecourse <- c(Cdh1_timecourse, Cdh1)
			SBF_timecourse <- c(SBF_timecourse, SBF)
			Cln2_timecourse <- c(Cln2_timecourse, Cln2)
			Clb5_timecourse <- c(Clb5_timecourse, Clb5)
			Clb2G_timecourse <- c(Clb2G_timecourse, Clb2G)
			Clb2M_timecourse <- c(Clb2M_timecourse, Clb2M)
			Cdc20_timecourse <- c(Cdc20_timecourse, Cdc20)
			Size_timecourse <- c(Size_timecourse, Size.tick)
			
			tick <- tick+1
			Size.tick = Size.tick*exp(mu*1)
			if(tick > tick_max) 
			{
				#print(tick); 
				break;
			}
	    	}
	    
		
	  	# cell division event
	    	if (Clb2G==1 && Clb2Gnew==0 && y5==1)
	    	{
	        	f = rlnorm(1, log(f.mean), f.CV)
	        
	        	# record initial conditions of the new daughter cell
	        	Size.D = Size*exp(mu*delt)*(1-f)
	        	initials.D <- rbind(initials.D, c(Size.D, tr))
	        
	        	# size at birth of the mother cell
	        	Size = Size*exp(mu*delt)*f

	        	S0=rlnorm(1,log(S0.mean), S0.CV) ## set a new S0
	    	} else
	    	{
	        	Size = Size*exp(mu*delt)
	    	}	
	
	    	# Update variables
	    	Cdh1 = Cdh1 + y1*delCdh1
		SBF  = SBF  + y2*delSBF
		Cln2 = Cln2 + y3*delCln2
		Clb5 = Clb5 + y4*delClb5
		Clb2G = Clb2G + y5*delClb2G
		Clb2M = Clb2M + y6*delClb2M
		Cdc20 = Cdc20 + y7*delCdc20
		  
	    	new.state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
		colnames(new.state)[10] <- 'Phase'
	    	state <- rbind(state, new.state)	
	    
	    	if(tick > tick_max) {break;} ## break from current repeat.
	    
	} ## end single trajectory
	
	Cdh1_array <- rbind(Cdh1_array, Cdh1_timecourse)
	SBF_array <- rbind(SBF_array, SBF_timecourse)
	Cln2_array <- rbind(Cln2_array, Cln2_timecourse)
	Clb5_array <- rbind(Clb5_array, Clb5_timecourse)
	Clb2G_array <- rbind(Clb2G_array, Clb2G_timecourse)
	Clb2M_array <- rbind(Clb2M_array, Clb2M_timecourse)
	Cdc20_array <- rbind(Cdc20_array, Cdc20_timecourse)
	Size_array <- rbind(Size_array, Size_timecourse)
	
} ## end repeat

colnames(initials.D) <- c('Size', 't.record')

state -> state1 # keep the state of a mother cell for plotting purpose

####################################
## PART 2 = Follow daughter cells ##
####################################

iD=0
while(1)
{	
	iD <- iD+1
	if(iD > dim(initials.D)[1]) { break }
	
	Size <- initials.D[iD,1]
	tr <- initials.D[iD,2]
	if(tr>399) { next; } # Skip dauther cells born after tr>399 (otherwise it will cause an array dimension problem)
	
	S0=rlnorm(1,log(S0.mean), S0.CV)
	
	# Initial condition
	Cdh1=1
	SBF=0
	Cln2=0
	Clb5=0
	Clb2G=0
	Clb2M=0
	Cdc20=1
			
	tick <- ceiling(tr)
	tick_timecourse <- c( 0:( ceiling(tr)-1 ) )
	tr_timecourse <- c( rep( NA, ceiling(tr) ) )
	Cdh1_timecourse <- c( rep( NA, ceiling(tr) ) )
	SBF_timecourse <- c( rep( NA, ceiling(tr) ) )
	Cln2_timecourse <- c( rep( NA, ceiling(tr) ) )
	Clb5_timecourse <- c( rep( NA, ceiling(tr) ) )
	Clb2G_timecourse <- c( rep( NA, ceiling(tr) ) )
	Clb2M_timecourse <- c( rep( NA, ceiling(tr) ) )
	Cdc20_timecourse <- c( rep( NA, ceiling(tr) ) )
	Size_timecourse <- c( rep( NA, ceiling(tr) ) )
	
	state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
	colnames(state)[10] <- 'Phase'
	
		
	while(1) # will break out when tick > tick_max
	{
	    	# Calculate new values by Boolean functions
	    	Cdh1new = !( Cln2 || Clb5 || Clb2G || Clb2M ) || ( Cdc20 && !( Cln2 || Clb5 || Clb2M ) )
		SBFnew = ( SBF && !( Clb2G || Clb2M ) ) || 
			( ( Cdh1 && !SBF && !Cln2 && !Clb5 && !Clb2G && !Clb2M && !Cdc20 ) &&
			  ( Size > S0 ) && ( runif(1) < (Size-S0)^2 )
			)
		Cln2new = SBF
		Clb5new = ( Clb5 || SBF ) && !( Cdh1 || Cdc20 )
		Clb2Gnew = Clb2M || ( ( Clb5 || Clb2G ) && !Cdh1 )
		Clb2Mnew = ( Clb2G || Clb2M ) && !Cdh1 && ( !Cdc20 || Clb5 ) && !Cln2
		Cdc20new = Clb2M || ( Clb2G && Cdc20 )

	    	# Which variables change? del...=0 if no change.
		delCdh1 = Cdh1new - Cdh1
		delSBF  = SBFnew - SBF
		delCln2 = Cln2new - Cln2
		delClb5 = Clb5new - Clb5
		delClb2G = Clb2Gnew - Clb2G
		delClb2M = Clb2Mnew - Clb2M
		delCdc20 = Cdc20new - Cdc20
	
	    	# Propensity
	    	x1 = abs(delCdh1)*pCdh1
		x2 = x1 + abs(delSBF)*pSBF
		x3 = x2 + abs(delCln2)*pCln2
		x4 = x3 + abs(delClb5)*pClb5
		x5 = x4 + abs(delClb2G)*pClb2G
		x6 = x5 + abs(delClb2M)*pClb2M
		x7 = x6 + abs(delCdc20)*pCdc20
	
		# If no variables change, x7=0 so set propensity = pG1
	    	if(x7>0) { x8=x7 } else 
	    	{ x8=pG1 }
	
		# Determine which variable actually changes?
	    	s2 = runif(1)*x8
		y1 = (s2<x1)
		y2 = (s2>=x1)&&(s2<x2)
		y3 = (s2>=x2)&&(s2<x3)
		y4 = (s2>=x3)&&(s2<x4)
		y5 = (s2>=x4)&&(s2<x5)
		y6 = (s2>=x5)&&(s2<x6)
		y7 = (s2>=x6)&&(s2<x7)
		y8 = (s2>=x7)
		
		# Update real time given the probability per unit time of the possible changes
		delt = -log(runif(1))/x8
		
		# delt during Clb2M and Cdc20 activation are drawn from a lognormal distribution
		# This will overwrite above delt
		if(y6==1 && Clb2M==0)
		{
			delt <- rlnorm(1,log(tM.mean), tM.CV)
		}
		if(y7==1 && Cdc20==0)
		{
			delt <- rlnorm(1,log(tM.mean), tM.CV)
		}
		
		tr = tr + delt
		
		# Before updating vairables, check if tr hits the tick    
		Size -> Size.tick
	    	while(tr>tick) 
	    	{
			tick_timecourse <- c(tick_timecourse, tick)
	    		tr_timecourse <- c(tr_timecourse, tr)
			Cdh1_timecourse <- c(Cdh1_timecourse, Cdh1)
			SBF_timecourse <- c(SBF_timecourse, SBF)
			Cln2_timecourse <- c(Cln2_timecourse, Cln2)
			Clb5_timecourse <- c(Clb5_timecourse, Clb5)
			Clb2G_timecourse <- c(Clb2G_timecourse, Clb2G)
			Clb2M_timecourse <- c(Clb2M_timecourse, Clb2M)
			Cdc20_timecourse <- c(Cdc20_timecourse, Cdc20)
			Size_timecourse <- c(Size_timecourse, Size.tick)
			
			tick <- tick+1
			Size.tick = Size.tick*exp(mu*1)
			if(tick > tick_max) 
			{
				#print(tick); 
				break;
			}
	    	}

	    	if (Clb2G==1 && Clb2Gnew==0 && y5==1)
	    	{
	        	f = rlnorm(1, log(f.mean), f.CV)
	        
	        	# record initial conditions of the new daughter cell
	        	Size.D = Size*exp(mu*delt)*(1-f)
	        	initials.D <- rbind(initials.D, c(Size.D, tr))
	        
	        	# size at birth of the mother cell
	        	Size = Size*exp(mu*delt)*f

	        	## Set a new S0
	        	S0=rlnorm(1,log(S0.mean), S0.CV)
	    	} else
	    	{
	        	Size = Size*exp(mu*delt)
	    	}	
	
	    	# Update variables
	    	Cdh1 = Cdh1 + y1*delCdh1
		SBF  = SBF  + y2*delSBF
		Cln2 = Cln2 + y3*delCln2
		Clb5 = Clb5 + y4*delClb5
		Clb2G = Clb2G + y5*delClb2G
		Clb2M = Clb2M + y6*delClb2M
		Cdc20 = Cdc20 + y7*delCdc20
		  
	    	new.state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
		colnames(new.state)[10] <- 'Phase'
	    	state <- rbind(state, new.state)		    

		if(tick > tick_max) {break;} ## break from current repeat.
		
	} ## end single trajectory
	
	Cdh1_array <- rbind(Cdh1_array, Cdh1_timecourse)
	SBF_array <- rbind(SBF_array, SBF_timecourse)
	Cln2_array <- rbind(Cln2_array, Cln2_timecourse)
	Clb5_array <- rbind(Clb5_array, Clb5_timecourse)
	Clb2G_array <- rbind(Clb2G_array, Clb2G_timecourse)
	Clb2M_array <- rbind(Clb2M_array, Clb2M_timecourse)
	Cdc20_array <- rbind(Cdc20_array, Cdc20_timecourse)
	Size_array <- rbind(Size_array, Size_timecourse)	
} ## end following all daughter cells from the list initials.D


#######################
## PART 3 = Plotting ##
#######################

pdf(width=8, height=8, file='Fig2.pdf')
par(mfrow = c(2,2), mar=c(4, 4, 2, 2)) #it goes c(bottom, left, top, right)

########################################
## plot single-cell dynamics
########################################

plot(state1$tr, state1$Cdh1, type='l', lwd=2, ylim=c(0,2.5), xlim=c(0,300), xlab='Time (min)', ylab='Boolean value')
lines(state1$tr, state1$Cln2+0.25, lwd=1, col="red")
lines(state1$tr, state1$Clb2G+state1$Clb2M, lwd=1, col="blue")
lines(state1$tr, state1$Cdc20+0.5, lwd=1, col="purple")
par(xpd=T)
text(-60, 2.8, expression(bold(a)), cex=1.5)
par(xpd=F)

legend('topright',
	legend=c(expression(paste(Cdh1)), expression(paste(Cln2+0.25)), expression(paste(Clb2[G]+Clb2[M])), expression(paste(Cdc20+0.5))),
	lty=c(1,1,1,1),
	lwd=c(2,1,1,1),
	col=c('black', 'red', 'blue', 'purple'),
	xjust=0, yjust=1,
	ncol=2
)
plot(state1$tr, state1$Size, xlim=c(0,300), type='l', lwd=2, ylim=c(0,1.5), xlab='Time (min)', ylab='Size')
par(xpd=T)
text(-60, 1.68, expression(bold(b)), cex=1.5)
par(xpd=F)

########################################
## plot population-level dynamics
########################################

plot(0:tick_max, colMeans(Cdh1_array, na.rm=T), type='l', lwd=2, ylim=c(0,2.5), xlim=c(0,300), xlab='Time (min)', ylab='Average Boolean value')
lines(0:tick_max, colMeans(Cln2_array, na.rm=T), col='red', lwd=1)
lines(0:tick_max, colMeans(Clb2G_array+Clb2M_array, na.rm=T), col='blue', lwd=1)
lines(0:tick_max, colMeans(Cdc20_array, na.rm=T), col='purple', lwd=1, lty=1)
par(xpd=T)
text(-60, 2.8, expression(bold(c)), cex=1.5)
par(xpd=F)

legend('topright',
	legend=c(expression(paste(Cdh1)), expression(paste(Cln2)), expression(paste(Clb2[G]+Clb2[M])), expression(paste(Cdc20))),
	lty=c(1,1,1,1),
	lwd=c(2,1,1,1),
	col=c('black', 'red', 'blue', 'purple'),
	xjust=0, yjust=1,
	ncol=2
)

plot(0:tick_max, colMeans(Size_array, na.rm=T), type='l', lwd=2, ylim=c(0,1.5), xlim=c(0,300), xlab='Time (min)', ylab='Average size')
par(xpd=T)
text(-60, 1.68, expression(bold(d)), cex=1.5)
par(xpd=F)

dev.off()
options(warn=0) # turn off warnings treated as errors
