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

#########################################################
## PART 1 = Generate initial mother and daughter cells ##
#########################################################

collect <- c()

for(type in c('M','D'))
{
	tr_max=2000;
	sim.who=type
		
	j=65 ## j=64 is for state '1000000' (G1); j=65 is for state '1000001' (cell at birth)
		
	for(rep in 1:10)
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
		Size=0.65;
		tr=0;
		
		state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
		colnames(state)[10] <- 'Phase'
			 
		S0=rlnorm(1,log(S0.mean), S0.CV)
		
		while(1) # simulation step loop # will break out when tr >= tr_max
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
			
			# cell division event
			if (Clb2G==1 && Clb2Gnew==0 && y5==1)
		    {
		        if(sim.who=='M')
		        {
		        	f = rlnorm(1, log(f.mean), f.CV)
		        } else
		        {
		        	f = rlnorm(1, log(1-f.mean), f.CV)
		        }
		        
		        Size = Size*exp(mu*delt)*f
		        S0=rlnorm(1,log(S0.mean), S0.CV) ## set a new S0
		        
		    } else
		    {
		        Size = Size*exp(mu*delt)
		    }
		
		    # Update variables
		    Cdh1 = Cdh1 + y1*delCdh1
			SBF = SBF  + y2*delSBF
			Cln2 = Cln2 + y3*delCln2
			Clb5 = Clb5 + y4*delClb5
			Clb2G = Clb2G + y5*delClb2G
			Clb2M = Clb2M + y6*delClb2M
			Cdc20 = Cdc20 + y7*delCdc20
		    
		    new.state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
			colnames(new.state)[10] <- 'Phase'
		    state <- rbind(state, new.state)
		
			if( tr >= tr_max)
			{
				collect <- rbind(collect, c(Size, 0, 0, 0, paste0('s', Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20), sim.who, 0))
		    	break;
			}	    
		} ## single trajectory		
	} ## end repeat
} ## end simulating mother or daughter cells

colnames(collect) <- c('Size', 'last.tdiv', 'last.TG1', 'last.Size_at_birth', 'initial.state', 'D or M', 't.record')
write.csv(collect, file='initials.csv', row.names=F)

# dev.new(width=8, height=4)
# par(mfrow=c(1,2))
# hist(collect[,1], main='Daughter cells', xlab='Initial sizes')
# hist(collect[,1], main='Mother cells', xlab='Initial sizes')

#####################################################
## PART 2 = Track all mother and daughter lineages ##
#####################################################

initials <- read.csv("initials.csv", header=T, stringsAsFactors=F)

tr_max=600;
Size_at_birth.M <- c()
Division_time.M <- c()
TG1_time.M <- c()
Tb_time.M <- c()
Size_at_birth.D <- c()
Division_time.D <- c()
TG1_time.D <- c()
Tb_time.D <- c()
Joint_dist.M <- c()
Joint_dist.D <- c()

ID=0
while(1)
{	
	ID=ID+1
	if(ID > dim(initials)[1]) { break }
	
	Size= as.numeric(initials$Size[ID])
	last.tdiv= as.numeric(initials$last.tdiv[ID])
	last.TG1= as.numeric(initials$last.TG1[ID])
	last.Size_at_birth= as.numeric(initials$last.Size_at_birth[ID])
	tr=as.numeric(initials$t.record[ID])
	sim.who=initials$D.or.M[ID]
	strsplit(initials$initial.state[ID], "") -> initial.State
	
	Cdh1=as.numeric(unlist(initial.State)[2])
	SBF=as.numeric(unlist(initial.State)[3])
	Cln2=as.numeric(unlist(initial.State)[4])
	Clb5=as.numeric(unlist(initial.State)[5])
	Clb2G=as.numeric(unlist(initial.State)[6])
	Clb2M=as.numeric(unlist(initial.State)[7])
	Cdc20=as.numeric(unlist(initial.State)[8])
	
	state <- data.frame(tr, Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20, Size, paste0(Cdh1, SBF, Cln2, Clb5, Clb2G, Clb2M, Cdc20))
	colnames(state)[10] <- 'Phase'

	S0=rlnorm(1,log(S0.mean), S0.CV)
	
	while(1)
	{
	    # Calculate new values by Boolean functions
	    Cdh1new = !( Cln2 || Clb5 || Clb2G || Clb2M ) || ( Cdc20 && !( Cln2 || Clb5 || Clb2M ) )
		SBFnew = ( SBF && !( Clb2G || Clb2M ) ) || 
			( 
				( Cdh1 && !SBF && !Cln2 && !Clb5 && !Clb2G && !Clb2M && !Cdc20 ) &&
			  	( Size	 > S0 ) && ( runif(1) < (Size-S0)^2 )
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
	    if(x7>0) {  x8=x7 } else 
	    {  x8=pG1 }
	
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
		
		# START event
		if( SBF==0 && y2==1 )
		{
			if(last.tdiv>0)
			{
				if(sim.who=='M')
				{
					TG1_time.M <- c(TG1_time.M, tr-last.tdiv)
					Joint_dist.M <- rbind(Joint_dist.M, c(last.Size_at_birth, tr-last.tdiv))
				} else if(sim.who=='D')
				{
					TG1_time.D <- c(TG1_time.D, tr-last.tdiv)
					Joint_dist.D <- rbind(Joint_dist.D, c(last.Size_at_birth, tr-last.tdiv))
				} else { stop('Error #1') }
				last.TG1 <- tr
			}
		}
		 
		# cell division event
	    if (Clb2G==1 && Clb2Gnew==0 && y5==1)
	    {
	    	f = rlnorm(1, log(f.mean), f.CV)
	        Size.D = Size*exp(mu*delt)*(1-f)
	        Size = Size*exp(mu*delt)*f
	        
	        Size_at_birth.D <- c(Size_at_birth.D, Size.D)
	        Size_at_birth.M <- c(Size_at_birth.M, Size)
	        
	        last.Size_at_birth=Size
	        
	        if(last.tdiv>0)
			{
				if(sim.who=='M')
				{
					Division_time.M <- c(Division_time.M, tr-last.tdiv)
				} else if(sim.who=='D') 
				{
					Division_time.D <- c(Division_time.D, tr-last.tdiv)
				} else { stop('Error #2') }
			}
			last.tdiv <- tr
			
			if( last.TG1 >0 )
			{
				if(sim.who=='M')
				{
					Tb_time.M <- c(Tb_time.M, tr-last.TG1)	
				} else if(sim.who=='D')
				{
					Tb_time.D <- c(Tb_time.D, tr-last.TG1)	
				} else { stop('Error #3') }
			}
	        
	        S0=rlnorm(1,log(S0.mean), S0.CV) ## set a new S0
	        
	        initials <- rbind(initials, c(Size.D, last.tdiv, last.TG1, Size.D, 's1000001', 'D', tr))
	        
	        sim.who='M' # if it's a daughter cell, it becomes a mother cell
	        
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
	    
	    if(tr>=tr_max) { break };
	    
	} ## end a single trajectory
} ## end all cells


#################################
## PART 3 Plot stat comparison ##
#################################

stat_exp <- read.csv('stat_exp.csv', header=F)
## calculate fL per 1 size unit
fl_ratio = stat_exp[6,1] / mean(Size_at_birth.M)

#dev.new(width=8, height=4)
pdf(file='Fig3 (Barplots).pdf', width=8, height=8)
par(mfrow=c(2,2))
m <- matrix(c( stat_exp[1,1], mean(Division_time.M),
				stat_exp[4,1], mean(TG1_time.M),
				stat_exp[5,1], mean(Tb_time.M),
				stat_exp[6,1], mean(Size_at_birth.M)*fl_ratio ),
				nrow=2, ncol=4, byrow=F )
barplot(m, beside=T, names.arg=c(expression(paste(italic(T)[c])), expression(paste(italic(T)[G1])), expression(paste(italic(T)[bud])), 'Size at birth'),
	main='Mother cells', ylim=c(0,120), ylab='Mean value (min or fL)', col=c('blue3','cadetblue3'))
par(xpd=T)
text(-1.5, 1.2*120, expression(bold(a)), cex=1.5)
par(xpd=F)

legend(
	"topright",
	legend=c('Experiment','Simulation'),
	fill=c('blue3', 'cadetblue3')
)

m <- matrix(c( stat_exp[7,1], mean(Division_time.D),
				stat_exp[10,1], mean(TG1_time.D),
				stat_exp[11,1], mean(Tb_time.D),
				stat_exp[12,1], mean(Size_at_birth.D)*fl_ratio ),
				nrow=2, ncol=4, byrow=F )
barplot(m, beside=T, names.arg=c(expression(paste(italic(T)[c])), expression(paste(italic(T)[G1])), expression(paste(italic(T)[bud])), 'Size at birth'),
	main='Daughter cells', ylim=c(0,120), ylab='Mean value (min or fL)', col=c('blue3','cadetblue3'))
par(xpd=T)
text(-1.5, 1.2*120, expression(bold(b)), cex=1.5)
par(xpd=F)


legend(
	"topright",
	legend=c('Experiment','Simulation'),
	fill=c('blue3', 'cadetblue3')
)

m <- matrix(c( stat_exp[1,2], sd(Division_time.M)/mean(Division_time.M),
				stat_exp[4,2], sd(TG1_time.M)/mean(TG1_time.M),
				stat_exp[5,2], sd(Tb_time.M)/mean(Tb_time.M),
				stat_exp[6,2], sd(Size_at_birth.M)/mean(Size_at_birth.M)),
				nrow=2, ncol=4, byrow=F )
barplot(m, beside=T, names.arg=c(expression(paste(italic(T)[c])), expression(paste(italic(T)[G1])), expression(paste(italic(T)[bud])), 'Size at birth'),
	main='Mother cells', ylim=c(0,1.0), ylab='CV', col=c('blue3','cadetblue3')) # , legend.text=c('Exp', 'Sim')
par(xpd=T)
text(-1.5, 1.2*1.0, expression(bold(c)), cex=1.5)
par(xpd=F)


legend(
	"topright",
	legend=c('Experiment','Simulation'),
	fill=c('blue3', 'cadetblue3')
)

m <- matrix(c( stat_exp[7,2], sd(Division_time.D)/mean(Division_time.D),
				stat_exp[10,2], sd(TG1_time.D)/mean(TG1_time.D),
				stat_exp[11,2], sd(Tb_time.D)/mean(Tb_time.D),
				stat_exp[12,2], sd(Size_at_birth.D)/mean(Size_at_birth.D)),
				nrow=2, ncol=4, byrow=F )
barplot(m, beside=T, names.arg=c(expression(paste(italic(T)[c])), expression(paste(italic(T)[G1])), expression(paste(italic(T)[bud])), 'Size at birth'),
	main='Daughter cells', ylim=c(0,1.0), ylab='CV', col=c('blue3','cadetblue3')) # , , legend.text=c('Exp', 'Sim')
par(xpd=T)
text(-1.5, 1.2*1.0, expression(bold(d)), cex=1.5)
par(xpd=F)

legend(
	"topright",
	legend=c('Experiment','Simulation'),
	fill=c('blue3', 'cadetblue3')
)

dev.off()

#####################################
## PART 4 Plot joint distributions ##
#####################################

Joint_dist.M -> Joint_dist.M_original
Joint_dist.D -> Joint_dist.D_original

# Sampling only some cells
Joint_dist.M[sample( dim(Joint_dist.M)[1], 80, replace=F ),] -> Joint_dist.M
Joint_dist.D[sample( dim(Joint_dist.D)[1], 80, replace=F ),] -> Joint_dist.D

# dev.new(width=8, height=4)
pdf(width=8, height=4, file='Fig4 joint distribution.pdf')
par(mfrow=c(1,2))

##################
## Mother cells ##
##################

Joint_dist.M[,1]*fl_ratio -> Size_in_fL.M

plot(log(Size_in_fL.M/mean(Size_in_fL.M)), mu*Joint_dist.M[,2], xlab='ln(Size at birth)', ylab=expression(paste(mu,italic(T)[G1])), main='Mother cells', xlim=c(-1,0.6), ylim=c(-0.1,1), cex.axis=0.9)

bin.range <- seq(0,100,2)
bin.data <- c()
for( i in 1:(length(bin.range)-1) )
{
  #print(bin.range[i])
  which( Size_in_fL.M >= bin.range[i] ) -> index1
  which( Size_in_fL.M < bin.range[i+1] ) -> index2
  intersect(index1, index2) -> index3
  bin.data <- rbind( bin.data, c( (bin.range[i]+bin.range[i+1])/2,
                                  mean(Joint_dist.M[index3,2])
                                 ) )
} 
points(log(bin.data[,1]/mean(Size_in_fL.M)), mu*bin.data[,2], cex=1.8, col='red', pch=16)
log(bin.data[,1]/mean(Size_in_fL.M)) -> x
mu*bin.data[,2] -> y
lm(y~x) -> fit
x <- seq(-2,2,0.5)
y = fit$coefficients[2]*x +fit$coefficients[1]
lines(x,y, xpd=F, lwd=2, col="red")
text(-0.5,0.6, paste0('slope=', sprintf("%.4f",fit$coefficients[2]) ))
par(xpd=T)
text(-1.5, 1.4, expression(bold(a)), cex=1.5)
par(xpd=F)

####################
## Daughter cells ##
####################

Joint_dist.D[,1]*fl_ratio -> Size_in_fL.D

plot(log(Size_in_fL.D/mean(Size_in_fL.M)), mu*Joint_dist.D[,2], xlab='ln(Size at birth)', ylab=expression(paste(mu,italic(T)[G1])), main='Daughter cells', xlim=c(-1.0,0.6), ylim=c(-0.1,1), cex.axis=0.9)

bin.range <- seq(0,100,2)
bin.data <- c()
for( i in 1:(length(bin.range)-1) )
{
  #print(bin.range[i])
  which( Size_in_fL.D >= bin.range[i] ) -> index1
  which( Size_in_fL.D < bin.range[i+1] ) -> index2
  intersect(index1, index2) -> index3
  bin.data <- rbind( bin.data, c( (bin.range[i]+bin.range[i+1])/2,
                                  mean(Joint_dist.D[index3,2])
                                 ) )
} 
points(log(bin.data[,1]/mean(Size_in_fL.M)), mu*bin.data[,2], cex=1.8, col='red', pch=16)

which( log(bin.data[,1]/mean(Size_in_fL.M)) <= -0.4 ) -> index.small
which( log(bin.data[,1]/mean(Size_in_fL.M)) > -0.4 ) -> index.large

log(bin.data[index.small,1]/mean(Size_in_fL.M)) -> x
mu*bin.data[index.small,2] -> y
lm(y~x) -> fit
x <- seq(-2,-0.4,0.1)
y = fit$coefficients[2]*x +fit$coefficients[1]
lines(x,y, xpd=F, lwd=2, col="red")
text(-0.7,0.9, paste0('slope=', sprintf("%.4f",fit$coefficients[2]) ))
par(xpd=F)

log(bin.data[index.large,1]/mean(Size_in_fL.M)) -> x
mu*bin.data[index.large,2] -> y
lm(y~x) -> fit
x <- seq(-0.4,2,0.1)
y = fit$coefficients[2]*x +fit$coefficients[1]
lines(x,y, xpd=F, lwd=2, col="red")
text(0.3,0.3, paste0('slope=', sprintf("%.4f",fit$coefficients[2]) ))
par(xpd=T)
text(-1.5, 1.4, expression(bold(b)), cex=1.5)
par(xpd=F)


dev.off()

options(warn=0) # turn off warnings treated as errors
