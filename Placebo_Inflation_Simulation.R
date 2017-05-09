# PLacebo Inflation Simulation
pla_effect=8 # Rounded effect from LDOPA STUDY
pla_sd=14 # Rounded effect from LDOPA STUDY
n=70 
n_rep=10000
ps=rep(NA, n_rep)
df=data.frame(ps)

for (i in 1:n_rep) {
sampl=rnorm(n,pla_effect,pla_sd)
df$ps[i]=t.test(sampl,mu=0)$p.value
}

length(df$ps[df$ps<0.05])/length(df$ps)
ggplot(df,aes(ps))+geom_histogram(binwidth = 0.05)

for (i in 1:n_rep) {
  sampl=rnorm(n,pla_effect,pla_sd)
  df$ps[i]=t.test(sampl,mu=0)$p.value
  df$ps[df$ps==min(df$ps)]=NA
}

length(df$ps[df$ps<0.05])/length(df$ps)
ggplot(df,aes(ps))+geom_histogram(binwidth = 0.05)