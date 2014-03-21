#structure run analysis
library(ggplot2)

#first run the parse_struc.py script

#load chain
chain_k2 = read.table("chains_K2_sum.txt",header=T,na.strings='-')

#plot alpha chains across Reps
ggplot(chain_k2,aes(x=Step,y=Alpha,col=Rep))+geom_line()

#check alpha histograms
ggplot(chain_k2,aes(x=Alpha))+geom_histogram()+facet_grid(Rep~.)

#plot LnLike chain
ggplot(chain_k2,aes(x=Step,y=Ln_Like,col=Rep))+geom_line()

#check Ln_Like histograms
ggplot(chain_k2,aes(x=Ln_Like))+geom_histogram()+facet_grid(Rep~.)

#check correlation between F
ggplot(chain_k2,aes(x=F1,y=F2))+geom_point()+coord_fixed(0.5)+facet_grid(Rep~.)

################################################################
# plotting evanno data after applying structureHarvester.py

evanno_res = read.table("../second_session/evanno.txt",header=F,comment.char='#')
names(evanno_res) = c("K","reps","mean_LnPK",	"sd_LnPK",	"Ln1K",	"Ln2K",	"Delta_K")

ggplot(evanno_res, aes(x=K, y=mean_LnPK)) + geom_errorbar(aes(ymin=mean_LnPK-sd_LnPK, ymax=mean_LnPK+sd_LnPK), width=.1) + geom_line() + geom_point()

ggplot(evanno_res, aes(x=K, y=Delta_K)) + geom_line() + geom_point()
