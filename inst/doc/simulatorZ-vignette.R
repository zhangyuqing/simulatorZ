### R code from vignette source 'simulatorZ-vignette.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: simulatorZ-vignette.Rnw:57-58
###################################################
library(simulatorZ)


###################################################
### code chunk number 3: simulatorZ-vignette.Rnw:93-97 (eval = FALSE)
###################################################
## library(curatedOvarianData)
## source(system.file("extdata",
##                    "patientselection.config",package="curatedOvarianData"))
## source(system.file("extdata", "createEsetList.R", package="curatedOvarianData"))


###################################################
### code chunk number 4: simulatorZ-vignette.Rnw:104-109
###################################################
library(curatedOvarianData)
data(GSE17260_eset)
data(E.MTAB.386_eset)
data(GSE14764_eset)
esets.orig <- list(GSE17260=GSE17260_eset, E.MTAB.386=E.MTAB.386_eset, GSE14764=GSE14764_eset)


###################################################
### code chunk number 5: cleanEsets
###################################################
cleanEsets <- function(obj, keep.common.only=TRUE, meta.required=NULL){
  if(keep.common.only){
    intersect.features <- Reduce(intersect, lapply(obj, featureNames))
    for (i in 1:length(obj))
      obj[[i]] <- obj[[i]][intersect.features, ]
  }
  if(!is.null(meta.required)){
    obj <- lapply(obj, function(obj1){
      for (x in meta.required)
        obj1 <- obj1[, !is.na(obj1[[x]])]
      if(ncol(obj1) > 0){
        return(obj1)
      }else{
        return(NULL)
      }
    })
    return(obj[!sapply(obj, is.null)])
  }
}


###################################################
### code chunk number 6: simulatorZ-vignette.Rnw:143-146
###################################################
esets <- cleanEsets(esets.orig, meta.required=c("days_to_death", "vital_status"))
esets <- lapply(esets, function(eset) eset[1:1000, ])
sapply(esets, dim)


###################################################
### code chunk number 7: simulatorZ-vignette.Rnw:158-161
###################################################
set.seed(8)
sim.esets <- simData(esets, n.samples=500)
names(sim.esets)


###################################################
### code chunk number 8: simulatorZ-vignette.Rnw:183-184
###################################################
sim.esets <- simData(esets, 500, type="one-step")


###################################################
### code chunk number 9: simulatorZ-vignette.Rnw:194-195
###################################################
cleaned.esets <- cleanEsets(esets.orig, meta.required="tumorstage")


###################################################
### code chunk number 10: simulatorZ-vignette.Rnw:199-203
###################################################
sim.sets <- simData(cleaned.esets, 500, 
                    balance.variables="tumorstage")
sim.sets$prob.desired
sim.sets$prob.real[[1]]


###################################################
### code chunk number 11: simulatorZ-vignette.Rnw:209-210
###################################################
cleaned.esets <- cleanEsets(esets.orig, meta.required=c("tumorstage", "grade"))


###################################################
### code chunk number 12: simulatorZ-vignette.Rnw:214-218
###################################################
sim.sets <- simData(cleaned.esets, 500, 
                    balance.variables=c("tumorstage", "grade"))
sim.sets$prob.desired
sim.sets$prob.real[[1]]


###################################################
### code chunk number 13: simulatorZ-vignette.Rnw:228-231
###################################################
y.list <- lapply(esets, function(eset){
  return( Surv(eset$days_to_death, eset$vital_status=="deceased") )
})


###################################################
### code chunk number 14: simulatorZ-vignette.Rnw:233-234
###################################################
true.mod <- getTrueModel(esets, y.list, 100)


###################################################
### code chunk number 15: simulatorZ-vignette.Rnw:236-237
###################################################
names(true.mod)


###################################################
### code chunk number 16: simulatorZ-vignette.Rnw:240-242
###################################################
simmodel <- simData(esets, 500, y.list)
new.esets <- simTime(simmodel, true.mod)


###################################################
### code chunk number 17: simulatorZ-vignette.Rnw:249-252
###################################################
true.mod <- getTrueModel(esets, y.list, 100)
sim.sets <- simData(esets, 500, y.list)
sim.sets <- simTime(sim.sets, true.mod)


###################################################
### code chunk number 18: simulatorZ-vignette.Rnw:256-257
###################################################
sim.sets <- simBootstrap(esets, y.list, 500, 100)


###################################################
### code chunk number 19: simulatorZ-vignette.Rnw:267-281
###################################################
y.vars=y.list
balance.variables=c("tumorstage", "grade")
X.list <- lapply(esets, function(eset){
  newx <- t(exprs(eset))
  return(newx)
}) 
all.X <- do.call(rbind, X.list)
cov.list <- lapply(esets, function(eset){
  return(pData(eset)[, balance.variables])
})  
all.cov <- as(data.frame(do.call(rbind, cov.list)), "AnnotatedDataFrame")
rownames(all.cov) <- colnames(t(all.X))
all.yvars <- do.call(rbind, y.vars)
whole_eset <- ExpressionSet(t(all.X), all.cov)


###################################################
### code chunk number 20: simulatorZ-vignette.Rnw:286-304
###################################################
lp.index <- c()
for(i in 1:length(esets)){
  lp.index <- c(lp.index, rep(i, length(y.vars[[i]][, 1])))
}

truemod <- getTrueModel(list(whole_eset), list(all.yvars), parstep=100)

simmodels <- simData(esets, 150, y.vars)
beta <- grid <- survH <- censH <- lp <- list()
for(listid in 1:length(esets)){
  beta[[listid]] <- truemod$beta[[1]]
  grid[[listid]] <- truemod$grid[[1]]
  survH[[listid]] <- truemod$survH[[1]]
  censH[[listid]] <- truemod$censH[[1]]
  lp[[listid]] <- truemod$lp[[1]][which(lp.index==listid)]
}
res <- list(beta=beta, grid=grid, survH=survH, censH=censH, lp=lp)
simmodels <- simTime(simmodels, res)


###################################################
### code chunk number 21: simulatorZ-vignette.Rnw:323-325
###################################################
library(parathyroidSE)
data("parathyroidGenesSE")


###################################################
### code chunk number 22: simulatorZ-vignette.Rnw:331-334
###################################################
sim.sets <- simData(list(parathyroidGenesSE), 100)
names(sim.sets)
sim.sets$obj[[1]]   #The simulated SummarizedExperiment


###################################################
### code chunk number 23: simulatorZ-vignette.Rnw:339-340
###################################################
table(colData(parathyroidGenesSE)$run)


###################################################
### code chunk number 24: simulatorZ-vignette.Rnw:345-346
###################################################
table(colData(sim.sets$obj[[1]])$run)


###################################################
### code chunk number 25: simulatorZ-vignette.Rnw:357-362
###################################################
tr.size <- 450
simmodel.tr <- simBootstrap(esets[1], y.list[1], tr.size, 100)
tr.set <- simmodel.tr$obj.list[[1]]
X.tr <- t(exprs(tr.set))
y.tr <- simmodel.tr$y.vars.list[[1]]


###################################################
### code chunk number 26: simulatorZ-vignette.Rnw:365-370
###################################################
val.size <- 450
simmodel.val <- simBootstrap(esets[1], y.list[1], val.size, 100)
val.set <- simmodel.val$obj.list[[1]]
X.val <- t(exprs(val.set))
y.val <- simmodel.val$y.vars.list[[1]]


###################################################
### code chunk number 27: simulatorZ-vignette.Rnw:373-378
###################################################
#check C-Index for true lp
val.par <- getTrueModel(esets, y.list, 100)
lpboot <- val.par$lp[[1]][simmodel.val$indices[[1]]]
library(Hmisc)
c.ind <- rcorr.cens(-lpboot, y.val)[1]


###################################################
### code chunk number 28: simulatorZ-vignette.Rnw:380-381
###################################################
print(c.ind)


###################################################
### code chunk number 29: simulatorZ-vignette.Rnw:386-394
###################################################
library(superpc)
tr.data<- data<-list(x=t(X.tr),y=y.tr[,1], censoring.status=y.tr[,2], 
                     featurenames=colnames(X.tr))
fit.tr<-superpc.train(data=tr.data,type='survival')
#tuning
cv.tr<-superpc.cv(fit.tr,data=tr.data)
n.comp<-which.max(apply(cv.tr$scor, 1, max, na.rm = TRUE))
thresh<-cv.tr$thresholds[which.max(cv.tr$scor[n.comp, ])]


###################################################
### code chunk number 30: superpcbox
###################################################
lp.tr<- superpc.predict(fit.tr, tr.data, tr.data, threshold=thresh, 
                        n.components=n.comp)$v.pred.1df
boxplot(lp.tr)


###################################################
### code chunk number 31: simulatorZ-vignette.Rnw:406-410
###################################################
data.val<- data<-list(x=t(X.val),y=y.val[,1], censoring.status=y.val[,2], 
                      featurenames=colnames(X.tr))
lp.val<-superpc.predict(fit.tr, tr.data, data.val, threshold=thresh, 
                        n.components=n.comp)$v.pred.1df


###################################################
### code chunk number 32: simulatorZ-vignette.Rnw:413-415
###################################################
print('C-Index')
(c.ind<-rcorr.cens(-lp.val,y.val)[1])


###################################################
### code chunk number 33: simulatorZ-vignette.Rnw:418-420
###################################################
print('correlation to true lp')
(corlps<-cor(lp.val,lpboot,method='pearson'))


###################################################
### code chunk number 34: simulatorZ-vignette.Rnw:427-429
###################################################
z <- zmatrix(obj=esets, y.vars=y.list,
             fold=3,trainingFun=plusMinus)


###################################################
### code chunk number 35: simulatorZ-vignette.Rnw:431-432
###################################################
print(z)


###################################################
### code chunk number 36: simulatorZ-vignette.Rnw:435-451
###################################################
Z.list <- list()
CV <- CSV <- c()
for(b in 1:20){
  print(paste("iteration: ", b, sep=""))
  sim2.esets <- simBootstrap(obj=esets, n.samples=150, y.vars=y.list,
                             parstep=100, type="two-steps")
  Z.list[[b]] <- zmatrix(obj=sim2.esets$obj.list, 
                         y.vars=sim2.esets$y.vars.list, fold=4,
                         trainingFun=plusMinus)
  sum.cv <- 0
  for(i in 1:length(esets)){
    sum.cv <- sum.cv + Z.list[[b]][i, i]
  }
  CV[b] <- sum.cv / length(esets)
  CSV[b] <- (sum(Z.list[[b]]) - sum.cv) / (length(esets)*(length(esets)-1))
}


###################################################
### code chunk number 37: box
###################################################
average.Z <- Z.list[[1]]
for(i in 2:length(Z.list)){
  average.Z <- average.Z + Z.list[[i]]
}
average.Z <- average.Z / 20
print(average.Z)

resultlist <- list(CSV=CSV, CV=CV)
boxplot(resultlist, col=c("white", "grey"), ylab="C-Index", 
        boxwex = 0.25, xlim=c(0.5, 2.5))


###################################################
### code chunk number 38: simulatorZ-vignette.Rnw:469-470
###################################################
sessionInfo()
