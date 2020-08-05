

source("Acti Functions.R")


act_29pt <- read.csv("act_29pt.csv")
act_8pt <- read.csv("act_8pt.csv")
clinic_29pt_ahi <- read.csv("clinic_29pt_ahi.csv")
clinic_29pt_bmi <- read.csv("clinic_29pt_bmi.csv")
clinic_8pt <- read.csv("clinic_8pt.csv")
weekday <- read.csv("weekday.csv")




### fda.matchid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(act_29pt)
data(clinic_29pt_ahi)

data <- act_29pt
ahi <- clinic_29pt_ahi

### Example 1: Continuous Covariate
matchida  <- fda.matchid(data, ahi, "contin")


### Example 2: Categorical Covariate
ahi$ahicat <- as.factor(ifelse(ahi$AHI >= 0 & ahi$AHI <= 5, 1, 
				ifelse(ahi$AHI > 5 & ahi$AHI <= 15, 2,
						ifelse(ahi$AHI > 15 & ahi$AHI <= 30, 3,
								ifelse(ahi$AHI > 30, 4, 0)))))

matchidb  <- fda.matchid(data, ahi[,-2], "factor", 
		c("normal", "mild", "moderate", "severe"))


### fda.smoothdata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(act_29pt)
data(clinic_29pt_ahi)

data <- act_29pt
ahi <- clinic_29pt_ahi

matchid  <- fda.matchid(data, ahi, "contin")
FDcont <- fda.smoothdata(matchid)

### Smooth the Results
ts.plot(predict(FDcont$fd$fd, 1:1440), main="Smoothed Activity Data")


### flm_cate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(act_29pt)
data(clinic_29pt_ahi)

data <- act_29pt
ahi <- clinic_29pt_ahi

matchid  <- fda.matchid(data, ahi, "contin")
FDcont <- fda.smoothdata(matchid)

geftFDcont <- flm_cate(FDcont)


### cont_flm_plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(act_29pt)
data(clinic_29pt_ahi)

data <- act_29pt
ahi <- clinic_29pt_ahi

matchid  <- fda.matchid(data, ahi, "contin")
FDcont <- fda.smoothdata(matchid)

L <- nrow(data)
lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
xat <- c(0, L/4, L/2, 3*L/4, L)

geftFDcont <- flm_cate(FDcont)
predy <- as.vector(geftFDcont$freg$yhatfdobj$y)

xlim <- c(0, L) 
ylim <-  c(min(predy), max(predy) + 100)

legendx <- 0
legendy <- max(predy) - 100

cont.flm.results <- cont_flm_plot(FDcont, matchid, geftFDcont, xlim, 
		ylim, TRUE, 10, lb, xat, legendx, legendy, L)


### cat_flm_plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(act_29pt)
data(clinic_29pt_ahi)

data <- act_29pt
ahi <- clinic_29pt_ahi

ahi$ahicat <- as.factor(ifelse(ahi$AHI >= 0 & ahi$AHI <= 5, 1, 
				ifelse(ahi$AHI > 5 & ahi$AHI <= 15, 2,
				ifelse(ahi$AHI > 15 & ahi$AHI <= 30, 3,
				ifelse(ahi$AHI > 30, 4, 0)))))

matchidb <- fda.matchid(data, ahi[,-2] , "factor", 
		c("normal", "mild", "moderate", "severe"))
FDcatahi <- fda.smoothdata(matchidb)

L <- nrow(data)
lb <- c("Midnight", "6AM", "Noon", "6PM", "Midnight") 
xat <- c(0, L/4, L/2, 3*L/4, L)

geftFDcatahi <- flm_cate(FDcatahi)
predy <- as.vector(geftFDcatahi$freg$yhatfdobj$y)

xlim <- c(0, L) 
ylim <- c(min(predy), max(predy) + 100)

cat.flm.results <- cat_flm_plot(FDcatahi, matchidb, geftFDcatahi, 
		TRUE, 5, lb, xat, "AHI", 1:4, ylim, L)

