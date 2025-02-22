# model
exp.dacay.model <- nls(y ~ alpha * exp(beta * x) + theta , data = data.df, start = start)

# GFP expression
## model 1
df.modelI <- df %>% dplyr::filter(model == "I")
df.gfp.1 <- data.frame(x = df.modelI$Mean.live, y = df.modelI$CV2)

theta.gfp.I.0 <- min(df.gfp.1$y) * 0.5  
model.gfp.I.0 <- lm(log(y - theta.gfp.I.0) ~ x, data=df.gfp.1)  
alpha.gfp.I.0 <- exp(coef(model.gfp.I.0)[1])
beta.gfp.I.0 <- coef(model.gfp.I.0)[2]

start.gfp.1 <- list(alpha = alpha.gfp.I.0, beta = beta.gfp.I.0, theta = theta.gfp.I.0)

model.gfp.I <- nls(y ~ alpha * exp(beta * x) + theta , data = df.gfp.1, start = start.gfp.1)

## model 2
df.modelII <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "II")
df.gfp.2 <- data.frame(x = df.modelII$Mean.live, y = df.modelII$CV2)

theta.gfp.II.0 <- min(df.gfp.2$y) * 0.5  
model.gfp.II.0 <- lm(log(y - theta.gfp.II.0) ~ x, data=df.gfp.2)  
alpha.gfp.II.0 <- exp(coef(model.gfp.II.0)[1])
beta.gfp.II.0 <- coef(model.gfp.II.0)[2]

start.gfp.2 <- list(alpha = alpha.gfp.II.0, beta = beta.gfp.II.0, theta = theta.gfp.II.0)

model.gfp.II <- nls(y ~ alpha * exp(beta * x) + theta , data = df.gfp.2, start = start.gfp.2)

## model 3
df.modelIII <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "III")
df.gfp.3 <- data.frame(x = df.modelIII$Mean.live, y = df.modelIII$CV2)

theta.gfp.III.0 <- min(df.gfp.3$y) * 0.5  
model.gfp.III.0 <- lm(log(y - theta.gfp.III.0) ~ x, data=df.gfp.3)  
alpha.gfp.III.0 <- exp(coef(model.gfp.III.0)[1])
beta.gfp.III.0 <- coef(model.gfp.III.0)[2]

start.gfp.3 <- list(alpha = alpha.gfp.III.0, beta = beta.gfp.III.0, theta = theta.gfp.III.0)

model.gfp.III <- nls(y ~ alpha * exp(beta * x) + theta , data = df.gfp.3, start = list(beta = beta.gfp.III.0), alg = "plinear")

## model 4
df.modelIV <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "IV")
df.gfp.4 <- data.frame(x = df.modelIV$Mean.live, y = df.modelIV$CV2)

theta.gfp.IV.0 <- min(df.gfp.4$y) * 0.5  
model.gfp.IV.0 <- lm(log(y - theta.gfp.IV.0) ~ x, data=df.gfp.4)  
alpha.gfp.IV.0 <- exp(coef(model.gfp.IV.0)[1])
beta.gfp.IV.0 <- coef(model.gfp.IV.0)[2]

start.gfp.4 <- list(alpha = alpha.gfp.IV.0, beta = beta.gfp.IV.0, theta = theta.gfp.IV.0)

model.gfp.IV <- nls(y ~ alpha * exp(beta * x) + theta , data = df.gfp.4, start = list(beta = beta.gfp.IV.0, theta = theta.gfp.IV.0), alg = "plinear")

## model 5
df.modelV <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "V")
df.gfp.5 <- data.frame(x = df.modelV$Mean.live, y = df.modelV$CV2)

theta.gfp.V.0 <- min(df.gfp.5$y) * 0.5  
model.gfp.V.0 <- lm(log(y - theta.gfp.V.0) ~ x, data=df.gfp.5)  
alpha.gfp.V.0 <- exp(coef(model.gfp.V.0)[1])
beta.gfp.V.0 <- coef(model.gfp.V.0)[2]

start.gfp.5 <- list(alpha = alpha.gfp.V.0, beta = beta.gfp.V.0, theta = theta.gfp.V.0)

model.gfp.V <- nls(y ~ alpha * exp(beta * x) + theta , data = df.gfp.5, start = list(beta = beta.gfp.V.0), alg = "plinear")

ggplot(df.gfp, aes(x = x, y = y, color = model))+geom_point(size = 2.5)+geom_line(aes(x = x, y = predict(model.gfp.I, list(x = x))), data = df.gfp.1, color = "red")+geom_line(aes(x = x, y = predict(model.gfp.II, list(x = x))), data = df.gfp.2, color = "orange")+geom_line(aes(x = x, y = predict(model.gfp.III, list(x = x))), data = df.gfp.3, color = "#1589FF")+geom_line(aes(x = x, y = predict(model.gfp.IV, list(x = x))), data = df.gfp.4, color = "grey")+geom_line(aes(x = x, y = predict(model.gfp.V, list(x = x))), data = df.gfp.5, color = "#00CC33")+theme_bw()+scale_color_manual(values = c("red", "orange", "#1589FF", "grey", "#00CC33"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# sRNAs (qPCR)
## model 1
df.sRNAs.1 <- data.frame(x = df.modelI$qPCR_s, y = df.modelI$CV2)

theta.sRNAs.I.0 <- min(df.sRNAs.1$y) * 0.5  
model.sRNAs.I.0 <- lm(log(y - theta.sRNAs.I.0) ~ x, data=df.sRNAs.1)  
alpha.sRNAs.I.0 <- exp(coef(model.sRNAs.I.0)[1])
beta.sRNAs.I.0 <- coef(model.sRNAs.I.0)[2]

start.sRNAs.1 <- list(alpha = alpha.sRNAs.I.0, beta = beta.sRNAs.I.0, theta = theta.sRNAs.I.0)

model.sRNAs.I <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.1, start = start.sRNAs.1)

## model 2
#df.modelII <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "II")
df.sRNAs.2 <- data.frame(x = df.modelII$qPCR_s, y = df.modelII$CV2)

theta.sRNAs.II.0 <- min(df.sRNAs.2$y) * 0.5  
model.sRNAs.II.0 <- lm(log(y - theta.sRNAs.II.0) ~ x, data=df.sRNAs.2)  
alpha.sRNAs.II.0 <- exp(coef(model.sRNAs.II.0)[1])
beta.sRNAs.II.0 <- coef(model.sRNAs.II.0)[2]

start.sRNAs.2 <- list(alpha = alpha.sRNAs.II.0, beta = beta.sRNAs.II.0, theta = theta.sRNAs.II.0)

model.sRNAs.II <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.2, start = start.sRNAs.2)

## model 3
#df.modelIII <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "III")
df.sRNAs.3 <- data.frame(x = df.modelIII$qPCR_s, y = df.modelIII$CV2)

theta.sRNAs.III.0 <- min(df.sRNAs.3$y) * 0.5  
model.sRNAs.III.0 <- lm(log(y - theta.sRNAs.III.0) ~ x, data=df.sRNAs.3)  
alpha.sRNAs.III.0 <- exp(coef(model.sRNAs.III.0)[1])
beta.sRNAs.III.0 <- coef(model.sRNAs.III.0)[2]

start.sRNAs.3 <- list(alpha = alpha.sRNAs.III.0, beta = beta.sRNAs.III.0, theta = theta.sRNAs.III.0)

model.sRNAs.III <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.3, start = start.sRNAs.3)

## model 4
#df.modelIV <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "IV")
df.sRNAs.4 <- data.frame(x = df.modelIV$qPCR_s, y = df.modelIV$CV2)

theta.sRNAs.IV.0 <- min(df.sRNAs.4$y) * 0.5  
model.sRNAs.IV.0 <- lm(log(y - theta.sRNAs.IV.0) ~ x, data=df.sRNAs.4)  
alpha.sRNAs.IV.0 <- exp(coef(model.sRNAs.IV.0)[1])
beta.sRNAs.IV.0 <- coef(model.sRNAs.IV.0)[2]

start.sRNAs.4 <- list(alpha = alpha.sRNAs.IV.0, beta = beta.sRNAs.IV.0, theta = theta.sRNAs.IV.0)

model.sRNAs.IV <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.4, start = list(beta = beta.sRNAs.IV.0), alg = "plinear")

## model 5
#df.modelV <- Time.qPCR.dPCR.FACS.8clones.814 %>% dplyr::filter(model == "V")
df.sRNAs.5 <- data.frame(x = df.modelV$qPCR_s, y = df.modelV$CV2)

theta.sRNAs.V.0 <- min(df.sRNAs.5$y) * 0.5  
model.sRNAs.V.0 <- lm(log(y - theta.sRNAs.V.0) ~ x, data=df.sRNAs.5)  
alpha.sRNAs.V.0 <- exp(coef(model.sRNAs.V.0)[1])
beta.sRNAs.V.0 <- coef(model.sRNAs.V.0)[2]

start.sRNAs.5 <- list(alpha = alpha.sRNAs.V.0, beta = beta.sRNAs.V.0, theta = theta.sRNAs.V.0)

model.sRNAs.V <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.5, start = list(beta = beta.sRNAs.V.0, theta = theta.sRNAs.V.0), alg = "plinear")

ggplot(df.sRNAs, aes(x = x, y = y, color = model))+geom_point(size = 2.5)+geom_line(aes(x = x, y = predict(model.sRNAs.I, list(x = x))), data = df.sRNAs.1, color = "red")+geom_line(aes(x = x, y = predict(model.sRNAs.II, list(x = x))), data = df.sRNAs.2, color = "orange")+geom_line(aes(x = x, y = predict(model.sRNAs.III, list(x = x))), data = df.sRNAs.3, color = "#1589FF")+geom_line(aes(x = x, y = predict(model.sRNAs.IV, list(x = x))), data = df.sRNAs.4, color = "grey")+geom_line(aes(x = x, y = predict(model.sRNAs.V, list(x = x))), data = df.sRNAs.5, color = "#00CC33")+theme_bw()+ylim(0, 21)+scale_color_manual(values = c("red", "orange", "#1589FF", "grey", "#00CC33"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# asRNAs (qPCR)
## model 1
df.asRNAs.1 <- data.frame(x = df.modelI$qPCR_as, y = df.modelI$CV2)

theta.asRNAs.I.0 <- min(df.asRNAs.1$y) * 0.5  
model.asRNAs.I.0 <- lm(log(y - theta.asRNAs.I.0) ~ x, data=df.asRNAs.1)  
alpha.asRNAs.I.0 <- exp(coef(model.asRNAs.I.0)[1])
beta.asRNAs.I.0 <- coef(model.asRNAs.I.0)[2]

start.asRNAs.1 <- list(alpha = alpha.asRNAs.I.0, beta = beta.asRNAs.I.0, theta = theta.asRNAs.I.0)

model.asRNAs.I <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.1, start = list(beta = beta.asRNAs.I.0, theta = theta.asRNAs.I.0), alg = "plinear")

## model 2
df.asRNAs.2 <- data.frame(x = df.modelII$qPCR_as, y = df.modelII$CV2)

theta.asRNAs.II.0 <- min(df.asRNAs.2$y) * 0.5  
model.asRNAs.II.0 <- lm(log(y - theta.asRNAs.II.0) ~ x, data=df.asRNAs.2)  
alpha.asRNAs.II.0 <- exp(coef(model.asRNAs.II.0)[1])
beta.asRNAs.II.0 <- coef(model.asRNAs.II.0)[2]

start.asRNAs.2 <- list(alpha = alpha.asRNAs.II.0, beta = beta.asRNAs.II.0, theta = theta.asRNAs.II.0)

model.asRNAs.II <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.2, start = start.asRNAs.2)

## model 3
df.asRNAs.3 <- data.frame(x = df.modelIII$qPCR_as, y = df.modelIII$CV2)

theta.asRNAs.III.0 <- min(df.asRNAs.3$y) * 0.5  
model.asRNAs.III.0 <- lm(log(y - theta.asRNAs.III.0) ~ x, data=df.asRNAs.3)  
alpha.asRNAs.III.0 <- exp(coef(model.asRNAs.III.0)[1])
beta.asRNAs.III.0 <- coef(model.asRNAs.III.0)[2]

start.asRNAs.3 <- list(alpha = alpha.asRNAs.III.0, beta = beta.asRNAs.III.0, theta = theta.asRNAs.III.0)

model.asRNAs.III <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.3, start = start.asRNAs.3)

## model 4
df.asRNAs.4 <- data.frame(x = df.modelIV$qPCR_as, y = df.modelIV$CV2)

theta.asRNAs.IV.0 <- min(df.asRNAs.4$y) * 0.5  
model.asRNAs.IV.0 <- lm(log(y - theta.asRNAs.IV.0) ~ x, data=df.asRNAs.4)  
alpha.asRNAs.IV.0 <- exp(coef(model.asRNAs.IV.0)[1])
beta.asRNAs.IV.0 <- coef(model.asRNAs.IV.0)[2]

start.asRNAs.4 <- list(alpha = alpha.asRNAs.IV.0, beta = beta.asRNAs.IV.0, theta = theta.asRNAs.IV.0)

model.asRNAs.IV <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.4, start = list(beta = beta.asRNAs.IV.0), alg = "plinear")

## model 5
df.asRNAs.5 <- data.frame(x = df.modelV$qPCR_as, y = df.modelV$CV2)

theta.asRNAs.V.0 <- min(df.asRNAs.5$y) * 0.5  
model.asRNAs.V.0 <- lm(log(y - theta.asRNAs.V.0) ~ x, data=df.asRNAs.5)  
alpha.asRNAs.V.0 <- exp(coef(model.asRNAs.V.0)[1])
beta.asRNAs.V.0 <- coef(model.asRNAs.V.0)[2]

start.asRNAs.5 <- list(alpha = alpha.asRNAs.V.0, beta = beta.asRNAs.V.0, theta = theta.asRNAs.V.0)

model.asRNAs.V <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.5, start = list(beta = beta.asRNAs.V.0, theta = theta.asRNAs.V.0), alg = "plinear")

ggplot(df.asRNAs, aes(x = x, y = y, color = model))+geom_point(size = 2.5)+geom_line(aes(x = x, y = predict(model.asRNAs.I, list(x = x))), data = df.asRNAs.1, color = "red")+geom_line(aes(x = x, y = predict(model.asRNAs.II, list(x = x))), data = df.asRNAs.2, color = "orange")+geom_line(aes(x = x, y = predict(model.asRNAs.III, list(x = x))), data = df.asRNAs.3, color = "#1589FF")+geom_line(aes(x = x, y = predict(model.asRNAs.IV, list(x = x))), data = df.asRNAs.4, color = "grey")+geom_line(aes(x = x, y = predict(model.asRNAs.V, list(x = x))), data = df.asRNAs.5, color = "#00CC33")+theme_bw()+ylim(0, 21)+scale_color_manual(values = c("red", "orange", "#1589FF", "grey", "#00CC33"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# sRNAs (dPCR)
## model 1
df.sRNAs.1d <- data.frame(x = df.modelI$dPCR_s, y = df.modelI$CV2)

theta.sRNAs.Id.0 <- min(df.sRNAs.1d$y) * 0.5  
model.sRNAs.Id.0 <- lm(log(y - theta.sRNAs.Id.0) ~ x, data=df.sRNAs.1d)  
alpha.sRNAs.Id.0 <- exp(coef(model.sRNAs.Id.0)[1])
beta.sRNAs.Id.0 <- coef(model.sRNAs.Id.0)[2]

start.sRNAs.1d <- list(alpha = alpha.sRNAs.Id.0, beta = beta.sRNAs.Id.0, theta = theta.sRNAs.Id.0)

model.sRNAs.Id <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.1d, start = start.sRNAs.1d)

## model 2
df.sRNAs.2d <- data.frame(x = df.modelII$dPCR_s, y = df.modelII$CV2)

theta.sRNAs.IId.0 <- min(df.sRNAs.2d$y) * 0.5  
model.sRNAs.IId.0 <- lm(log(y - theta.sRNAs.IId.0) ~ x, data=df.sRNAs.2d)  
alpha.sRNAs.IId.0 <- exp(coef(model.sRNAs.IId.0)[1])
beta.sRNAs.IId.0 <- coef(model.sRNAs.IId.0)[2]

start.sRNAs.2d <- list(alpha = alpha.sRNAs.IId.0, beta = beta.sRNAs.IId.0, theta = theta.sRNAs.IId.0)

model.sRNAs.IId <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.2d, start = start.sRNAs.2d)

## model 3
df.sRNAs.3d <- data.frame(x = df.modelIII$dPCR_s, y = df.modelIII$CV2)

theta.sRNAs.IIId.0 <- min(df.sRNAs.3d$y) * 0.5  
model.sRNAs.IIId.0 <- lm(log(y - theta.sRNAs.IIId.0) ~ x, data=df.sRNAs.3d)  
alpha.sRNAs.IIId.0 <- exp(coef(model.sRNAs.IIId.0)[1])
beta.sRNAs.IIId.0 <- coef(model.sRNAs.IIId.0)[2]

start.sRNAs.3d <- list(alpha = alpha.sRNAs.IIId.0, beta = beta.sRNAs.IIId.0, theta = theta.sRNAs.IIId.0)

model.sRNAs.IIId <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.3d, start = list(beta = beta.sRNAs.IIId.0), alg = "plinear")

## model 4
df.sRNAs.4d <- data.frame(x = df.modelIV$dPCR_s, y = df.modelIV$CV2)

theta.sRNAs.IVd.0 <- min(df.sRNAs.4d$y) * 0.5  
model.sRNAs.IVd.0 <- lm(log(y - theta.sRNAs.IVd.0) ~ x, data=df.sRNAs.4d)  
alpha.sRNAs.IVd.0 <- exp(coef(model.sRNAs.IVd.0)[1])
beta.sRNAs.IVd.0 <- coef(model.sRNAs.IVd.0)[2]

start.sRNAs.4d <- list(alpha = alpha.sRNAs.IVd.0, beta = beta.sRNAs.IVd.0, theta = theta.sRNAs.IVd.0)

model.sRNAs.IVd <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.4d, start = list(beta = beta.sRNAs.IVd.0), alg = "plinear")

## model 5
df.sRNAs.5d <- data.frame(x = df.modelV$dPCR_s, y = df.modelV$CV2)

theta.sRNAs.Vd.0 <- min(df.sRNAs.5d$y) * 0.5  
model.sRNAs.Vd.0 <- lm(log(y - theta.sRNAs.Vd.0) ~ x, data=df.sRNAs.5d)  
alpha.sRNAs.Vd.0 <- exp(coef(model.sRNAs.Vd.0)[1])
beta.sRNAs.Vd.0 <- coef(model.sRNAs.Vd.0)[2]

start.sRNAs.5d <- list(alpha = alpha.sRNAs.Vd.0, beta = beta.sRNAs.Vd.0, theta = theta.sRNAs.Vd.0)

model.sRNAs.Vd <- nls(y ~ alpha * exp(beta * x) + theta , data = df.sRNAs.5d, start = list(beta = beta.sRNAs.Vd.0), alg = "plinear")

ggplot(df.sRNAs.d, aes(x = x, y = y, color = model))+geom_point(size = 2.5)+geom_line(aes(x = x, y = predict(model.sRNAs.Id, list(x = x))), data = df.sRNAs.1d, color = "red")+geom_line(aes(x = x, y = predict(model.sRNAs.IId, list(x = x))), data = df.sRNAs.2d, color = "orange")+geom_line(aes(x = x, y = predict(model.sRNAs.IIId, list(x = x))), data = df.sRNAs.3d, color = "#1589FF")+geom_line(aes(x = x, y = predict(model.sRNAs.IVd, list(x = x))), data = df.sRNAs.4d, color = "grey")+geom_line(aes(x = x, y = predict(model.sRNAs.Vd, list(x = x))), data = df.sRNAs.5d, color = "#00CC33")+theme_bw()+ylim(0, 21)+scale_color_manual(values = c("red", "orange", "#1589FF", "grey", "#00CC33"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))

# asRNAs (dPCR)
## model 1
df.asRNAs.1d <- data.frame(x = df.modelI$dPCR_as, y = df.modelI$CV2)

theta.asRNAs.Id.0 <- min(df.asRNAs.1d$y) * 0.5  
model.asRNAs.Id.0 <- lm(log(y - theta.asRNAs.Id.0) ~ x, data=df.asRNAs.1d)  
alpha.asRNAs.Id.0 <- exp(coef(model.asRNAs.Id.0)[1])
beta.asRNAs.Id.0 <- coef(model.asRNAs.Id.0)[2]

start.asRNAs.1d <- list(alpha = alpha.asRNAs.Id.0, beta = beta.asRNAs.Id.0, theta = theta.asRNAs.Id.0)

model.asRNAs.Id <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.1d, start = list(beta = beta.asRNAs.Id.0), alg = "plinear")

## model 2
df.asRNAs.2d <- data.frame(x = df.modelII$dPCR_as, y = df.modelII$CV2)

theta.asRNAs.IId.0 <- min(df.asRNAs.2d$y) * 0.5  
model.asRNAs.IId.0 <- lm(log(y - theta.asRNAs.IId.0) ~ x, data=df.asRNAs.2d)  
alpha.asRNAs.IId.0 <- exp(coef(model.asRNAs.IId.0)[1])
beta.asRNAs.IId.0 <- coef(model.asRNAs.IId.0)[2]

start.asRNAs.2d <- list(alpha = alpha.asRNAs.IId.0, beta = beta.asRNAs.IId.0, theta = theta.asRNAs.IId.0)

model.asRNAs.IId <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.2d, start = start.asRNAs.2d)

## model 3
df.asRNAs.3d <- data.frame(x = df.modelIII$dPCR_as, y = df.modelIII$CV2)

theta.asRNAs.IIId.0 <- min(df.asRNAs.3d$y) * 0.5  
model.asRNAs.IIId.0 <- lm(log(y - theta.asRNAs.IIId.0) ~ x, data=df.asRNAs.3d)  
alpha.asRNAs.IIId.0 <- exp(coef(model.asRNAs.IIId.0)[1])
beta.asRNAs.IIId.0 <- coef(model.asRNAs.IIId.0)[2]

start.asRNAs.3d <- list(alpha = alpha.asRNAs.IIId.0, beta = beta.asRNAs.IIId.0, theta = theta.asRNAs.IIId.0)

model.asRNAs.IIId <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.3d, start = list(beta = beta.asRNAs.IIId.0), alg = "plinear")

## model 4
df.asRNAs.4d <- data.frame(x = df.modelIV$dPCR_as, y = df.modelIV$CV2)

theta.asRNAs.IVd.0 <- min(df.asRNAs.4d$y) * 0.5  
model.asRNAs.IVd.0 <- lm(log(y - theta.asRNAs.IVd.0) ~ x, data=df.asRNAs.4d)  
alpha.asRNAs.IVd.0 <- exp(coef(model.asRNAs.IVd.0)[1])
beta.asRNAs.IVd.0 <- coef(model.asRNAs.IVd.0)[2]

start.asRNAs.4d <- list(alpha = alpha.asRNAs.IVd.0, beta = beta.asRNAs.IVd.0, theta = theta.asRNAs.IVd.0)

model.asRNAs.IVd <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.4d, start = list(beta = beta.asRNAs.IVd.0, theta = theta.asRNAs.IVd.0), alg = "plinear")

## model 5
df.asRNAs.5d <- data.frame(x = df.modelV$dPCR_as, y = df.modelV$CV2)

theta.asRNAs.Vd.0 <- min(df.asRNAs.5d$y) * 0.5  
model.asRNAs.Vd.0 <- lm(log(y - theta.asRNAs.Vd.0) ~ x, data=df.asRNAs.5d)  
alpha.asRNAs.Vd.0 <- exp(coef(model.asRNAs.Vd.0)[1])
beta.asRNAs.Vd.0 <- coef(model.asRNAs.Vd.0)[2]

start.asRNAs.5d <- list(alpha = alpha.asRNAs.Vd.0, beta = beta.asRNAs.Vd.0, theta = theta.asRNAs.Vd.0)

model.asRNAs.Vd <- nls(y ~ alpha * exp(beta * x) + theta , data = df.asRNAs.5d, start = list(beta = beta.asRNAs.Vd.0, theta = theta.asRNAs.Vd.0), alg = "plinear")

ggplot(df.asRNAs.d, aes(x = x, y = y, color = model))+geom_point(size = 2.5)+geom_line(aes(x = x, y = predict(model.asRNAs.Id, list(x = x))), data = df.asRNAs.1d, color = "red")+geom_line(aes(x = x, y = predict(model.asRNAs.IId, list(x = x))), data = df.asRNAs.2d, color = "orange")+geom_line(aes(x = x, y = predict(model.asRNAs.IIId, list(x = x))), data = df.asRNAs.3d, color = "#1589FF")+geom_line(aes(x = x, y = predict(model.asRNAs.IVd, list(x = x))), data = df.asRNAs.4d, color = "grey")+geom_line(aes(x = x, y = predict(model.asRNAs.Vd, list(x = x))), data = df.asRNAs.5d, color = "#00CC33")+theme_bw()+ylim(0, 21)+scale_color_manual(values = c("red", "orange", "#1589FF", "grey", "#00CC33"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 45,vjust = 1, hjust = 1), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
