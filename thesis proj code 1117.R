# Libraries and data loading----
rm(list=ls())

#install.packages(data.table)
#install.packages(calibrate)
library(data.table)
library(calibrate)
library(ggplot2)
library(MASS)  # Implements LDA & QDA
library(caret) # confusionMatrix
library(pROC)  # ROC & AUC
library(class) # Implements KNN
library(robCompositions)
library(ToolsForCoDa)

ECU.MF.raw <- read.csv("D:/24 Winter UW/PCA/ECU-MF.csv")

colnames(ECU.MF.raw) <-  c("sex","age","leukocytes","neutrophilsP","lymphocytesP",
                           "monocytesP","eosinophilsP","basophilsP","neutrophils",
                           "lymphocytes","NLR","monocytes","eosinophils",
                           "basophils","redbloodcells","mcv","mch","mchc",
                           "rdwP","hemoglobin","hematocritP","platelets","mpv","Condition")

dim(ECU.MF.raw)
sapply(ECU.MF.raw, is.numeric)
ECU.MF.raw$sex <- as.factor(ECU.MF.raw$sex)

sum(is.na(ECU.MF.raw))

# Summary box plots----
# Positive_df <- ECU.MF.raw[ECU.MF.raw$Condition == 1,]
# Negative_df <- ECU.MF.raw[ECU.MF.raw$Condition == 0,]
# 
# par(mfrow = c(3,7))
# 
# for (i in 3:23) {
#   boxplot(ECU.MF.raw[,i] ~ ECU.MF.raw[,24],
#           xlab = colnames(ECU.MF.raw)[i], ylab = NA)
# }
# ----

# Data processing----
ECU.MF.counts.raw <- ECU.MF.raw[, c(9,10,12:14)]

diff_wbc <- ECU.MF.raw$leukocytes - rowSums(ECU.MF.counts.raw)
rdiff_wbc <- abs(diff_wbc)/rowSums(ECU.MF.counts.raw)

exclude_index <- which(rdiff_wbc > 0.1)
length(exclude_index)

  # difference behavior before filtering
  par(mfrow = c(1,2))
  boxplot(diff_wbc, xlab = "Leukocytes counts - (sum of 5 WBC)")
  boxplot(rdiff_wbc, xlab = "|Leukocytes counts/(sum of 5 WBC) - 1| ")
  par(mfrow = c(1, 1))
  par(pty = "s")
  plot(x = ECU.MF.raw$leukocytes, y = rowSums(ECU.MF.counts.raw), 
       xlab = "leukocytes counts", ylab = "wbc counts", 
       main = "leukocytes counts vs (sum of 5 WBC) raw")
  lines(x = ECU.MF.raw$leukocytes, y = ECU.MF.raw$leukocytes, col="red", lty=2)
  legend(800, 90000, legend=c("line 1: y=x"), col=c("red"), lty=1:2, cex=0.8)

  # Filtering
  ECU.MF.counts <- ECU.MF.counts.raw[-exclude_index, ]
  ECU.MF <- ECU.MF.raw[-exclude_index, ]
  diff_wbc_filtered <- ECU.MF$leukocytes - rowSums(ECU.MF.counts)
  rdiff_wbc_filtered <- abs(diff_wbc_filtered)/rowSums(ECU.MF.counts)
  
  # difference behavior after filtering
  par(mfrow = c(1,2))
  boxplot(diff_wbc_filtered, xlab = "Leukocytes counts - (sum of 5 WBC) after filtering")
  boxplot(rdiff_wbc_filtered, xlab = "|Leukocytes counts/(sum of 5 WBC) - 1| after filtering")
  par(mfrow = c(1, 1))
  par(pty = "s")
  plot(x = ECU.MF$leukocytes, y = rowSums(ECU.MF.counts), 
       xlab = "leukocytes counts", ylab = "wbc counts", 
       main = "leukocytes counts vs (sum of 5 WBC) filtered")
  lines(x = ECU.MF$leukocytes, y = ECU.MF$leukocytes, col="red", lty=2)
  legend(2000, 13000, legend=c("line 1: y=x"), col=c("red"), lty=1:2, cex=0.8)
  
  # combine boxplots
  par(mfrow = c(2,2))
  boxplot(diff_wbc, xlab = "Leukocytes counts - (sum of 5 WBC)")
  title("Before filtering")
  boxplot(rdiff_wbc, xlab = "|Leukocytes counts/(sum of 5 WBC) - 1| ")
  title("Before filtering")
  boxplot(diff_wbc_filtered, xlab = "Leukocytes counts - (sum of 5 WBC)")
  title("After filtering")
  boxplot(rdiff_wbc_filtered, xlab = "|Leukocytes counts/(sum of 5 WBC) - 1|")
  title("After filtering")
  par(mfrow =c(1,1))


colvec <- rep(NA,nrow(ECU.MF))
colvec[ECU.MF$Condition==0] <- "green"
colvec[ECU.MF$Condition==1] <- "red"

condition <- factor(ECU.MF$Condition)

ECU.MF <- ECU.MF[, c(-1,-3,-4:-8,-11,-24)]

n <- nrow(ECU.MF)
p <- ncol(ECU.MF)

# PCA with all variables (use corr)----

    out.pca <- princomp(ECU.MF, cor=T)
    
    plot(out.pca, type = "lines", main = "Scree plot")
    
    biplot(out.pca)
    
    
    # choose consistent definition (use out.pca or manually)
    
    Fp <- out.pca$scores
    Gs <- out.pca$loadings
    
    # or manually do PCA, remember to be consistent with our use as PCA e-vectors doesn't care about directions (signals)
    # S <- cor(ECU.MF)
    # V <- eigen(S)$vectors
    # Dl <- diag(eigen(S)$values)
    # 
    # Xc <- scale(ECU.MF,scale=T) # use corr instead of cov
    # 
    # Fp <- Xc%*%V
    # Gs <- V
    
    la <- apply(Fp,2,var)
    fr <- la/sum(la)
    cu <- cumsum(fr)
    
    # bplot(Fp,Gs,cex.rowlab = 0.25,rowch=1)
    # arrows(0,0,out.pca$loadings[,1],out.pca$loadings[,2], length = 2)
    
    bplot(Fp,8*Gs,cex.rowlab = 0.25,colch=NA,collab=colnames(ECU.MF),rowch=1, rowcol=colvec,
          cex.collab=0.5,main="Form biplot")

    # remove outliers
    plot(Fp[,1],Fp[,2],asp=1, col=colvec)
    ol <- which(Fp[,1]> 8)

    plot(Fp[,1][-ol],Fp[,2][-ol],asp=1,col=colvec[-ol],
         xlab="First principal axis (20.21%)",
         ylab="Second principal axis (15.25%)",main="CBC")
    abline(h=0,lty="dotted")
    abline(v=0,lty="dotted")
    

    bplot(Fp[-ol,],8*Gs,cex.rowlab = 0.25,colch=NA,collab=colnames(ECU.MF),rowch=1, rowcol=colvec[-ol],
          cex.collab=0.8,main="Form biplot",
          xlab = "First principal axis (20.21%)", ylab = "Second principal axis (15.25%)",
          cex.axis = 1)
    legend(-5, 6, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.8)

# PCA with 5 counts variables----
    
    # colnames(ECU.MF.counts) <- c("n","l","m","e","b")
    
    # pca with 5 counts variables (use corr)
    
    out.pca.counts <- princomp(ECU.MF.counts,cor=F)
    
    plot(out.pca.counts, type = "lines", main = "Scree plot")
    
    biplot(out.pca.counts)
    
    # choose consistent definition (use out.pca or manually)
    
    Fp.counts <- out.pca.counts$scores
    Gs.counts <- out.pca.counts$loadings
    
    la.counts <- apply(Fp.counts,2,var)
    fr.counts <- la.counts/sum(la.counts)
    cu.counts <- cumsum(fr.counts)
    
    # or manually do PCA, remember to be consistent with our use as PCA e-vectors doesn't care about directions (signals)
    
    # bplot(Fp,Gs,cex.rowlab = 0.25,rowch=1)
    # arrows(0,0,out.pca$loadings[,1],out.pca$loadings[,2], length = 2)
    
    bplot(Fp.counts,3000*Gs.counts,cex.rowlab = 0.25,colch=NA,collab=colnames(ECU.MF.counts),rowch=1, rowcol=colvec,
          cex.collab=0.8,main="Form biplot",
          xlab = "First principal axis (76.47%)", ylab = "Second principal axis (21.40%)",
          cex.axis = 1)
    legend(-3500, 4000, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.8)
    
    
    # remove outliers, we dont really have any after filtering
    # plot(Fp.counts[,1],Fp.counts[,2],asp=1, col=colvec)
    # ol.counts <- which(Fp.counts[,1]>20000)
    # 
    # plot(Fp.counts[,1][-ol.counts],Fp.counts[,2][-ol.counts],asp=1,col=colvec[-ol.counts],
    #      xlab="First principal axis",
    #      ylab="Second principal axis",main="CBC (5 major white blood cell counts)")
    # abline(h=0,lty="dotted")
    # abline(v=0,lty="dotted")
    # 
    # bplot(Fp.counts[-ol.counts,],5000*Gs.counts,cex.rowlab = 0.25,colch=NA,collab=colnames(ECU.MF.counts),rowch=1, rowcol=colvec[-ol.counts],
    #       cex.collab=1,main="Form biplot")
  
    
    # cov biplot:
    S <- cov(ECU.MF.counts)
    V <- eigen(S)$vectors
    Dl <- diag(eigen(S)$values)
    
    Xc <- scale(ECU.MF.counts, scale=FALSE)
    
    Fp <- Xc%*%V
    Gs <- V
    
    Fs <- Fp%*%sqrt(solve(Dl))
    Gp <- Gs%*%sqrt(Dl)
    
    plot(Fs[,1],Fs[,2],asp=1, col=colvec)
    
    bplot(Fs,0.0025*Gp,cex.rowlab=0.25,colch=NA,collab=colnames(ECU.MF.counts),rowch=1, rowcol=colvec,
          cex.collab= 0.5,main="Covariance biplot",
          xlab = "First principal axis", ylab = "Second principal axis",
          cex.axis = 1)
    legend(-2.5, 3.5, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.6)
    
    # PCA with WBC (using correlation-based SVD)
    # form biplot:
    out.pca.counts <- princomp(ECU.MF.counts,cor=T)
    
    plot(out.pca.counts, type = "lines", main = "Scree plot")
    
    biplot(out.pca.counts)
    
    Fp.counts <- out.pca.counts$scores
    Gs.counts <- out.pca.counts$loadings
    
    la.counts <- apply(Fp.counts,2,var)
    fr.counts <- la.counts/sum(la.counts)
    cu.counts <- cumsum(fr.counts)
    
    bplot(Fp.counts,5*Gs.counts,cex.rowlab = 0.25,colch=NA,collab=colnames(ECU.MF.counts),rowch=1, rowcol=colvec,
          cex.collab=0.8,main="Form biplot",
          xlab = "First principal axis (43.83%)", ylab = "Second principal axis (19.96%)",
          cex.axis = 1)
    legend(-3.5, 4, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.8)
    
    
# GLM and LDA----
    # set train dataset and test dataset
    set.seed(0)
    
    df_raw <- ECU.MF
    df_scale <- data.frame(scale(ECU.MF))
    
    df_raw$Condition <- factor(ECU.MF.raw[-exclude_index, 24])
    df_scale$Condition <- factor(ECU.MF.raw[-exclude_index, 24])
    
    df_full <- df_scale
    df_WBC <- df_raw[,c(2:6, 16)]

    #############################################
    ################# LDA #######################
    #############################################
    
    # Full data----
    lda.model.full = lda(Condition ~ ., data = df_full, center = TRUE)
    lda.model.full.predict = predict(lda.model.full, newdata = df_full)
    tt.lda.full = table(Prediction = lda.model.full.predict$class, True = df_full$Condition)
    tt.lda.full
    mean(lda.model.full.predict$class == df_full$Condition)
    
    x0 <- apply(df_full[df_full$Condition == "0", -16], 2, mean)
    x1 <- apply(df_full[df_full$Condition == "1", -16], 2, mean)
    
    X_neg <- as.matrix(df_full[df_full$Condition == "0", -16])
    n1 <- nrow(X_neg)
    X_pos <- as.matrix(df_full[df_full$Condition == "1", -16])
    n2 <- nrow(X_pos)
    
    Sp <- (n1 - 1)/(n1 + n2 - 2) * cov(X_neg) + (n2 - 1)/(n1 + n2 - 2) * cov(X_pos)
    scaling <- solve(as.matrix(Sp)) %*% (x1 - x0)
    
    classifier <- (1/2)*sum(x0*scaling + x1*scaling) + 
                  log(lda.model.full$prior[1] / lda.model.full$prior[2])
    
    X = as.matrix(df_full[,-16]) %*% scaling
    
    sum(X[df_full$Condition=="0"] > classifier)
    sum(X[df_full$Condition=="1"] < classifier)
    
    plot(x = X, y = jitter(as.numeric(df_full$Condition, amount=0.1) - 1), 
         ylim = c(-0.5, 1.5), 
         yaxp = c(0, 1, 1),
         col = colvec,
         xlab = "Linear discriminant", ylab = "Group",
         main = "LDA using full dataset")
    abline(v = classifier,lty="dotted")
    abline(h = c(0,1), lty="dotted")
    points(x = mean(X[df_full$Condition=="0"]), y=0, pch = 19, 
           col = "darkgreen", cex = 1.5)
    points(x = mean(X[df_full$Condition=="1"]), y=1, pch = 19, 
           col = "darkred", cex = 1.5)
    legend(x = -6.5, y = 1.5, legend=c("PCR Positive", "PCR Negative", 
                                       "PCR Positive Mean", "PCR Negative Mean"),  
           fill = c("red","green", "darkred", "darkgreen"), cex = 0.6)
    legend(x = classifier+0.2, y = -0.3, legend=c("Threshold = 1.077882"),  
           col = c("grey"), lty = c(2), cex = 0.6)
    
    
    coef_full <- data.frame(x = rownames(scaling), y = scaling)
    ggplot(coef_full, aes(x = x, y=y)) + 
      geom_bar(stat = "identity", width = 0.5) +
      coord_flip() +
      labs( y = "LDA coefficient", x = "")
    
    roc_score_full = roc(response = df_full$Condition, predictor = lda.model.full.predict$posterior[, 2])
    roc_score_full
    roc_full = ggroc(roc_score_full, linetype=1, size = 1) + 
      xlab("FPR") + ylab("TPR") + 
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_full
    
    
    lda.model.full.cv = lda(Condition ~ ., data = df_full, center = TRUE, CV = TRUE)
    tt.lda = table(Prediction = lda.model.full.cv$class, True = df_full$Condition)
    tt.lda
    mean(lda.model.full.cv$class == df_full$Condition)
    
    roc_score_full_cv = roc(response = df_full$Condition, predictor = lda.model.full.cv$posterior[, 2])
    roc_score_full_cv
    
    
    # WBC data----
    lda.model.WBC = lda(Condition ~ ., data = df_WBC, center = TRUE)
    lda.model.WBC.predict = predict(lda.model.WBC, newdata = df_WBC)
    tt.lda.WBC = table(Prediction = lda.model.WBC.predict$class, True = df_WBC$Condition)
    tt.lda.WBC
    mean(lda.model.WBC.predict$class == df_WBC$Condition)
    
    x0 <- apply(df_WBC[df_WBC$Condition == "0", -6], 2, mean)
    x1 <- apply(df_WBC[df_WBC$Condition == "1", -6], 2, mean)
    
    X_neg <- as.matrix(df_WBC[df_WBC$Condition == "0", -6])
    n1 <- nrow(X_neg)
    X_pos <- as.matrix(df_WBC[df_WBC$Condition == "1", -6])
    n2 <- nrow(X_pos)
    
    Sp <- (n1 - 1)/(n1 + n2 - 2) * cov(X_neg) + (n2 - 1)/(n1 + n2 - 2) * cov(X_pos)
    scaling <- solve(as.matrix(Sp)) %*% (x1 - x0)
    
    classifier <- (1/2)*sum(x0*scaling + x1*scaling) + 
      log(lda.model.WBC$prior[1] / lda.model.WBC$prior[2])
    
    X = as.matrix(df_WBC[, -6]) %*% scaling
    
    sum(X[df_WBC$Condition=="0"] > classifier)
    sum(X[df_WBC$Condition=="1"] < classifier)
    
    plot(x = X, y = jitter(as.numeric(df_full$Condition, amount=0.1) - 1), 
         ylim = c(-0.5, 1.5), 
         yaxp = c(0, 1, 1),
         col = colvec,
         xlab = "Linear discriminant", ylab = "Group",
         main = "LDA using WBC dataset")
    abline(v = classifier,lty="dotted")
    abline(h = c(0,1), lty="dotted")
    points(x = mean(X[df_full$Condition=="0"]), y=0, pch = 19, 
           col = "darkgreen", cex = 1.5)
    points(x = mean(X[df_full$Condition=="1"]), y=1, pch = 19, 
           col = "darkred", cex = 1.5)
    legend(x = -9.5, y = 1.5, legend=c("PCR Positive", "PCR Negative", 
                                       "PCR Positive Mean", "PCR Negative Mean"),  
           fill = c("red","green", "darkred", "darkgreen"), cex = 0.6)
    legend(x = classifier - 3, y = -0.3, legend=c("Threshold = -2.963664 "),  
           col = c("grey"), lty = c(2), cex = 0.6)
    
    coef_full <- data.frame(x = rownames(scaling), y = scaling)
    ggplot(coef_full, aes(x = x, y=y)) + 
      geom_bar(stat = "identity", width = 0.3) +
      coord_flip() +
      labs( y = "LDA coefficient", x = "")
    
    roc_score_WBC = roc(response = df_WBC$Condition, predictor = lda.model.WBC.predict$posterior[, 2])
    roc_score_WBC
    roc_WBC = ggroc(roc_score_WBC, linetype=1, size = 1) + 
      xlab("FPR") + ylab("TPR") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_WBC
    
    roc_combined = ggroc(list("Full" = roc_score_full, "WBC" = roc_score_WBC), 
                         linetype=1, size = 1) + 
      xlab("FPR") + ylab("TPR") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_combined
    
    posterior_prefect <- rep(0, length(condition))
    posterior_prefect[condition == "1"] <- 1
    posterior_random <- rep(0.5, length(condition))
    roc_score_eg_perfect = roc(response = condition, predictor = posterior_prefect)
    roc_score_eg_random = roc(response = condition, predictor = posterior_random)
    roc_eg <- ggroc(list("Example ROC curve" = roc_score_WBC, 
                         "Ideal ROC curve" = roc_score_eg_perfect,
                         "Random classification" = roc_score_eg_random), 
                    linetype=1, size = 0.5) + 
              xlab("FPR") + ylab("TPR")
    roc_eg
    
    
    lda.model.WBC.cv = lda(Condition ~ ., data = df_WBC, center = TRUE, CV = TRUE)
    tt.lda = table(Prediction = lda.model.WBC.cv$class, True = df_WBC$Condition)
    tt.lda
    mean(lda.model.WBC.cv$class == df_WBC$Condition)
    
    roc_score_WBC_cv = roc(response = df_WBC$Condition, predictor = lda.model.WBC.cv$posterior[, 2])
    roc_score_WBC_cv
    

    
    

    
    #############################################
    ################# GLM #######################
    #############################################
    
    # glm.model = glm(formula = Condition ~ ., family = binomial(link = "logit"), data = train_df)
    # summary(glm.model)
    # 
    # contrasts(train_df$Condition)
    # coef(glm.model)
    # 
    # confint(glm.model)
    # 
    # glm.prob.train = predict(glm.model, type = "response")
    # 
    # glm.label.train = ifelse(glm.prob.train > .5, "1", "0")
    # glm.label.train = factor(glm.label.train, levels = c("0","1"))
    # 
    # # Misclassification error rate
    # 1-mean(glm.label.train == train_df$Condition)
    # 
    # # Confusion matrix
    # tt.glm.train = table(True = train_df$Condition, Predicted = glm.label.train)
    # tt.glm.train
    # 
    # # Alternative approach
    # confusionMatrix(glm.label.train, train_df$Condition,
    #                 positive = "1")
    # 
    # # Test set performances
    # glm.prob.test = predict(glm.model,type = "response", newdata = test_df)
    # 
    # glm.label.test = ifelse(glm.prob.test > .5, "1", "0")
    # glm.label.test = factor(glm.label.test, levels = c("0","1"))
    # 
    # tt.glm.test = table(True = test_df$Condition, Predicted = glm.label.test)
    # tt.glm.test
    # 1-mean(glm.label.test == test_df$Condition)
    # 
    # # ROC curve (test set)
    # 
    # roc_score = roc(response = test_df$Condition, predictor = glm.prob.test)  # AUC score
    # roc_score # get AUC
    # 
    # roc_1 = ggroc(roc_score, linetype=1, size = 1) + 
    #   xlab("FPR") + ylab("TPR") +
    #   ggtitle("Test ROC curve")
    

    

    
# LR transformation ----
    # transform to compositional data----
    Co_df <- ECU.MF.counts
    zero_baso_id <- which(Co_df$basophils == 0)
      # remove zeros----
      library(zCompositions)
      Co_df <- cmultRepl(Co_df)

    Co_df$sum <- rowSums(Co_df)
    for (i in 1:5) {
      Co_df[,i] <- Co_df[,i] / Co_df[,6]
    }
    
    # # try remove 8 outliers----
    # Co_df <- ECU.MF.counts
    # zero_baso_id <- which(Co_df$basophils == 0)
    # Co_df <- Co_df[-zero_baso_id,]
    # Co_df$sum <- rowSums(Co_df)
    # for (i in 1:5) {
    #   Co_df[,i] <- Co_df[,i] / Co_df[,6]
    # }
    

    # clr----
    clr_df <- Co_df[,1:5]
    clr_df$geo_mean <- apply(clr_df, 1, function(x){exp(mean(log(x)))})
    for (i in 1:5) {
      clr_df[,i] <- log(clr_df[,i] / clr_df[,6])
    }
    clr_df <- clr_df[,1:5]
    
    
    # apply PCA----
    out.pca.counts <- princomp(clr_df,cor=F)
    
    plot(out.pca.counts, type = "lines", main = "Scree plot")
    
    biplot(out.pca.counts)
    
    # choose consistent definition (use out.pca or manually)
    
    Fp.counts <- out.pca.counts$scores
    Gs.counts <- out.pca.counts$loadings
    
    la.counts <- apply(Fp.counts,2,var)
    fr.counts <- la.counts/sum(la.counts)
    cu.counts <- cumsum(fr.counts)
    
    # or manually do PCA, remember to be consistent with our use as PCA e-vectors doesn't care about directions (signals)
    
    # bplot(Fp,Gs,cex.rowlab = 0.25,rowch=1)
    # arrows(0,0,out.pca$loadings[,1],out.pca$loadings[,2], length = 2)
    
    bplot(Fp.counts,1*Gs.counts,cex.rowlab = 0.25,colch=NA,collab=colnames(clr_df),rowch=1, rowcol=colvec,
          cex.collab=0.8,main="Form biplot",
          xlab = "First principal axis (46.42%)", ylab = "Second principal axis (37.05%)",
          cex.axis = 1)
    legend(-2, 2, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.6)
    
    
    # cov biplot:
    S <- cov(clr_df)
    V <- eigen(S)$vectors
    Dl <- diag(eigen(S)$values)
    Dl[5,5] <- 2e-18
    
    Xc <- scale(clr_df, scale=FALSE)
    
    Fp <- Xc%*%V
    Gs <- V
    
    Fs <- Fp%*%sqrt(solve(Dl, tol = 1e-20))
    Gp <- Gs%*%sqrt(Dl)
    
    plot(Fs[,1],Fs[,2],asp=1, col=colvec)
    
    bplot(Fs,1.3*Gp,cex.rowlab=0.25,colch=NA,collab=colnames(clr_df),rowch=1, rowcol=colvec,
          cex.collab= 0.8,main="Covariance biplot",
          xlab = "First principal axis", ylab = "Second principal axis",
          cex.axis = 1)
    legend(-3, 3, legend=c("PCR Positive", "PCR Negative"),  
           fill = c("red","green"), cex = 0.6)

    
# LR-LDA ----
    library(mvtnorm)
    
    # library(ToolsForCoDa)
    # source("D:/24 Winter UW/PCA/LR-LDA code/lrlda.R")
    # source("D:/24 Winter UW/PCA/LR-LDA code/tab2vec.R")
    
    # source("D:/24 Winter UW/PCA/LR-LDA code/ToolsForCoDa/R/lrlda.R")
    # source("D:/24 Winter UW/PCA/LR-LDA code/ToolsForCoDa/R/tab2vec.R")
    
    lrlda_df <- Co_df[,1:5]
    out.lrlda <- lrlda(lrlda_df, condition)
    
    out.lrlda$gsizes
    out.lrlda$Mclr
    head(out.lrlda$LD)
    out.lrlda$CM
    head(round(out.lrlda$prob.posterior,4))
    
    # Fp <- out.lrlda$Fp
    # Fc <- out.lrlda$Fc
    LD <- out.lrlda$LD
    
    # library(Correlplot) ----
    # lims <- jointlim(Fc,Fp)
    # opar <- par(bty="n",xaxt="n",yaxt="n")
    # plot(Fp[,1],Fp[,2],asp=1,pch=17,xlab="",ylab="",col=c("red","green","blue"),
    #      xlim=lims$xlim,ylim=lims$ylim,cex=1.25)
    # points(Fc[,1],Fc[,2],col=colvec)
    # origin()
    # arrows(0,0,10*Gs[,1],10*Gs[,2],angle = 10, length = 0.1)
    # textxy(10*Gs[,1],10*Gs[,2],colnames(Oxides))
    # par(opar)
    # legend("topleft",c("G","NF","W"),col=c("red","green","blue"),pch=1,cex=0.5)
    # 
    # lrMgOAl2O2 <- Oxides[,c("MgO")]/Oxides[,c("Al2O3")]
    # boxplot(lrMgOAl2O2~site,col=c("red","green","blue"))
    # ----
    
    plot(x = LD[-zero_baso_id], y = jitter(as.numeric(df_full$Condition[-zero_baso_id], amount=0.1) - 1), 
         ylim = c(-0.5, 3), 
         yaxp = c(0, 1, 1),
         col = colvec[- zero_baso_id],
         xlab = "LD", ylab = "Group",
         main = "LR-LDA")
    abline(h = c(0,1), lty="dotted")
    abline(v = 0, lty = "dotted")
    points(x = LD[zero_baso_id], y = as.numeric(condition[zero_baso_id]) - 1, pch = 2)
    points(x = mean(LD[df_full$Condition=="0"]), y=0, pch = 19, 
           col = "darkgreen", cex = 1.5)
    points(x = mean(LD[df_full$Condition=="1"]), y=1, pch = 19, 
           col = "darkred", cex = 1.5)
    
    arrows(0, 2.4, 3*out.lrlda$Gs[1], length = 0.1, col = "lightskyblue4")
    text(3*out.lrlda$Gs[1], 2.4, labels = "neutrophils", pos = 4, col = "lightskyblue4", cex = 0.7)
    
    arrows(0, 2.2, 3*out.lrlda$Gs[2], length = 0.1, col = "lightskyblue3")
    text(3*out.lrlda$Gs[2], 2.2, labels = "lymphocytes", pos = 4, col = "lightskyblue3", cex = 0.7)
    
    arrows(0, 2.0, 3*out.lrlda$Gs[3], length = 0.1, col = "lightskyblue2")
    text(3*out.lrlda$Gs[3], 2.0, labels = "monocytes", pos = 4, col = "lightskyblue2", cex = 0.7)
    
    arrows(0, 1.8, 3*out.lrlda$Gs[4], length = 0.1, col = "lightskyblue1")
    text(3*out.lrlda$Gs[4], 1.8, labels = "eosinophils", pos = 2, col = "lightskyblue1", cex = 0.7)
    
    arrows(0, 1.6, 3*out.lrlda$Gs[5], length = 0.1, col = "lightskyblue")
    text(3*out.lrlda$Gs[5], 1.6, labels = "basophils", pos = 2, col = "lightskyblue", cex = 0.7)
    
    legend(x = -3, y = 3.1, legend=c("PCR Positive", "PCR Negative", 
                                       "PCR Positive Mean", "PCR Negative Mean",
                                     "Zero basophils cases"),
           pch = c(1, 1, 19, 19, 2),
           col = c("red","green", "darkred", "darkgreen","black"), cex = 0.6)
    
    roc_score_lrlda = roc(response = condition, predictor = out.lrlda$prob.posterior[,2])
    roc_score_lrlda
    roc_lrlda = ggroc(roc_score_lrlda, linetype=1, size = 1) + 
      xlab("FPR") + ylab("TPR") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_lrlda
    

    
    # try remove 8 outliers----
    # Co_df <- ECU.MF.counts
    # zero_baso_id <- which(Co_df$basophils == 0)
    # Co_df <- Co_df[-zero_baso_id,]
    # Co_df$sum <- rowSums(Co_df)
    # for (i in 1:5) {
    #   Co_df[,i] <- Co_df[,i] / Co_df[,6]
    # }
    # 
    # lrlda_df <- Co_df[,1:5]
    # out.lrlda <- lrlda(lrlda_df, condition[-zero_baso_id])
    # out.lrlda$CM
    # out.lrlda$Gs
    # accuracy rate 75.078%----
    
    
    # try replace zero with small values----
    # Co_df <- ECU.MF.counts
    # zero_baso_id <- which(Co_df$basophils == 0)
    # Co_df[zero_baso_id, 5] <- 1
    # Co_df$sum <- rowSums(Co_df)
    # for (i in 1:5) {
    #   Co_df[,i] <- Co_df[,i] / Co_df[,6]
    # }
    # 
    # lrlda_df <- Co_df[,1:5]
    # out.lrlda <- lrlda(lrlda_df, condition)
    # out.lrlda$CM
    # accuracy rate 75.988%----
    
#----
    
    
# Modified LR-LDA----
  
  # ilr transformation
    ilr_df <- pivotCoord(cmultRepl(ECU.MF.counts))
    total_wbc <- rowSums(ECU.MF.counts)
    
  # just the ilr transformation 
    df0 <- ilr_df
    df0$Condition <- condition
    
    lda.model0 = lda(Condition ~ ., data = df0, center = TRUE)
    lda.model0.predict = predict(lda.model0, newdata = df0)
    tt.lda.model0 = table(Prediction = lda.model0.predict$class, True = df0$Condition)
    mean(lda.model0.predict$class == df0$Condition)
    
    lda.model0.cv = lda(Condition ~ ., data = df0, center = TRUE, CV = TRUE)
    tt.model0.cv = table(Prediction = lda.model0.cv$class, True = df0$Condition)
    mean(lda.model0.cv$class == df0$Condition)
    
    roc(response = df0$Condition, predictor = lda.model0.predict$posterior[, 2])
    
    
  # add total  
    df1 <- ilr_df
    df1$total <- total_wbc
    df1$Condition <- condition
  
    lda.model1 = lda(Condition ~ ., data = df1, center = TRUE)
    lda.model1.predict = predict(lda.model1, newdata = df1)
    tt.lda.model1 = table(Prediction = lda.model1.predict$class, True = df1$Condition)
    mean(lda.model1.predict$class == df1$Condition)
    
    lda.model1.cv = lda(Condition ~ ., data = df1, center = TRUE, CV = TRUE)
    tt.model1.cv = table(Prediction = lda.model1.cv$class, True = df1$Condition)
    mean(lda.model1.cv$class == df1$Condition)
  
  # stepwise
    df2 <- data.frame(ilr_df, ECU.MF)
    df2$total <- total_wbc
    df2$Condition <- condition
    df2 <- df2[,-6:-10]
    
    library(caret)
    # lda.model2 <- train(Condition~., data=df2, method = "stepLDA")
    # lda.model2$finalModel
    
    df_step <- data.frame(ECU.MF$platelets, "Condition" = condition)
    lda.model_step = lda(Condition ~ ., data = df_step, center = TRUE)
    lda.model_step.predict = predict(lda.model_step, newdata = df_step)
    tt.lda.model_step = table(Prediction = lda.model_step.predict$class, True = df_step$Condition)
    mean(lda.model_step.predict$class == df_step$Condition)
    
    lda.model.step.cv = lda(Condition ~ ., data = df_step, center = TRUE, CV = TRUE)
    tt.model.step.cv = table(Prediction = lda.model.step.cv$class, True = df_step$Condition)
    mean(lda.model.step.cv$class == df_step$Condition)
    
    # step-wise full for reference----
    
    df2.ref <- ECU.MF
    df2.ref$Condition <- condition
    # lda.model2.ref <- train(Condition~., data=df2.ref, method = "stepLDA")
    # lda.model2.ref$finalModel
    
    df_step_ref <- data.frame("lymphocytes" = ECU.MF$lymphocytes,
                              "hematocritP" = ECU.MF$hematocritP,
                              "Condition" = condition)
    lda.model_step_ref = lda(Condition ~ ., data = df_step_ref, center = TRUE)
    lda.model_step_ref.predict = predict(lda.model_step_ref, newdata = df_step_ref)
    tt.lda.model_step_ref = table(Prediction = lda.model_step_ref.predict$class,
                                  True = df_step_ref$Condition)
    mean(lda.model_step_ref.predict$class == df_step_ref$Condition)
    
    lda.model.step.ref.cv = lda(Condition ~ ., data = df_step_ref, center = TRUE, CV = TRUE)
    tt.model.step.ref.cv = table(Prediction = lda.model.step.ref.cv$class, True = df_step_ref$Condition)
    mean(lda.model.step.ref.cv$class == df_step_ref$Condition)
    
    # ----
    
  # manually
    df3 <- data.frame(ilr_df,ECU.MF[, c(7, 14)])
    df3$total <- total_wbc
    df3$Condition <- condition
    
    lda.model3 = lda(Condition ~ ., data = df3, center = TRUE)
    lda.model3.predict = predict(lda.model3, newdata = df3)
    tt.lda.model3 = table(Prediction = lda.model3.predict$class, True = df3$Condition)
    mean(lda.model3.predict$class == df3$Condition)
    
    lda.model3.cv = lda(Condition ~ ., data = df3, center = TRUE, CV = TRUE)
    tt.model3.cv = table(Prediction = lda.model3.cv$class, True = df3$Condition)
    mean(lda.model3.cv$class == df3$Condition)
    
  # adding variables with counts data  
    df3_ref <- data.frame(ECU.MF.counts, ECU.MF[, c(7, 14)])
    # df3_ref$total <- total_wbc
    df3_ref$Condition <- condition
    
    lda.model3_ref = lda(Condition ~ ., data = df3_ref, center = TRUE)
    lda.model3_ref.predict = predict(lda.model3_ref, newdata = df3_ref)
    tt.lda.model3_ref = table(Prediction = lda.model3_ref.predict$class, True = df3_ref$Condition)
    mean(lda.model3_ref.predict$class == df3_ref$Condition)
    
    lda.model3_ref.cv = lda(Condition ~ ., data = df3_ref, center = TRUE, CV = TRUE)
    tt.model3.cv = table(Prediction = lda.model3_ref.cv$class, True = df3_ref$Condition)
    mean(lda.model3_ref.cv$class == df3_ref$Condition)
    
    
  # only three added variables
    df4 <- data.frame(ECU.MF[, c(7, 14)])
    df4$total <- total_wbc
    df4$Condition <- condition
    
    lda.model4 = lda(Condition ~ ., data = df4, center = TRUE)
    lda.model4.predict = predict(lda.model4, newdata = df4)
    tt.lda.model4 = table(Prediction = lda.model4.predict$class, True = df4$Condition)
    mean(lda.model4.predict$class == df4$Condition)
    
    lda.model4.cv = lda(Condition ~ ., data = df4, center = TRUE, CV = TRUE)
    tt.model4.cv = table(Prediction = lda.model4.cv$class, True = df4$Condition)
    mean(lda.model4.cv$class == df4$Condition)
    
  # ROC curves
    roc_score_1 = roc(response = df1$Condition, predictor = lda.model1.predict$posterior[, 2])
    roc_score_2 = roc(response = df_step$Condition, predictor = lda.model_step.predict$posterior[, 2])
    roc_score_3 = roc(response = df3$Condition, predictor = lda.model3.predict$posterior[, 2])
    roc_score_ref = roc(response = df4$Condition, predictor = lda.model4.predict$posterior[, 2])
    roc_score_fullsteplda = roc(response = df_step_ref$Condition, predictor = lda.model_step_ref.predict$posterior[, 2])
    roc_score_WBCext = roc(response = df3_ref$Condition, predictor = lda.model3_ref.predict$posterior[, 2])
    
    roc_FS = ggroc(list("LDA Full" = roc_score_full, 
                        "LDA Full with stepwise" = roc_score_fullsteplda,
                        "LDA WBC" = roc_score_WBC,
                        "LDA with extra variables" = roc_score_ref,
                        "LDA with WBC + extra variables" = roc_score_WBCext,
                        "LR-LDA with 4LR" = roc_score_lrlda, 
                        "LR-LDA with 4LR + total" = roc_score_1,
                        "LR-LDA stepwise 4LR + full + total" = roc_score_2, 
                        "LR-LDA with 4LR + total + extra " = roc_score_3), 
                   linetype=1, size = 0.5) + 
      xlab("FPR") + ylab("TPR") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_FS
    
    roc_FS_Co = ggroc(list("LR-LDA with 4LR" = roc_score_lrlda, 
                           "LR-LDA with 4LR + total" = roc_score_1,
                           "LR-LDA stepwise 4LR + full + total" = roc_score_2, 
                           "LR-LDA with 4LR + total + extra " = roc_score_3), 
                   linetype=1, size = 1) + 
      xlab("FPR") + ylab("TPR") +
      geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed")
    roc_FS_Co
    