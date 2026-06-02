library(quantreg)

result <- data.frame()
for (i in 1:length(unique(all_SES_pair_long_subset$plotID_full))){
  subset_data <- subset(all_SES_pair_long_subset,plotID_full==unique(all_SES_pair_long_subset$plotID_full)[[i]])
  
  if (nrow(subset_data) < 100) next
  m <- rq(obs.cor~tl_dist*phy.cor,data=subset_data,tau=c(0.25,0.5,0.75))
  m_summary <- summary(m,se="nid")
  
  result_table<- rbind(cbind(m_summary[[1]]$coefficients,tau=m_summary[[1]]$tau),
                       cbind(m_summary[[2]]$coefficients,tau=m_summary[[2]]$tau),
                       cbind(m_summary[[3]]$coefficients,tau=m_summary[[3]]$tau))
  
  result_table <- as.data.frame(result_table)
  result_table <- cbind(result_table,factor=rep(rownames(m_summary[[1]]$coefficients),length(m_summary)))
  result <- rbind(result,data.frame(result_table,
                                    plotID = unique(all_SES_pair_long_subset$plotID_full)[[i]],
                                    n=nrow(subset_data)))
}

#result$sig <- ifelse(sign(result$lower.bd) == sign(result$upper.bd),1,0)
result$sig <- ifelse(result$Pr...t.. < 0.05,1,0)

library(tidyverse)
summary_result <- result %>% group_by(factor,tau) %>% summarize(sig=sum(sig)/length(unique(result$plotID)))  

library(lqmm)
lqmm(obs.cor~tl_dist*phy.cor,data=all_SES_pair_long_subset,group=~site_code,random=~tl_dist*phy.cor)

int_df <- subset(result,factor=="tl_dist:phy.cor")

summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(int_df,tau == 0.5)))
summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(int_df,tau == 0.25)))
summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(int_df,tau == 0.75)))
