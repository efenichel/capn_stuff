#R code to recreate the Fenichel and Abbott 2014 figure for Nick Hanley
#Eli Fenichel, Yale University, August 2019

rm(list = ls()) #clear memory
#get libraries
if (!require("pacman")) install.packages("pacman")
p_load(capn, ggplot2, cowplot)


#pull in Gulf of Mexico Data, which is included with capn library.
#you will need to click through a few times. 
demo("GOM")

devAskNewPage(ask = FALSE) #turn off a pause feature in the demo

#gets the data out of the demo results and organizes it for graphing in ggplot
gom.out <- as.data.frame(cbind(GOMSimV$stock/1e6,  GOMSimV$shadowp, GOMSimV$vfun))
colnames(gom.out)<-c("stock", "sp", "vfun")

annuity <- as.data.frame(cbind(
  gom.out$stock,
  sapply(gom.out$stock*1e6, dwds,  Z = param)/ (10*param$delta)
))
colnames(annuity) <- c("stock", "sp")

msy.data<-as.data.frame(rbind(c(0.5*param$k/1e6,0), c(0.5*param$k/1e6,14)))
colnames(msy.data) <- c("stock", "sp")

cal.data<-as.data.frame(rbind(c(8.633e7/1e6,0), c(8.633e7/1e6,14)))
colnames(cal.data) <- c("stock", "sp")


#plot just the accounting price in ggplot
ggplot() +
  geom_line(data = gom.out, aes(x = stock, y = sp),
            color = 'black') +
  geom_line(data = annuity, aes(x = stock, y = sp),
            color = 'gray', linetype = "dashed") +
  geom_line(data = msy.data, aes(x = stock, y = sp),
            color = 'black', linetype = 'dotted') +
  geom_line(data = cal.data, aes(x = stock, y = sp),
            color = 'black', linetype = 'dotted') +
  geom_text() + 
    annotate("text", 
             label = "shadow price", 
             x = 300, 
             y = 2.7, 
             size = 3, 
             color = "black")+
  geom_text() + 
  annotate("text", 
           label = "Annuity based 
accounting price
in units $10", 
           x = 300, 
           y = 9.6, 
           size = 3, 
           color = "black")+
  geom_text()+
    annotate("text", 
           label = "MSY stock level", 
           x = 225, 
           y = 12, 
           size = 3, 
           color = "black")+
  geom_text()+
    annotate("text", 
           label = 
                "Calibration 
stock level", 
           x = 115, 
           y = 12, 
           size = 3, 
           color = "black")+
  labs(
    x= "Fish stock in millions of pounds",
    y = "Natural capital accounting price, USD")  +
  #xlim(0,25) +
  ylim(0, 14)+
  theme(  #http://ggplot2.tidyverse.org/reference/theme.html
    axis.line = element_line(color = "black"), 
    panel.background = element_rect(fill = "transparent",color = NA),
    plot.background = element_rect(fill = "transparent",color = NA)
  )
