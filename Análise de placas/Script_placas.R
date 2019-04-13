######################## II Avaliação MLG #############################
##################### Arthur C. Rocha RA 94361 ########################

#Pacotes (todos necessários para o funcionamento das funções)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(hnp)

#------------Função de diagnóstico geral (binomial)--------------------

diag_binom=function(fit.model){
  X <- model.matrix(fit.model)
  n <- nrow(X)
  p <- ncol(X)
  w <- fit.model$weights
  W <- diag(w)
  H <- solve(t(X)%*%W%*%X)
  H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
  h <- diag(H)
  ts <- resid(fit.model,type="pearson")/sqrt(1-h)
  td <- resid(fit.model,type="deviance")/sqrt(1-h)
  di <- (h/(1-h))*(ts^2)
  a <- max(td)
  b <- min(td)
  x11(title = paste("Diagnosticos do",deparse(substitute(fit.model))))
  par(mfrow=c(2,2))
  plot(fitted(fit.model),h,xlab="Valores Ajustados", ylab="Medida h",
       main="Pontos de Alavanca", pch=16,panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  # identify(fitted(fit1.model), h, n=1)
  #
  plot(di,xlab="Índice", ylab="Distância de Cook",
       main="Pontos Influentes",pch=16,panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  # identify(di, n=1)
  #
  plot(td,xlab="Índice", ylab="Resíduo Componente do Desvio",
       main="Pontos Aberrantes", ylim=c(b-1,a+1), pch=16,
       panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  abline(2,0,lty=2,col="orange")
  abline(-2,0,lty=2,col="orange")
  # identify(td, n=1)
  #
  plot(predict(fit.model),td,xlab="Preditor Linear", 
       ylab="Residuo Componente do Desvio",
       main="Função de Ligação", ylim=c(b-1,a+1), pch=16,
       panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
       frame.plot = F,col.axis="gray70")
  axis(1,col = "gray70",col.ticks ="gray70")
  axis(2,col = "gray70",col.ticks ="gray70")
  abline(2,0,lty=2,col="orange")
  abline(-2,0,lty=2,col="orange")
  
  #envelope
  x11(title = paste("Envelope do",deparse(substitute(fit.model))))
  par(mfrow=c(1,1))
  fit.model %>% hnp(halfnormal = F,xlab="Quantis Teóricos",
                  ylab='Componente de Desvio',pch=16)
  grid(lty=1)
  
  #resultados
  cat("-------------------Resultados do modelo-------------------")
  sumario=fit.model %>% summary()
  x11(title = paste("Resultados do",deparse(substitute(fit.model))))
  grid.arrange(bottom=paste(sumario$call[2],""),top=paste("Resultados do",deparse(substitute(fit.model))),
               tableGrob(rbind(sumario$coefficients,
                               deviance=c(sumario$deviance,NA,NA,NA),
                               AIC=c(sumario$aic,NA,NA,NA))%>% round(digits = 3)))
  
  
}

#---------------- Fim função diagnóstico ------------------------------


#----------------------Leitura de dados--------------------------------

dados=read.table(header = T,text = "
calor imersao inadequado teste
                  7 1 5 19
                 14 1 4 16
                 27 1 9 22
                 51 1 1 5
                 7 1.7 2 8
                 14 1.7 2 8
                 27 1.7 9 21
                 51 1.7 5 16
                 7 2.2 0 7
                 14 2.2 7 27
                 27 2.2 2 9
                 51 2.2 1 6
                 7 2.8 0 6
                 14 2.8 0 7
                 27 2.8 2 14
                 51 4 3 8
                 7 4 0 10
                 14 4 3 9
                 27 4 1 8
                 ")

head(dados)
plot(dados)

#######Criando variável de proporção
dados$prop=dados$inadequado/dados$teste

#-------------------------Fim leitura de dados------------------------

#--------------------------Análise exploratória------------------------

##Obsevando a quantidade de placas inadequadas
plot(table(dados$inadequado),type="h",
     panel.first=grid(col = rgb(.2,.2,.2,.1),lty=1),
     frame.plot = F,lwd=3.3,col="orange",ylab="Frequência",
     xlab="Placas Inadequadas")
axis(1,col = "gray70",col.ticks ="gray70",at = 0:9)
axis(2,col = "gray70",col.ticks ="gray70",at = 0:4)

#OBS: Bimodal, com bastante 0


### Densidade da proporção de placas inadequadas
#Gráfico de densidade
p1=ggplot(dados,aes(prop))+geom_density(colour=NA,fill="orange",alpha=.5)+
  theme_minimal()+ylab("Densidade")+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank())
#Boxplot
p2=ggplot(dados,aes(y=prop,x=""))+
  geom_boxplot(alpha=.4,fill="orange")+theme_minimal()+
  theme(legend.background = element_blank(),
        legend.position = "none",axis.title.y =element_blank(),
        axis.text.y = element_blank(),axis.ticks = element_blank(),
        plot.margin=margin(t=-0.3,b=.4,unit="cm"))+
  coord_flip()+ylab(expression(pi))

#Juntando densidade e boxplot
p3=plot_grid(p1,p2,nrow = 2,align = "v",
          rel_heights = c(3/4,1/4))


### Observando proporção de acordo com o calor

ggplot(dados,aes(factor(calor),prop))+geom_point()+theme_minimal()+
  xlab("Calor")+ylab(expression(pi))


### Observando proporção de acordo com a imersão

ggplot(dados,aes(factor(imersao),prop))+geom_point()+theme_minimal()+
  xlab("Imersão")+ylab(expression(pi))


### Procurando Interação entre covariáveis

ggplot(dados,aes(factor(imersao),prop,col=factor(calor)))+geom_point()+
  theme_minimal()+geom_line(aes(group=factor(calor)))+xlab("Imersão")+
  ylab(expression(pi))+scale_color_discrete(name = "Calor")

##------------------------Fim análise exploratória ---------------------


#-----------------------------MODELAGEM------------------------------- 

####Modelos binomiais

##aditivo com ligação canônica
#modelo
modelo1=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
      dados$calor+dados$imersao,family = binomial(link = "logit"))

#diagnóstico
diag_binom(modelo1)  


##Com interação
#modelo
modelo1.1=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
                dados$calor*dados$imersao,family = binomial(link = "logit"))

#diagnóstico
diag_binom(modelo1.1)

  
#########Testando função de ligação (canônica)

lp=modelo1$linear.predictors^2  ##preditor linear ao quadrado

#modelo com lp somado
modelo2=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$imersao+dados$calor+lp,family = binomial(link = "logit"))
#diagnóstico
diag_binom(modelo2)

########## Testando outras funções de ligação ######################

##Modelo com ligação probit
modelo3=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$calor*dados$imersao,family = binomial(link = "probit"))
#diagnóstico
diag_binom(modelo3)


##Modelo com ligação Cauchit
modelo4=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$calor*dados$imersao,family = binomial(link = "cauchit"))
#diagnóstico
diag_binom(modelo4)

##Modelo com ligação log
modelo5=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$calor*dados$imersao,family = binomial(link = "log"))

#diagnóstico
diag_binom(modelo5)


#Modelo com ligação complemento log-log
modelo6=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$calor*dados$imersao,family = binomial(link = "cloglog"))

#diagnóstico
diag_binom(modelo6)

###Em todos os gráficos o comportamento do resíduo é similar

###################### Modelo de Quasi binomial ######################
#modelo
modelo7=glm(cbind(dados$inadequado,dados$teste-dados$inadequado)~
              dados$calor*dados$imersao,family = quasibinomial)
#diagnóstico
diag_binom(modelo7)

#--------------------------Fim modelagem -----------------------------