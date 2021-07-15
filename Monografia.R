###### Carregando Pacotes ######
library(ggplot2)
library(dplyr)
library(plm)
library(splm)
library(pder)
library(lmtest)
library(spNetwork)
library(spdep)
library(maptools)
library(fields)
library(Matrix)
library(MuMIn)
###### Criando Funções #####
SpatWeight <- function(long,lat,max) {
  W = rdist.earth(cbind(long,lat),miles=F) #Gives km distance matrix
  W[W > max] <- 0 #makes observations beyond max dist zero
  W <- ifelse(W!=0,1/W,W) #inverts non-zero distances
  for(i in 1:dim(W)[1]) {W[i,i] = 0} #because same observation is not always exactly zero, makes diagonal zero
  W <- Matrix(W,sparse=T) #Makes matrix sparse
  return(W)}

intervaloci = function(x,ci){
  obj=summary(x)
  A=NULL
  for (i in 1:length(x$coefficients)){
    A=c(A,par.avg(obj$CoefTable[i,1],obj$CoefTable[i,2],1,level = ci))
  }
  coni=t(matrix(A,ncol = 10))
  rownames(coni)=names(x$coefficients)
  colnames(coni)=names(A)[1:5]
  return(coni)
}

###### Gráficos Qualitativos ######
ev_escolas=read.csv2("C:\\Users\\evani\\Documents\\Pastas\\UEFS\\Monografia\\Dados\\ev_escolas.csv")

ev_escolas %>%
  filter(!(Característica %in% "Total"),
         Ano >= 2008) %>%
  ggplot() +
  aes(x = Ano, fill = Característica, weight = Escolas) +
  geom_bar() +
  labs(x = "", y = "", title = "", subtitle = "", caption = "", fill = "") +
  scale_fill_grey() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_text(aes(y=Escolas, label=Escolas),size=3.5, color="white",position = position_stack(0.9))

###### Coletando dados ######
coord=read.csv2("C:\\Users\\evani\\Documents\\Pastas\\UEFS\\Monografia\\CONCLUSÃO MONO\\coord.csv")
coord=coord[,c(1,3,4)]
colnames(coord)[1]="COD"

dados=read.csv("C:\\Users\\evani\\Documents\\Pastas\\UEFS\\Monografia\\CONCLUSÃO MONO\\dados_finais.csv")
dados=cbind(dados[,4],dados[,1:3],dados[,5:36])
colnames(dados)[1]="CIDADE"

dados=merge(dados,coord, by="COD")
dados=pdata.frame(dados,index = c('CIDADE','ANO'))
pdim(dados)

###### Coletando matriz de pesos #######
wpe=read.csv2("C:\\Users\\evani\\Documents\\Pastas\\UEFS\\Monografia\\CONCLUSÃO MONO\\wpe.csv")
rownames(wpe)=wpe[,1]
wpe=wpe[,-1]
colnames(wpe)=rownames(wpe)
wpe=as.matrix(wpe)


coords=data.frame(cbind(as.numeric(dados$latitude[dados$ANO==2000]),as.numeric(dados$longitude[dados$ANO==2000])))
wped=as.matrix(SpatWeight(coords$X2,coords$X1,500))
wped=wped/apply(wped,1,sum)
colnames(wped)=colnames(wpe)
rownames(wped)=rownames(wpe)

##### Correlação Espacial #####
pcdtest(dados$HOMIC, w=wped)
rwtest(dados$HOMIC, w=wped, replications = 999)

formula=HOMIC~T+TEMPO+DID+Dens+PIB_PC+SALARIO_REAL+GINI_SALARIO+EMPPOP+QUALI_EDUC
modelo=plm(formula, dados, model="within")
summary(modelo)

pcdtest(modelo$residuals, w=wped)
rwtest(modelo$residuals, w=wped, replications = 999)

##### Estimações #####

##### WDIST #####
dados$wp <- kronecker(wped, diag(1,13)) %*% dados$DID

#Lag espacial
#Efeitos fixos
felsmo=spml(formula, dados, listw = wped,lag = TRUE, effect = "individual", spatial.error = "none")
obj=summary(felsmo)
summary(felsmo)
obj$rsqr
intervaloci(felsmo,0.95)