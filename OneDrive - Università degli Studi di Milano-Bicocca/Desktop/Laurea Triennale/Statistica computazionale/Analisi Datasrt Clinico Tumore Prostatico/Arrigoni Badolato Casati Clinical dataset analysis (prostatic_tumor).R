rm(list=ls())
#Carichiamo il dataset. 

load('prostate_cancer.RData')
View(pros_canc)

#Di seguito sono riportate le library utilizzate:

library(tidyverse)
library(ggcorrplot)
library(GGally)
library(gridExtra)
library(mclust)
library(caret)
library(Rmixmod)
library(flexmix)

#Sistemiamo le variabili e visualizziamo il dataset.

colnames(pros_canc) <- c("Age","Weight","Perf.rating","Cardio.dis.hist","Syst.blood.pres","Diast.blood.pres","Electroc.code","Serum.haemoglob","Size.prim.tumor","Histologic.grade","Serum.prostatic.acid.phosphatase","Bone.metastases","Stage","Observation","SurvStat")

#Metadata:
#-Age: età del paziente al momento della rilevazione.
#-Weight: peso del paziente in chilogrammi al momento della rilevazione.
#-Performance.rating: stato di salute del paziente al momento della rilevazione. Nei dati analizzati i valori variano fra 0 e 3 ed assumono il seguente significato:
#   0 = Attività normale 
#   1 = A letto < 50% del tempo diurno 
#   2 = A letto > 50% del tempo diurno 
#   3 = Costretto a letto 

#-Cardiovascular.disease.history: storico di patologie cardiovascolari (0 assente, 1 presente).
#-Systolic.blood.pressure: pressione arteriosa sistolica al momento della rilevazione (in kPa).
#-Diastolic.blood.pressure: pressione arteriosa diastolica al momento della rilevazione (in kPa).
#-Electrocardiogram.code: codice di classificazione della salute del cuore del paziente, dove valori bassi sono associati a cuori "più sani". I valori variano fra 0 a 6 ed assumono i seguenti significati:
#   0 = Normale. 
#   1 = Benigno. 
#   2 = Disturbi ritmici e alterazioni elettrolitiche. 
#   3 = Blocchi cardiaci o difetti di conduzione.
#   4 = Sforzo cardiaco.
#   5 = Vecchio infarto miocardico.
#   6 = Infarto miocardico recente.
#-Serum.haemoglobine: livelli di emoglobina nel sangue (in g/L).
#-Size.primary.tumor: dimensione tumore primario, cioè di quello prostatico (in mm).
#Index.of.tumour.stage.and.Histologic.grade: è un indice composto dato dalla somma dello stadio e dell'indice di aggressività del tumore in scala Gleason. I valori nei dati analizzati variano fra 5 e 15 (si noti che il valore massimo è associato, per costruzione, al valore massimo di stage).
#-Serum.prostatic.acid.phosphatase: livello di fosfatasi acida sierica all'interno del sangue (in U/L). Tale acido è normalmente presente nella prostata, la sua diffusione nel sangue indica la rottura della barriera di protezione prostatica e la metastasi del tumore.
#-Bone.metastases: presenza di metastasi ossea (0 assente 1 presente).
#-Stage: stadio clinico del tumore. Nel dataset sono presenti pazienti che si trovano nello stadio 3 o 4.
#-SurvStat: definisce il codice di stato di sopravvivenza che assume valori da 0 a 9. Ogni valore dell'indice assume un particolare significato, che riportiamo per completezza:
# 0 = Vivo. 
# 1 = Morto per cancro alla prostata.
# 2 = Morto per malattia cardiaca/vascolare.
# 3 = Morto per ictus. 
# 4 = Morto per embolia polmonare. 
# 5 = Morto per altro cancro.
# 6 = Morto per malattia respiratoria. 
# 7 = Morto per altra causa specifica non-oncologica. 
# 8 = Morto per altra causa non specificata non-oncologica. 
# 9 = Morto per causa sconosciuta. 

pros_canc <- as_tibble(pros_canc)
str(pros_canc)

#Trasformiamo in factor le variabili che ne hanno le caratteristiche.

pros_canc$Perf.rating <- as.factor(pros_canc$Perf.rating)
pros_canc$Electroc.code <- as.factor(pros_canc$Electroc.code)
pros_canc$Histologic.grade <- as.factor(pros_canc$Histologic.grade)
pros_canc$Bone.metastases <- as.factor(pros_canc$Bone.metastases)
pros_canc$Stage <- as.factor(pros_canc$Stage)
pros_canc$SurvStat <- as.factor(pros_canc$SurvStat)
pros_canc$Cardio.dis.hist   <- as.factor(pros_canc$Cardio.dis.hist)
str(pros_canc)

#Rimuoviamo la variabile Observation, essa rappresenta un'ID che non aggiunge alcuna informazione ai fini dell'analisi. 

pros_canc <- select(pros_canc, -Observation)

#Creiamo un secondo dataset selezionando le variabili unicamente quantitative. 

pros_canc_quant <- pros_canc %>% select_if(is.numeric)


############################# ANALISI ESPLORATIVA################################

# Siamo interessati alla variabile Serum.prostatic.acid.phosphatase, che è una sostanza prodotta solo dalla ghiandola prostatica. Se presente nel sangue può indicare la metastasi del cancro.

summary(pros_canc$Serum.prostatic.acid.phosphatase)
sum(pros_canc$Serum.prostatic.acid.phosphatase == 0)

p <- ggplot(pros_canc)

p + geom_boxplot(aes(y = Serum.prostatic.acid.phosphatase)) +
  labs(y = "Fosfatasi acida sierica")+
  labs(title = "Boxplot fosfatasi acida sierica")

p +  geom_density(aes(x=Serum.prostatic.acid.phosphatase), fill='orange', alpha=0.2)+
  scale_x_continuous(breaks=seq(0,8, by=1))+
  scale_y_continuous(breaks=seq(0,1, by=0.05))+
  theme_minimal()+
  labs(x =" Fosfatasi acida sierica",
       y = 'Densità', title = 'Densità non parametrica della distribuzione')


# I valori assunti dalla variabile 'Serum.prostatic.acid.phosphatase' sono molto dispersi assumendo valori in un range che varia da 1 a 9999, evidenziano inoltre una grande asimmetria della distribuzione. Per questo motivo applichiamo la trasformazione logaritmica, così da rendere l'analisi più robusta.

# Boxplot e grafico sulla densità non parametrica della distribuzione (Kernel density).

p + geom_boxplot(aes(y=log(Serum.prostatic.acid.phosphatase)))+
  labs( x = "",
        y = "Log fosfatasi acida sierica",
        title = "Distribuzione del logaritmo della fosfatasi acida sierica") +
  theme_minimal()

p +  geom_density(aes(x=log(Serum.prostatic.acid.phosphatase)), fill='orange', alpha=0.2)+
  scale_x_continuous(breaks=seq(0,8, by=1))+
  scale_y_continuous(breaks=seq(0,1, by=0.05))+
  theme_minimal()+
  labs(x ="Logaritmo della fosfatasi acida sierica",
       y = 'Densità', title = 'Densità non parametrica della distribuzione')

# Da questo punto in poi, quando si parlerà di fosfatasi acida sierica, ci si riferirà al suo logaritmo.

# Condizionato rispetto allo stadio clinico del tumore:

p + geom_boxplot(aes(x = Stage,
                     y = log(Serum.prostatic.acid.phosphatase), fill =Stage
                     ), alpha = 0.5) +
  theme_minimal() +
  labs(x = "Stage",
       y = "Fosfatasi acida sierica",
       title = "Confronto dell'Fosfatasi acida sierica per stadio tumorale")

p +
  geom_density(aes(x = log(Serum.prostatic.acid.phosphatase),
                   fill = Stage),
               alpha = 0.4) +
  theme_minimal() +
  labs(
    x = "Fosfatasi acida sierica",
    y = "Densità",
    fill = "Stage",
    title = "Distribuzione dell'Fosfatasi acida sierica per stadio tumorale")

#Entrambi i grafici evidenziano una netta differenza in distribuzione della fosfatasi acida in base allo stadio tumorale del paziente. Per i pazienti allo stadio 3 la variabilità è di gran lunga minore, ed anche moda e mediana hanno un valore inferiore rispetto ai pazienti allo stadio 4, che presentano anche variabilità maggiore.

# Analisi grafica condizionata alla presenza di metastasi ossea:

p + geom_boxplot(aes(x = factor(Bone.metastases),
                     y = log(Serum.prostatic.acid.phosphatase), 
                     fill= Bone.metastases), alpha = 0.5) +
  theme_minimal() +
  labs(x = "Metastasi ossee",
       y = "Fosfatasi acida sierica",
       title = "Fosfatasi acida sierica e metastasi ossee")

p +
  geom_density(aes(x = log(Serum.prostatic.acid.phosphatase),
                   fill = factor(Bone.metastases)),
               alpha = 0.4) +
  theme_minimal() +
  labs(
    x = "Fosfatasi acida sierica",
    y = "Densità",
    fill = "Metastasi ossea",
    title = "Distribuzione della Fosfatasi acida sierica per presenza metastasi ossea")  +
  scale_fill_discrete(
    labels = c("Assente", "Presente"))

# Similmente a quanto accadeva condizionando la distribuzione della fosfatasi acida allo stadio del tumore, anche conizionando rispetto alla assenza o presenza di metastasi ossea sembra esserci un'evidente differenzazione tra le due distrubuzioni. In questo caso quando la metastasi è assente la variabilità è minore e valori mediani e modali dell'acido sono più bassi. 
# Si noti però che in questa situazione, rispetto al condizionamento in base allo stadio, i valori a discostarsi dal baffo superiore del primo boxplot sono decisamente più numerosi. Questo intuitivamente indica che vi sono più pazienti che nonostante presentino valori della fosfatasi molto grandi non hanno ancora una metastasi ossea diagnosticata.

# Istogramma che rappresenta la distribuzione dell'indice di aggressività Gleason per presenza di metastasi ossea:

table(pros_canc$Bone.metastases, pros_canc$Histologic.grade)

ggplot(pros_canc, aes(x = Histologic.grade)) +
  geom_bar(aes(fill = Bone.metastases, y =after_stat(count)/sum(after_stat(count))), col = 'grey',  position = position_dodge(preserve = "single")) + 
  theme_minimal() +
  labs(
    x = "Indice di aggressività del tumore",
    y = "presenza percentuale",
    fill = "Metastasi ossea",
    title = "Distribuzione dell'indice per presenza metastasi ossea")+
  scale_fill_discrete(
    labels = c("Assente", "Presente"))

# L'indice di aggressività è ottenuto dalla somma della scala "Gleason", che quantifica l'aggressività del tumore, con il rispettivo stadio.
# Se il nuovo tessuto è molto differente dal tessuto sano, le cellule abbandonano le loro funzioni principali e prioritarizzano l'espansione e la riproduzione delle cellule cancerogene, rendendo il tumore più aggressivo, ottenendo un punteggio della scala molto elevato. E' logico quindi che un valore di aggressività più alto comporti una maggiore probabilità della presenza di metastasi e l'analisi empirica sembra confermarlo.

# Scatterplot fra dimensione del tumore e fosfatasi acida sierica, colorato in base allo stage.

(ggp <- ggplot(data=pros_canc_quant, mapping = aes(y=log(Serum.prostatic.acid.phosphatase), x= Size.prim.tumor, color = pros_canc$Stage)) + geom_point() +
    theme_minimal() +
    labs(x= "Dimensione tumore primario",
         y="Fosfatasi acida sierica nel sangue",
         color="Stage",
         title="Relazione tra dimensione del tumore primario e fosfatasi acida sierica"))

# Dal grafico non sembra emergere una chiara ed evidente associazione tra dimensione del tumore primario e concentrazione di fosfatasi acida, che si paleserebbe con due zone ad alta densità di unità statistiche. Dal punto di vista clinico, ciò si può spiegare dal fatto che la fosfatasi acida aumenta molto in caso di metastasi ed esse si verificano anche se il tumore primario non è particolarmente grande ma in base ad altre condizioni biologiche. Inoltre potrebbe essere che il soggetto in metastasi abbia curato il tumore nella sede primaria, cioè la prostata.

# Pairs plot per stage del tumore:

pros_canc_log <- pros_canc
pros_canc_log$Serum.prostatic.acid.phosphatase <- log(pros_canc$Serum.prostatic.acid.phosphatase)
print(pros_canc_log, width = Inf)

ggpairs(pros_canc_log, 
        columns = c("Age", "Weight", "Syst.blood.pres", "Diast.blood.pres",
                    "Serum.haemoglob", "Size.prim.tumor", 
                    "Serum.prostatic.acid.phosphatase"),
        mapping = aes(color = Stage, alpha = 0.5),
        title = "Matrice di grafici a coppie per stadio del tumore") +
  theme_bw()

# Lo stadio del tumore sembra caratterizzare solo in parte una diversa distribuzione della dimensione del tumore primario. In concreto ciò accade in quanto lo stadio riflette principalmente la diffusione ed aggressività del tumore, non solo la crescita locale.

# L'unica variabile che sembra caratterizzare in maniera significativa lo stadio del carcinoma è nei fatti la fosfatasi acida sierica.
# Le altre variabili non sembrano spiegare particolarmente lo stadio. QUesta prima analisi grafica ci porta a pensare, in maniera esplorativa, che molte delle variabili a nostra disposizione non saranno di grande aiuto nel processo di classificazione e clusering.

# Si procede osservando quali sono le caratteristiche empiriche dei due gruppi :

pros_canc %>%
  group_by(Stage) %>%
  summarize(
    n = n(),
    mean = mean(log(Serum.prostatic.acid.phosphatase), na.rm = TRUE),
    sd = sd(log(Serum.prostatic.acid.phosphatase), na.rm = TRUE),
    min = min(log(Serum.prostatic.acid.phosphatase), na.rm = TRUE),
    max = max(log(Serum.prostatic.acid.phosphatase), na.rm = TRUE),
    median = median(log(Serum.prostatic.acid.phosphatase), na.rm = TRUE),
    Q1 = quantile(log(Serum.prostatic.acid.phosphatase), 0.25, na.rm = TRUE),
    Q3 = quantile(log(Serum.prostatic.acid.phosphatase), 0.75, na.rm = TRUE))



########### ANALISI DELLE CORRELAZIONI E DELLE COMPONENTI PRINCIPALI ###########

#Segue il grafico delle correlazioni; attenzione, questo grafico è disposto in modo da raggruppare le variabili con correlazione positiva alta (e bassa) vicine. Per questo motivo si perde la simmetria ed il grafico deve essere interpretato con questa consapevolezza.

corr <- cor(pros_canc_quant)

ggcorrplot(corr, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE,
           colors = c("darkblue", "white", "red3")) + 
  labs(title = "Matrice di correlazione")

# Osservando il grafico non sono presenti correlazioni elevate che potrebbero essere problematiche per le tecniche da applicare successivamente, se non per la pressione sistolica e diastolica (0.63). Altre variabili con correlazioni relativamente alte sono peso e pressione diastolica (0.23), peso e pressione sistolica (0.19) e peso e emoglobina nel siero (0.26).

# Si è poi svolta l'analisi delle componenti principali.

pca <- princomp(pros_canc_quant, cor=T)
pca_summary <- summary(pca)

#Varianza spiegata dalle prime 5 componenti:

sum((pca$sdev[1:5])^2) / 7 

#Contributo di ogni variabile a ciascuna componente principale:

pca$loadings[,1:7]

#Variabile dominante per ciascuna comp. principale:

names(pros_canc_quant)[apply(pca$loadings[,1:7], 2, function(x) which(x^2==max(x^2)))]

#Informazione contenuta da ciascuna componente:

y <- (pca_summary$sdev^2)/sum(pca_summary$sdev^2)
variance <- tibble(x = c(1:7), var = (as.numeric(y)))

#Scree plot e varianza cumulata:

scree <- ggplot(variance, aes(x, (var))) +
  geom_point() +
  geom_line(col="blue") +
  labs(x = "Componente principale",
       y = "Varianza spiegata",
       title = "Scree plot")

comulative <- ggplot(variance, aes(x, cumsum(var))) +
  geom_point() +
  geom_line(col="red") +
  labs(x = "Componente principale",
       y = "Varianza spiegata cumulata",
       title = "Varianza Cumulata")

grid.arrange(comulative, scree , ncol=2, nrow=1)

# Analizzando le componenti principali del dataset emerge un decadimento graduale della varianza spiegata dalle componenti, senza una chiara struttura gerarchica. Le prime sei componenti spiegano tra il 10% e il 26% della varianza totale. A fronte di ciò, si è deciso di non escludere alcuna variabile dalle analisi, anche in virtù della moderata quantità a disposizione (7 rispetto a 475 u.s.). In questa situazione la PCA è comunque utile in ambito esplorativo per avere contezza di come le variabili interagiscano insieme, e non prese a coppie come nel caso della correlazione che risulta invece limitante.


########################### MODEL BASED CLUSTERING #############################

# Si effettua la model based clustering analysis. Tale tecnica è 'unsupervised', cioè l'etichetta dei gruppi (nel nostro caso stadio del tumore) è "nascosta".    

# Si sceglie in un primo momento di usare la fosfatasi acida sierica come unica variabile per descrivere i gruppi, in quanto è quella che dall'analisi esplorativa caratterizza meglio il fenomeno.

# Per semplicità di codice si trasforma la variabile acid in modo permanente:

pros_canc_quant$Serum.prostatic.acid.phosphatase <- log(pros_canc_quant$Serum.prostatic.acid.phosphatase)

#Si procede con la clustering analysis: 

mclust.options("emModelNames")

Mclust.Acid <- Mclust(pros_canc_quant$Serum.prostatic.acid.phosphatase ,model="V",G=2)
summary(Mclust.Acid, parameters = TRUE)

# Si osserva che vi sono differenze evidenti tra i cluster sia in termini di media che di variabilità
# Grafico densità:

waiting.density <- densityMclust(pros_canc_quant$Serum.prostatic.acid.phosphatase, model="V",G=2) 
plot(waiting.density, what ="density", xlab = "Acido Prostatico",main = "Densità")


#Si salvano gli oggetti necessari e si procede a mostrare il disegno della mistura e delle densità delle componenti:

(mean1<-Mclust.Acid$parameters$mean[1])
(mean2<-Mclust.Acid$parameters$mean[2])
(sd1<-sqrt(Mclust.Acid$parameters$variance$sigmasq[1]))
(sd2<-sqrt(Mclust.Acid$parameters$variance$sigmasq[2]))
(p1 <- Mclust.Acid$parameters$pro [1])
(p2 <- Mclust.Acid$parameters$pro [2])
(n<-Mclust.Acid$n)   

comp.dens <- function(x, mu, sd, p){p*dnorm(x, mu, sd)}  
mixt.dens <- function(x, mu1, mu2, sd1, sd2, p){
  p*dnorm(x, mu1, sd1) + (1-p)*dnorm(x, mu2, sd2)}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
(plot <-p + stat_function(fun = mixt.dens, args = list(mu1=mean1,mu2=mean2,sd1=sd1,sd2=sd2,p=p1),lwd=1.3, col ='black')+
    xlim(0,10) + labs(title="Densità 1-dimensionale di una mistura di normali 2-componenti")+
    stat_function(fun = comp.dens, args = list(mu=mean1,sd=sd1,p=p1),color="red",lwd=1, alpha = 0.5)+
    stat_function(fun = comp.dens, args = list(mu=mean2,sd=sd2,p=p2),color="green",lwd=1, alpha = 0.5)+
    stat_function(fun = dnorm, args = list(mean1,sd1),color="orange",lwd=0.4, lty='dashed')+
    stat_function(fun = dnorm, args = list(mean2,sd2),color="green",lwd=0.4, lty='dashed'))

#Grafici incertezza e classification:

par(mfrow = c(1, 2))
plot(Mclust.Acid, what = "classification", xlab ="Fosfatasi acida sierica")
plot(Mclust.Acid, what = "uncertainty", xlab = "Fosfatasi acida sierica")
title(main = "Grafici di Incertezza e Classificazione")
par(mfrow = c(1, 1))

# La variabilità elevata della seconda componente porta il modello stimato ad attribuire al secondo gruppo le u.s. aventi come valore 0, che a livello empirico appartengono al primo. Tale situazione è un limite strutturale del modello considerato. Per questa ragione possono sorgere dei dubbi circa la validità di utilizzare modelli con distribuzioni gaussiane.

cluster.acid <- Mclust.Acid$classification

cluster_true <- pros_canc$Stage
cluster_true <- ifelse(pros_canc$Stage == 3, 1, 2)

#Conversione a factor con livelli uguali:

cluster.acid <- factor(cluster.acid, levels = c(1, 2))
stage_class <- factor(cluster_true, levels = c(1, 2))

#Matrice di confusione:

confusionMatrix(stage_class, cluster.acid)

# Si osserva che l'accuracy della classificazione, ovvero la frequenza relativa delle osservazioni assegnate al gruppo corretto di appartenenza, è del 90.5% utilizzando unicamente la variabile fosfatasi acida sierica. Il modello sembra, quindi, identificare in modo soddisfacente i due gruppi.
# Un dato che vale la pena osservare è il p-value ottenuto dal test che confronta la classificazione effettuata con il caso "banale", ovvero quello dove viene assegnata a tutte le osservazioni l'etichetta della classe più frequente (in questa analisi il gruppo 1, presente nel 57.5% dei casi). Il valore del p-value è numericamente nullo (< 2.2e-16), confermando la significatività statistica della classificazione effettuata.
# Altre informazioni importanti per l'analisi sono la sensitività e la specificità. La sensitività (96.7%) indica quanto bene il modello riconosce i pazienti del gruppo 1 (classe positiva), mentre la specificità (82.2%) misura la capacità di identificare correttamente il gruppo 2. Notiamo una specificità relativamente più bassa rispetto alla sensitività, suggerendo che il modello è eccellente nell'identificare i pazienti del primo gruppo (solo 9 falsi negativi su 273), ma tende a commettere più errori nel riconoscimento del secondo gruppo, producendo 36 falsi positivi (pazienti del gruppo 2 erroneamente classificati come gruppo 1). Il test di McNemar (p = 0.0001) conferma statisticamente questa asimmetria negli errori in una direzione. 

# Si calcola ora la distanza Fullback - Lieber

mean_a <- as.vector(Mclust.Acid$parameters$mean[1])
mean_b <- as.vector(Mclust.Acid$parameters$mean[2])
variance_a <- Mclust.Acid$parameters$variance$sigma[1]
variance_b <- Mclust.Acid$parameters$variance$sigma[2]

full_KL <-  function(mean1, mean2, variance1, variance2){
  kls <- 1/2 * t(mean1 - mean2) %*% (solve(variance1) + solve(variance2)) %*% (mean1 - mean2) + 
    1/2*(sum(diag(variance1 %*% solve(variance2) + solve(variance1) %*% variance2))) - length(mean_a)}

distanza_kl_ncluster <- as.numeric(full_KL(mean_a, mean_b, variance_a, variance_b))

# Implementando e calcolando la distanza di KL emerge che i clusters in questo caso sembrano essere ben identificabili.

# Ora si procede con la clustering n-variata, tenendo tutte le variabili.

# Si vedano tutti i modelli disponibili:

Mclust.Canc <- mclustICL((pros_canc_quant))
plot(Mclust.Canc, ylim = c( -21000,-20000))

summary(Mclust.Canc, parameters = TRUE)
Mclust.Canc

# E' stato selezionato il modello VVE, che riconosce (correttamente) 2 gruppi, utilizzando il criterio ICL. Tale criterio si ottiene dalla somma del BIC e dell'Entropia di Shannon.

# Si utilizza allora la funzione Mclust con il modello selezionato.

mod <- Mclust(pros_canc_quant, G = 2, modelNames = 'VVE')
clusters <- mod$classification
clusters <- ifelse(clusters==2, 1,2)
clusters <- as.factor(clusters)

# Si convertono sia gli esiti della nostra cluster, sia i valori empirici a factor degli stessi livelli:

stage_class <- factor(cluster_true, levels = c(1, 2))
#Matrice di confusione:

confusionMatrix(stage_class, clusters)


# Si osserva che l'accuracy della classificazione utilizzando tutte e 7 le variabili quantitative è del 73.05%, sostanzialmente inferiore rispetto al modello univariato basato sulla sola fosfatasi acida sierica (90.5%).
# Il p-value del test che confronta la classificazione con il caso banale, rimane altamente significativo (p = 1.335e-12), confermando che il modello multivariato è comunque statisticamente migliore della strategia naive. Tuttavia, il drastico calo del coefficiente Kappa, indice di accuratezza corretto per lo sbilanciamento delle classi, (da 0.80 a 0.46) rivela un peggioramento qualitativo sostanziale.
# Analizzando sensitività e specificità, notiamo un pattern opposto rispetto al modello univariato: la sensitività scende drasticamente dal 96.7% al 69.2% (84 falsi negativi contro i precedenti 9), mentre la specificità peggiora leggermente dal 82.2% al 78.2% (44 falsi positivi contro i precedenti 36). Il test di McNemar rimane significativo (p = 0.0006), ma ora l'asimmetria degli errori è invertita e meno marcata. Il modello multivariato appare più "equilibrato" ma complessivamente meno performante.
# Questi risultati, apparentemente controintuitivi, trovano un'interpretazione nell'analisi delle componenti principali. Da tale analisi la fosfatasi acida sierica  contribuisce solo modestamente alla varianza totale del dataset, ma quella porzione di variabilità è altamente concentrata nella varianza tra i gruppi, rendendola estremamente discriminante. Le altre sei variabili, pur spiegando collettivamente la gran parte della variabilità totale, catturano principalmente la varianza entro i gruppi, che non contribuisce alla separazione, e anzi introduce "rumore" nella clustering. Nel modello univariato, il 100% dell'informazione utilizzata è discriminante, mentre nel modello multivariato questa informazione si perde a causa delle variabili che aumentano la dimensionalità senza migliorare la separazione tra gruppi.

mean_a <- as.vector(mod$parameters$mean[,1])
mean_b <- as.vector(mod$parameters$mean[,2])
variance_a <- mod$parameters$variance$sigma[,,1]
variance_b <- mod$parameters$variance$sigma[,,2]

distanza_kl_ncluster <- as.numeric(full_KL(mean_a, mean_b, variance_a, variance_b))

# La distanza ottenuta non rispetta, come sappiamo, le regole di distanza metrica (per esempio non rispetta la triangolarità) e non è nemmeno standardizzata. A livello empirico però, si può affermare che un valore ottenuto sembra essere in linea con con quello ottenuto in precedenza, evidenziando una discreta separazione fra i due cluster

# Con 475 osservazioni e 7 variabili debolmente correlate tra loro, come confermato dalla PCA, si rende impossibile ridurre la dimensionalità senza perdere informazioni potenzialmente cruciali. Nei fatti però, l'utilizzo di tutte le variabili porta a instabilità numerica e variabilità, considerando che sei di esse apportano limitate informazioni discriminanti.
# In situazioni di alta dimensionalità, le distribuzioni gaussiane stimate tendono a sovrapporsi maggiormente, le distanze diventano meno discriminanti, e l'algoritmo EM fatica a identificare componenti ben separate.
# In conclusione, la fosfatasi acida sierica emerge come la variabile che cattura la vera separazione tra i due gruppi, rendendo il modello univariato non solo più performante, ma anche più stabile numericamente e interpretabile clinicamente. In ambito clinico, l'identificazione di singoli "biomarker" che catturano la separazione tra gruppi può quindi rivelarsi più efficace rispetto a complessi modelli multivariati.

# Per scrupolo si è allora pensato di provare a fare cluster utilizzando solo le variabili che riscontrano una bimodalità evidente: cioè la dimensione del tumore primario e la presenza di fosfatasi acida sierica, come illustrano dai grafici sottostanti:

ggpairs(pros_canc_quant[, c(6, 7)],
        mapping = aes(color = pros_canc$Stage, alpha = 0.5),
        title = "Grafici di dimensione del tumore e fosfatasi per stadio del tumore")

# Visualizziamo meglio le osservazioni.

ggplot(pros_canc_quant) +
  geom_point(
    aes(y = Serum.prostatic.acid.phosphatase,x = Size.prim.tumor,
        color = pros_canc$Stage)) +
  theme_minimal() +labs(
    x = "Fosfatasi acida sierica nel sangue",
    y = "Dimensione del tumore primario",
    color = "Stage",
    title = "Relazione tra fosfatasi acida sierica e dimensione del tumore primario")

#Si procede con mclust:

Mclust.Acid.Size <- mclustICL(pros_canc_quant[,c(6,7)], G=2)
Mclust.Acid.Size

#Si seleziona il modello migliore secondo il criterio ICL, cioè  di tipo "VEI":

mod.S <- Mclust(pros_canc_quant[,c(6,7)], G = 2, modelNames = 'VEI')
clusters.s <- mod.S$classification
mod.S$parameters

# Si calcola la distanza Fullback-Lieber

mean_a <- as.vector(mod.S$parameters$mean[,1])
mean_b <- as.vector(mod.S$parameters$mean[,2])
variance_a <- mod.S$parameters$variance$sigma[,,1]
variance_b <- mod.S$parameters$variance$sigma[,,2]

distanza_kl_ncluster <- as.numeric(full_KL(mean_a, mean_b, variance_a, variance_b))

# La distanza di Fullback Lieber risulta ancora nello stesso ordine delle precedenti, con un valore di 21

# Conversione a factor con livelli uguali:

clusters.s <- factor(clusters.s, levels = c(1, 2))

# Matrice di confusione:

confusionMatrix(stage_class, clusters.s)

# Anche la combinazione delle 2 variabili che sembravano spiegare in maniera migliore la bimodalità dei nostri dati non ci ha dato un'ottima cluster (rispetto alle precedenti).

# Si riporta il grafico di densità delle 2 normali bivariate, in particolare si nota che la seconda (con medie 21.6 e 3.5) ha un'elevatissima dispersione dei dati, di conseguenza è impercettibile graficamente.

Mclust.Acid.Size.D <- densityMclust(pros_canc_quant[,c(6,7)], model="VEI", G=2)

# A causa della bassa densità del secondo cluster è necessario aumentare di molto il numero di livelli così da evidenziare il gruppo molto disperso, con densità particolarmente più bassa del primo gruppo:

plotDensityMclust2(Mclust.Acid.Size.D, data = pros_canc_quant[,c(6,7)], nlevels = 40,points.pch = 1, points.col = as.factor(Mclust.Acid.Size.D$classification), points.cex = 0.8)

ggplot(pros_canc_quant, aes(x = Size.prim.tumor, 
                            y = Serum.prostatic.acid.phosphatase, 
                            col = clusters)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c("1" = "red", "2" = "blue")) +
  labs(x = "Dimensione del tumore primario",
       y = "Fosfatasi acida sierica",
       col="Cluster") +
  theme_minimal()

plot(mod.S, what = "uncertainty", xlab = "Dimensione del tumore", ylab = "Fosfatasi acida sierica")
title("Plot uncertanty")

# Di seguito visualizziamo cosa sta effettivamente accadendo, rispetto alla classificazione reale:

par(mfrow=c(1,2))
coordProj( as.data.frame(pros_canc_quant[, c(6,7)]),
           what = "classification",
           classification = stage_class,
           symbols = c(13, 18),
           col = adjustcolor(c("red", "dodgerblue2"), alpha.f = 0.5),
           sub = "(a) classificazione reale")

#Plot model based mclust classification:

coordProj(as.data.frame(pros_canc_quant[, c(6,7)]), 
          what="classification",
          classification=clusters.s,
          col=adjustcolor(c("red","dodgerblue2"), alpha.f=0.5), symbols=c(13, 18),
          sub="(b) Model-Based Clustering")

title("Confronto fra le vere etichette e con la clustering effettuata")
par(mfrow=c(1,1))

# L'elevata dispersione del secondo cluster riflette il fatto che, mentre la fosfatasi acida separa nettamente i due gruppi lungo una direzione ben definita, la dimensione del tumore primario presenta una bimodalità molto leggera che non contribuisce sostanzialmente alla discriminazione. Come già emerso dall'analisi esplorativa, la dimensione del tumore mostra sovrapposizione considerevole tra i due stadi, introducendo variabilità che "allarga" la seconda componente senza migliorare la separazione.
# Il motivo del peggioramento rispetto al modello univariato è geometricamente evidente: nello spazio bidimensionale, la separazione tra i due cluster avviene esclusivamente lungo la direzione della fosfatasi acida, mentre la direzione della dimensione del tumore aggiunge principalmente dispersione all'interno dei cluster piuttosto che separazione tra di essi. 
# Questo risultato conferma ulteriormente dhe la fosfatasi acida sierica è il biomarker dominante per la separazione tra gli stadi del tumore, e l'aggiunta di altre variabili, anche quelle che mostrano una qualche bimodalità, non migliora sostanzialmente l'analisi. La scelta di un modello univariato, che è anche più semplice, emerge quindi non solo come la soluzione più performante, ma anche come la più stabile numericamente e interpretabile clinicamente.

########################## MODEL BASED CLASSIFICATION ##########################

# Si effettua ora la model based classification, che è  una tecnica di supervised learning (le etichette sono note al modello).
# Per effettuare la classificazione si è deciso di dividere il dataset in 3 gruppi: un gruppo con circa il 60% delle u.s., il 'training set', per stimare i classificatori (modelli), un secondo gruppo con circa il 20% delle u.s., il 'validation set', su cui si utilizzano i modelli, si stima il MER e si seleziona il modello più performante, ed un rimanente 20%, il 'test set', per verificare il funzionamento del modello selezionato. 

# Si preparano i dati selezionando le variabili predittive:

X <- pros_canc[, c("Age", "Weight", "Syst.blood.pres","Serum.haemoglob", "Size.prim.tumor","Serum.prostatic.acid.phosphatase")]
X$Serum.prostatic.acid.phosphatase <- log(X$Serum.prostatic.acid.phosphatase)

y <- ifelse(pros_canc$Stage==3, 1, 2)

# Si è costruita la matrice delle variabili X e la variabile risposta y, ricodificando lo stadio (Stage) come factor.

# Si fissa il seed, per garantire la riproducibilità dei risultati e si ottiene dal dataset il 'training set' (60%).

set.seed(24)
trainIndex <- createDataPartition(y, p = 0.6, list = FALSE)
X_train <- X[trainIndex, ]
y_train <- y[trainIndex]
X_post <- X[-trainIndex, ]
y_post <- y[-trainIndex]


#Sulle rimanenti u.s. (40%) è stata effettuata la divisione in 'validation set' (20%) e 'test set' (20%).
validIndex <- createDataPartition(y_post, p = 0.5, list = FALSE)
X_valid <- X_post[validIndex, ]
y_valid <- y_post[validIndex]
X_test <- X_post[-validIndex, ]
y_test <- y_post[-validIndex]

# Si definisce ora la lista contenente i vari modelli di utilizzabili.
listmod <- c("Gaussian_pk_L_I", "Gaussian_pk_Lk_I", 
             "Gaussian_pk_L_B", "Gaussian_pk_Lk_B",
             "Gaussian_pk_L_Bk", "Gaussian_pk_Lk_Bk",
             "Gaussian_pk_L_C", "Gaussian_pk_Lk_C",  
             "Gaussian_pk_L_D_Ak_D", "Gaussian_pk_Lk_D_Ak_D",
             "Gaussian_pk_L_Dk_A_Dk", "Gaussian_pk_Lk_Dk_A_Dk",
             "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")

# Essi differiscono nella struttura della matrice di varianze e covarianze. Ci sono modelli LDA, con matrici di covarianza comune a tutti i gruppi, modelli QDA, con matrici di covarianza specifiche per ogni gruppo e dei modelli intermedi. Questo passaggio serve ad evitare di fare assunzioni a priori sulla struttura di covarianza delle componenti lasciando che siano i dati a indicarci quale modello sia migliore.


# STEP 1: utilizzando il Training set (60%) si stima il modello ottimale atraverso il metodo della Cross-Validation (CV) ripetuta:

# In primo luogo si definiscono i parametri della procedura. Identifichiamo il numero di iterazioni per ottenere una stima del MER affidabile e il numero di gruppi in cui dividere il dataset per la CV.

n_iter <- 100  
n_cv_blocks <- 10   

#Matrice per salvare i risultati dei diversi classificatori (modelli).

cv_results_matrix <- matrix(NA, nrow = n_iter, ncol = length(listmod))
colnames(cv_results_matrix) <- listmod

#Ricodifica del vettore delle risposte.

y_train <- as.integer(y_train)

#Si avvia il ciclo for di CV ripetuta, selezionando un seed diverso per ogni modello così da garantire riproducibilità ed indipendenza dei risultati.

for(iter in 1:n_iter) {
  set.seed(100 + iter*1000)  
  
  #Per ogni iterazione vengono stimate le capacità predittive di tutti i modelli tramite CV.
  
  res <- mixmodLearn(as.data.frame(X_train), y_train,
                     models = mixmodGaussianModel(listModels=listmod),
                     criterion = "CV",
                     nbCVBlocks = n_cv_blocks)
  
  #Si estrae e si salva per ogni modello il MER ottenuto dall'iterazione i-esima del ciclo.
  
  for(i in 1:length(res@results)) {
    model_name <- res@results[[i]]@model
    cv_value <- res@results[[i]]@criterionValue  
    col_idx <- which(listmod == model_name)
    cv_results_matrix[iter, col_idx] <- cv_value
  }}

# Si analizzano i modelli per vedere il migliore in termini di minimizzazione e stabilità del MER. Si calcola, per ogni modello, media e deviazione standard dell'errore di classificazione.

cv_mean <- colMeans(cv_results_matrix, na.rm = TRUE)
cv_sd <- apply(cv_results_matrix, 2, sd, na.rm = TRUE)

#Segue una tabella riassuntiva delle statistiche calcolate:

results_summary <- data.frame(
  CV_mean = round(cv_mean, 4),
  CV_sd = round(cv_sd, 4),
  stringsAsFactors = FALSE)
results_summary

#Per selezionare il modello migliore si ordinano i modelli in base all'errore medio e si seleziona quello con il valore minore. In caso di modelli con valori simili verrà selezionato il modello con variabilità minore, e nel caso di ulteriore equivalenza verrà seleziomnato il modello più semplice.

results_summary_ord <- arrange(results_summary, by=results_summary$CV_mean)
best_model_name <- rownames(results_summary_ord)[c(1,2,3)]

# Nella nostra analisi il modello migliore è stato "Gaussian_pk_L_Bk" . Tale modello ha le seguenti caratteristiche:
#-pk: le proporzioni dei gruppi sono stimate autonomamente (non è una vera informazione in quanto deriva direttamente dal dataset di partenza)
#-L: il volume della covarianza è uguale tra classi.
#-Bk: la forma della matrice di covarianza di ogni cluster è specifica mentre la sua direzione resta costante lungo gli assi.
# Utilizzando una nomenclatura equivalente si tratta di un modello EVI.

# Questo modello può essere interpretato clinicamente dal fatto che i pazienti in stadio 3 e stadio 4 presentano una variabilità complessiva simile (stesso volume complessivo), ma differiscono nella struttura delle relazioni tra le variabili, in particolar modo come abbiamo constatato anche dall'analisi esplorativa una grande discrepanza emergeva evidentemente nella misurazione della fosfatasi acida che risulta essere il biomarker.

# Di seguito sono riportati dei grafici che riassumono la selezione del modello.

#Boxplot della CV:

par(mfrow=c(1,2))
boxplot(cv_results_matrix, las=2, 
        col="lightblue",
        main="Variabilità errore per modello",
        ylab="MER (CV)",
        names=substr(listmod, 10, 30),
        cex.axis=0.7)
abline(h=min(cv_mean), col="red", lty=2)

boxplot(cv_results_matrix[, c(5,6)], las=2, 
        col="lightblue",
        main="Variabilità erorre per modello",
        ylab="MER (CV)",
        names=substr(listmod[c(5,6)], 10, 30),
        cex.axis=0.7)
abline(h=min(cv_mean), col="red", lty=2)
par(mfrow=c(1,1))

#Barplot dell'errore medio stimato e dello standard error per i vari modelli:

bp <- barplot(cv_mean, las=2,
              col="lightgreen",
              main="Errore medio per modello",
              ylab="Errore medio",
              ylim=c(0,0.37),
              names.arg=substr(listmod, 10, 30),
              cex.names=0.7)

arrows(x0=bp, y0=cv_mean - cv_sd, 
       x1=bp, y1=cv_mean + cv_sd, 
       angle=90, code=3, length=0.05, 
       col="black", lwd=1.5)

abline(h=min(cv_mean), col="red", lty=2)

# Da entrambi i grafici emerge come il modello selezionato sia il migliore in termini di CV medio (sebbene il modello pk_Lk_Bk abbia una media del tutto simile) e variabilità (che risulta essere in linea con gli algtri modelli e minore del 'concorrente'). La scelta è corroborata dal fatto che il modello pk_L_Bk sia anche il più semplice.

# STEP 2: validazione

# Si utilizza il modello selezionato nello step precedente, addestrato sul training set per ottenere delle stime del MER su dati esterni all'allenamento dei modelli.

final_train_model <- mixmodLearn(data = X_train, knownLabels = y_train,
                                 models = mixmodGaussianModel(listModels = best_model_name))

#Si classificano le osservazioni del 'validation set':

pred_valid <- mixmodPredict(data = X_valid,
                            classificationRule = final_train_model@bestResult)

#Si calcola poi il MER sul 'validation set'.

y_pred_valid <- pred_valid@partition
(mer_valid <- mean(y_pred_valid != as.numeric(y_valid)))

#Questo valore è circa pari a 0.1158. Questo significa che il 11.58% dei pazienti è classificato in maniera errata. Il valore è in linea con le stime ottenute dal MER stimato tramite CV nel training set e ci indica che il modello ha una buona capacità di distinguere pazienti in base allo stadio del carcinoma, assegnando all'88.42% delle u.s. il corretto stadio di appartenenza.

# Matrice di confusione:

confusionMatrix(as.factor(y_valid), as.factor(y_pred_valid))

# Dalla matrice di confusione emerge come il modello non abbia mai predetto in maniera errata pazienti al 3° stadio allocandoli al 4°, ma tutti gli errori sono stati nell'altra direzione. Si può quindi affermare che la precisione per lo Stage 3 (sensitività) è molto alta (a causa della minore variabilità del fenomeno), invece il modello fatica a classificare correttamente pazienti allo Stadio 4. Questo comportamento potrebbe essere causato dalla forte variabilità della seconda componente per alcune variabili come fosfatasi acida sierica e dimensione del tumore primario.

# STEP 3: 'test set' (20%)

# Si applica il modello selezionato precedentemente ai dati del test set, i quali non sono stati utiulizzati per le stime del modello migliore:

pred_test <- mixmodPredict(data = X_test, classificationRule = final_train_model@bestResult)
y_pred_test <- pred_test@partition

#Si calcola il MER del modello sul 'test set'.

(mer_test <- mean(y_pred_test != as.numeric(y_test)))

#Questo MER, essendo applicato sul test set, rappresenta la stima più realistica di errore del modello su eventuali dati nuovi. Si noti che assume un valore in linea (circa 9.47%) con il MER calcolato nello step precedente, portando a credere che il modello generalizzi bene.

#Matrice di confusione finale:

confusionMatrix(as.factor(y_test), as.factor(y_pred_test))

#Questa seconda matrice di confusione fornisce un'indicazione sulla capacità globale del nostro modello, a differenza della precedente che valutava la prestazione durante la fase di sviluppo.
# I risultati sono del tutto in liena con la precentente analisi effettuata.

#############################################################################

# Si è poi effettuata la classification utilizzando come unica variabile la fosfatasi acida sierica. Si è deciso allora di dividere il dataset in 3 gruppi come per il caso multivariato.

# Si preparano i dati selezionando la sola variabile predittiva (la fosfatasi acida sierica):

X_acid <- data.frame(Serum.prostatic.acid.phosphatase = 
                       log(pros_canc$Serum.prostatic.acid.phosphatase))

y <- as.factor(pros_canc$Stage) 
y <- ifelse(pros_canc$Stage == 3, 1, 2)

# Si è costruita la matrice X_acid contenente la variabile predittiva e la variabile risposta Stage.

# Si fissa il seed e si ottiene dal dataset il 'training set' (60%).

set.seed(242)
trainIndex <- createDataPartition(y, p = 0.6, list = FALSE)
X_train_acid <- X_acid[trainIndex, , drop = FALSE]
y_train_acid <- y[trainIndex]
X_post_acid <- X_acid[-trainIndex, , drop = FALSE]
y_post_acid <- y[-trainIndex]

# Sulle rimanenti u.s. (40%) è stata effettuata la divisione in 'validation set' (20%) e 'test set' (20%).

validIndex <- createDataPartition(y_post_acid, p = 0.5, list = FALSE)
X_valid_acid <- X_post_acid[validIndex, , drop = FALSE]
y_valid_acid <- y_post_acid[validIndex]
X_test_acid <- X_post_acid[-validIndex, , drop = FALSE]
y_test_acid <- y_post_acid[-validIndex]

# Si definisce ora la lista contenente i vari modelli utilizzabili per il caso univariato. Con una sola variabile, molti modelli collassano a casi equivalenti.

listmod_acid <- c("Gaussian_pk_L_I", "Gaussian_pk_Lk_I", 
                  "Gaussian_pk_L_B", "Gaussian_pk_Lk_B")

# STEP 1: Utilizziamo il training set (60%) per la stima del modello 'migliore' con Cross-Validation ripetuta:

# Matrice per salvare i risultati dei diversi classificatori (modelli).

cv_results_matrix_acid <- matrix(NA, nrow = n_iter, ncol = length(listmod_acid))
colnames(cv_results_matrix_acid) <- listmod_acid

# Ricodifica del vettore delle risposte:

y_train_acid <- as.integer(y_train_acid)

# Si avvia il ciclo for di CV ripetuta:

for(iter in 1:n_iter) {
  set.seed(100 + iter*1000)  
  
  #Per ogni iterazione vengono stimati tutti i modelli gaussiani utilizzando la CV.
  
  res <- mixmodLearn(X_train_acid, y_train_acid,
                     models = mixmodGaussianModel(listModels = listmod_acid),
                     criterion = "CV",
                     nbCVBlocks = n_cv_blocks)
  
  #Si estrae e si salva per ogni modello il MER ottenuto dalla CV
  
  for(i in 1:length(res@results)) {
    model_name <- res@results[[i]]@model
    cv_value <- res@results[[i]]@criterionValue  
    col_idx <- which(listmod_acid == model_name)
    cv_results_matrix_acid[iter, col_idx] <- cv_value
  }}

#Si analizzano i modelli per vedere il migliore in termini di minimizzazione e stabilità del MER. Si calcola, per ogni modello, media e deviazione standard dell'errore di classificazione.

cv_mean_acid <- colMeans(cv_results_matrix_acid, na.rm = TRUE)
cv_sd_acid <- apply(cv_results_matrix_acid, 2, sd, na.rm = TRUE)

#Segue una tabella riassuntiva delle statistiche calcolate:

results_summary_acid <- data.frame(
  CV_mean = round(cv_mean_acid, 4),
  CV_sd = round(cv_sd_acid, 4),
  stringsAsFactors = FALSE)
results_summary_acid

results_summary_ord_acid <- arrange(results_summary_ord_acid)

#Osservando i modelli si nota che i due con MER minore hanno tale valore praticamente identico e anche la deviazione standard è praticamente uguale. Per questo motivo si è deciso di selezionare il modello pk_Lk_I in quanto è un modello più semplice.

best_model_names_acid <- rownames(results_summary_ord_acid)[c(1,2)]

#Di seguito sono riportati dei grafici che riassumono la selezione del modello.

#Boxplot della CV:

boxplot(cv_results_matrix_acid, las = 2, 
        col = "lightblue",
        main = "Variabilità errore per modello (solo acido)",
        ylab = "MER (CV)",
        names = substr(listmod_acid, 10, 30),
        cex.axis = 0.8)
abline(h = min(cv_mean_acid), col = "red", lty = 2)

# Barplot dell'errore medio stimato per i vari modelli:

bp <- barplot(cv_mean_acid, las = 2,
              col = "lightgreen",
              main = "Errore medio per modello (solo acido)",
              ylab = "Errore medio",
              ylim=c(0,0.16),
              names.arg = substr(listmod_acid, 10, 30),
              cex.names = 0.8)

arrows(x0=bp, y0=cv_mean_acid - cv_sd_acid, 
       x1=bp, y1=cv_mean_acid + cv_sd_acid, 
       angle=90, code=3, length=0.05, 
       col="black", lwd=1.5)

abline(h=min(cv_mean_acid), col="red", lty=2)

#Si noti che i due modelli migliori sono quelli con varianze diverse per classe (pk_Lk). Questo risultato era prevedibile dato che la distribuzione della fosfatasi acida mostra variabilità diversa tra i due stadi.

best_model_names_acid <-  'Gaussian_pk_Lk_I'

# STEP 2: Validazione (20%)
# Si utilizza il modello selezionato nello step precedente, addestrato sul training set.

final_train_model_acid <- mixmodLearn(data = X_train_acid, 
                                      knownLabels = y_train_acid,
                                      models = mixmodGaussianModel(
                                        listModels = best_model_names_acid))

#Si classificano le osservazioni del 'validation set':

pred_valid_acid <- mixmodPredict(data = X_valid_acid,
                                 classificationRule = final_train_model_acid@bestResult)

#Si calcola poi il MER sul 'validation set'.

y_pred_valid_acid <- pred_valid_acid@partition
(mer_valid_acid <- mean(y_pred_valid_acid != as.numeric(y_valid_acid)))

# Questo valore sembra essere generalmente minore rispetto alla classificazione effettuata con tutte le variabili, indicando che la sola fosfatasi acida è un biomarker efficace per la classificazione.

# Matrice di confusione:

confusionMatrix(
  as.factor(y_valid_acid),
  as.factor(y_pred_valid_acid))

# STEP 3: 'Test set' (20%)

# Si applica il modello selezionato ai dati del 'test set':

pred_test_acid <- mixmodPredict(data = X_test_acid, 
                                classificationRule = final_train_model_acid@bestResult)
y_pred_test_acid <- pred_test_acid@partition

# Si calcola il MER del modello sul 'test set'.

(mer_test_acid <- mean(y_pred_test_acid != as.numeric(y_test_acid)))

# Questo MER risulta inferiore rispetto al 16.84% ottenuto con tutte le variabili, confermando che il modello univariato generalizza meglio ed è molto più semplice rispetto a quello multivariato.

# Matrice di confusione finale:

confusionMatrix(
  as.factor(y_test_acid),
  as.factor(y_pred_test_acid))

# Visualizzazione delle distribuzioni e classificazione

#Rispetto alla classificazione non supervisionata, model bades clustering, la classificazione supervisionata, model based classification, subisce molto meno il rumore introdotto nel caso multivariato dalle molteplici variabili non informative presenti nel dataset. Questo comportamento è interessante ma non ci stupisce in quanto il conoscere l'etichetta ha permesso alla model based classification di pesare diversamente le informazione contenute delle variabili, ed in questo modo il modello 'si accorge' da solo che il biomarker è la fosfatasi acida sierica.
# Ciò deriva direttaemnte dalla formula di log verosomigliamnza da massimizzare per una normale multivariata f(x∣μk​,Σk​)=(2π)p/2∣Σk​∣1/21​exp(−21​(x−μk​)TΣk−1​(x−μk​)), il peso sarebbe (x-μk)

# Successivamente si è deciso di procedere ad effetuare una classificazion  con modelli MDA (mixture discriminant analysis).

# Segue l'analisi con il dataset completo:

n <- nrow(pros_canc_quant)
X_quant <- pros_canc_quant
X_full <- pros_canc_quant
X_full$Stage = pros_canc$Stage

#Si impostino i parametri per implementare la CV (cross-validation):

set.seed(24)
y <- ifelse(pros_canc$Stage == 3, 1, 2)

# Successivamente si identifica il 'training set', 60% delle osservazioni.

trainIndex <- createDataPartition(y, p = 0.6, list = FALSE)
X_train <- X_quant[trainIndex, ]
y_train <- y[trainIndex]
X_post <- X_quant[-trainIndex, ]
y_post <- y[-trainIndex]

# Si identificano anche il 'validation' ed il 'test set', composti dal 20% delle osservazioni ciascuno:

validIndex <- createDataPartition(y_post, p = 0.5, list = FALSE)
X_valid <- X_post[validIndex, ]
y_valid <- y_post[validIndex]
X_test <- X_post[-validIndex, ]
y_test <- y_post[-validIndex]

G <- 4               # Numero massimo di componenti da testare
V <- 5               # Numero di fold
iter <- 10           # Numero di ripetizioni CV per stabilizzare il risultato
n_train = n*60/100   # Numero di osservazioni nel 'training set'
B <- (n*60/100)/V    # Numero di osservazioni per fold

#Di seguito creiamo le matrici per salvare gli errori:

err_i = matrix(NA ,G,V)
err = matrix(NA ,G, iter)

#Segue implementato un triplo ciclo for, quello più esterno ripete la CV, poi vi è un ciclo sul numero di componenti e infine, più internamente, vi è un ciclo sui fold.

for (i in 1:iter) {
  perm <- sample(n_train)
  for (g in 1:G) {
    err_fold <- numeric(V)
    for (v in 1:V) {
      start_idx <- (v - 1) * B + 1
      end_idx <- min(v * B, n_train)
      val_indices <- perm[start_idx:end_idx]
      X_train_fold <- X_train[-val_indices, ]
      y_train_fold <- y_train[-val_indices]
      X_val_fold <- X_train[val_indices, ]
      y_val_fold <- y_train[val_indices]
      mod <- MclustDA(X_train_fold, y_train_fold, G = g, modelNames = 'VEE')
      pred <- predict(mod, X_val_fold)$class
      err_fold[v] <- sum(pred != y_val_fold) / length(y_val_fold)}
    err[g, i] <- mean(err_fold)
  }
}

# Al termine del ciclo, sono calcolati i MER:

MER_hat <- apply(err, 1, mean)
MER_sd <- apply(err, 1, sd)

#Emergono quindi due criticità: in primo luogo la scelta del tipo di modello. Si sta infatti utilizzando un modello VEE che è stato scelto dopo aver effettuato alcuni tentativi che restano però una scelta empirica e semiarbitraria. 
# Risultava infatti tedioso ripetere per 14 volte l'analisi testando tutti i modelli. Avremmo potuto non inserire il modello in modo esplicito e farlo scegliere al modello, ma, in quel caso avremmo avuto problemi nella ripetizione dell'esperimento per ottnere risultati consistenti. 
# Anche in virtù dei dati a disposizione, dove l'applicazione di modelli MDA non sembra essere utile in quanto si hanno poche unità statistiche a disposizione e relazioni che non presentano multimodalità marcate.
# Inoltre anche la scelta di G è problematica: si sta infatti supponendo che ogni cluster sia formato da una mistura di G elementi, ma nella realtà potrebbe esssere anche che cluster diversi necessitino un numero diverso di elementi, quindi in realtà si sta facendo un'approssimazione.
# In ultima istanza i modelli testati sembravano indicare che la scelta in corrispondenza di G=1, cioè la classification semplice, fosse paragonabile con i modelli più complessi (risultati diversi in base al seed).
# Per tutte queste ragioni si è deciso di non estendere ulteriormente l'analisi.

#Si sono selezionati i modelli migliori, cioè quelli con G = 1 e G = 2. Si valutino ora per il 'validation set'.

#Vengono create le matrici per salvare predizioni e MER sul 'validation set'.

pred_valid <- matrix(NA, nrow= 95, ncol= 2)
mer_valid <- matrix(NA, nrow= 1, ncol= 2)

# Di seguito si inizializza il ciclo sui due valori di G scelti:

for (g_opt in c(1,2)){
  final_model <- MclustDA(X_train, y_train, G = g_opt, modelNames = 'VEE')
  pred_valid[1:95,g_opt] <- predict(final_model, X_valid)$class}

#Viene riportato il calcolo del MER sul 'validation set'.

colnames(mer_valid) <- c('G=1','G=2')

mer_valid[1,1] <- mean(pred_valid[,1] != y_valid)
mer_valid[1,2] <- mean(pred_valid[,2] != y_valid)
mer_valid

# Il 'validation set' sembra portare a preferire il modello con G=2 anche se con un margine molto fine. Per questa ragione si è pensato di testare entrambi i modelli con il 'test set'.

#Ecco la classificazione sul 'test set', per entrambi i modelli.

#Modello con G=1:

final_model_1 <- MclustDA(X_train, y_train, G = 1, modelNames = 'VEE')
pred_test_1 <- predict(final_model_1, X_test)$class

#Calcolo MER sul 'test set'

mer_test_1 <- mean(pred_test_1 != y_test)

# Modello con G=2:

final_model_2 <- MclustDA(X_train, y_train, G = 2, modelNames = 'VEE')
pred_test_2 <- predict(final_model_2, X_test)$class

# Calcolo MER sul 'test set':

mer_test_2 <- mean(pred_test_2 != y_test)


#Matrice di confusione sul 'test set':

conf_test <- confusionMatrix(as.factor(y_test), as.factor(pred_test_2))
print(conf_test)

#A livello di MER i modelli con G=1 e G=2 sembrano produrre dei risultati del tutto simili (0.1073684 e 0.1150877 rispettivamente per questo seed nel training set). Dunque si è optato per scegliere il modello più semplice, G=1.
#In pratica ci siamo ricondotti al caso di classificazione semplice con un modello VEE. L'analisi effettuata è risultata interessante, ma, i dati a nostra disposizione non necessitano di una complicazione aggiuntiva che risulta anche meno interpretabile.




# MODEL BASED CLUSERING WITH COVARIATES#################### ---------------



#Si procede con l'implementazione della regressione MEM (mixture Expert Model). Per la selezione delle variabili su cui appliccare questa tecnica si è in un primo momento utilizzato il grafico di ggpair dall'analisi esplorativa. Si selezionano la coppia di variabili che presentava un'interessante differenziazione del valore della correlazione in base al gruppo di appartenenza (Stage). Si sono quindi scelte le variabili fosfatasi acida sierica e l'emoglobina nel siero, la coppia di variabili che mostrava correlazione parziale pari a 0.11 per il gruppo 1 (Stage = 3) e -0.21 per il gruppo 2 (Stage=4). 
# Si riporta di seguito lo scatter plot della coppia di variabili ed il calcolo della correlazione.

ggplot(pros_canc_quant) +
  geom_point(
    aes(y = Serum.prostatic.acid.phosphatase,x = Serum.haemoglob,
        color = pros_canc$Stage)) +
  theme_minimal() +labs(
    y = "log fosfatasi acida sierica nel sangue",
    x = "emoglobina nel siero",
    color = "Stage",
    title = "Relazione tra fosfatasi acida sierica e emoglobina nel siero")

# Osserviamo che le correlazioni sono diverse per i gruppi. In particolare, per il gruppo Stage = 4 risulta addirittura esserci una correlazione negativa fra la quantità di siero emoglobina e fosfatasi acida 

cor(pros_canc_quant$Serum.prostatic.acid.phosphatase[pros_canc$Stage==3], pros_canc_quant$Serum.haemoglob[pros_canc$Stage==3])

cor(pros_canc_quant$Serum.prostatic.acid.phosphatase[pros_canc$Stage==4], pros_canc_quant$Serum.haemoglob[pros_canc$Stage==4])

# Si adopera la funzione flexmix per la stima di un modello di regressione.

ModFlex <- flexmix(Serum.prostatic.acid.phosphatase~Serum.haemoglob ,data=pros_canc_quant ,k=2)
parameters(ModFlex)
ICL(ModFlex)

# La differenza fra la KLdiv del pacchetto Flexmix e la nostra "custom" distanza di Kullback leiber sembrano dare risultati su scala diversa, ma questo in realtà ha senso, poichè la nostra valuta le distanze fra i parametri di media tenendo conto delle varianze, mentre qui osserviamo le distanze fra le osservazioni (475 in tutto)

(fullback_library <- sum(KLdiv(ModFlex)))

summary(ModFlex) 
test.param<-refit(ModFlex)
summary(test.param)

# Il modello stimato evidenzia aspetti interessanti. come si poteva immaginare dall'analisi esplorativa il coefficente associato al parametro gamma(2,1) risulta essere essere largamente significativo, ad indicare quindi il fatto che la variabile fosfatasi acida sierica sia influenzata negativamente da un aumento dell'emoglobina. Si noti che il secondo cluster risulta essere molto più disperso come evidenziato anche dal fatto che lo standard error associato al parametro gamma(2,1) sia 4 volte maggiore rispetto a quello del primo cluster.

plot(ModFlex)

# Analizzando il grafivo risulta importnante notare che le probabilità di appartenenza al gruppo sono altamente influenzate dalla variabilità dei coefficienti. In particolare, vi sono diverse u.s. con probailità di appartenenza (zj) al secondo gruppi vicine ad 1 (cioè situazioni di sicurezza statistica rispetto alla cluster). Ciò invece non accade per la prima componente, infatti si notano moltissime u.s con probabilità di appartenenza al primo gruppo che sono alte ma non 'certe', e di conseguenza vi sono altrettante u.s. etichettate come molto poco probabili ma non nulle dalla seconda componente. Questo comportamento è spiegato dalla variabilità del secondo cluster, che avendo 'code' della distrubuzione molto pesanti tende ad avere valori di densità non nulle anche se ci si trova mlto lontani dal baricentro. Questa osservazione era stata del tutto simile qunando si era analizzata la model based clustering con la variabile acid.

ggplot(data=pros_canc_quant, mapping =  aes(y = Serum.prostatic.acid.phosphatase,x = Serum.haemoglob,
                                            color = pros_canc$Stage)) +
  geom_point()+ 
  geom_smooth(method="lm", se=F, size=1)


#In questo ultimo paragrafo, si vuole modellare la variabile fosfatasi acida sierica in funzione delle altre variabili quantitative. Si cre inizialmente un modello lineare (scelta effettuata in base all'analisi esplorativa della variabile acid) nullo.
# Si è deciso di utilizzare una procedure stepwise forward per l'individuazione delle tre variabili più significative. Abbiamo deciso di scegliere le prime tre variabili prima di osservare i valori dei p-value ottenuti dalla procedura, in quanto si è consapevoli dei limiti statistici legati alla riproducibilità dei risultati e dai valori non affidabili risultanti dalla ripetizione di test statistici, implicita nei criterion based selection processes.

mod0 <- lm(Serum.prostatic.acid.phosphatase ~ 1, data = pros_canc_quant)
summary(mod0)

# Successivamente si utilizza la procedura stepwise forward. Si opta per il test F in quanto il modello preso in considerazione è lineare. 

mod_acid <- step(object = mod0  , scope = list( lower = ~1, upper = ~  Age + Weight + Syst.blood.pres + Diast.blood.pres + Serum.haemoglob + Size.prim.tumor ), direction = "forward", steps = 3, trace = T, data = pros_canc_quant, test = "F")
summary(mod_acid)

# Dalla summary di quest'ulitmo modello si nota che la la variabile "Weight" non risulta statisticamente significativa con un test effettuato al 5%. Tuttavia si decide di mantenerla nel modello poiché l'interesse principale è valutarne l'effetto all'interno dei gruppi e non in generale.

# Dopo aver fissato il seed, si utilizza il modello individuato all'interno della funzione flexmix.

set.seed(111) 

fit <- flexmix(Serum.prostatic.acid.phosphatase ~ Size.prim.tumor + 
                 Serum.haemoglob + Weight, data = pros_canc_quant , k=2)
summary(fit)

# Visualizziamo i parametri stimalti del modello:

parameters(fit)

 # Da questi valori emerge in maniera chiara quanto osservato anche durante le precedenti analisi: è presente un cluster (comp. 1) con valori medi (intercetta) della risposta più elevati e un'elevata variabilità; troviamo poi un secondo gruppo (comp. 2) che ha valori medi più bassi e una variabilità minore. Il primo gruppo corrisponderà ai pazienti allo stadio 4, mentre il secondo allo stadio 3.

# Di seguito si riporta il valore dell'ICL. Naturalmente questo valore preso da solo non è molto esplicativo ma si utilizzerà successivamente per valutare l'inizializzazione dell'algoritmo EM.

ICL(fit)     


# Si analizzano i parametri stimati dal modello, diviso per cluster, con i p-value associati.

fit2 <- refit(fit)
summary(fit2)

# Dall’analisi dei modelli stimati emerge che il cluster che include prevalentemente pazienti con Stage 4, presenta una situazione in cui diverse variabili esplicative risultano essere significative statisticamente (size ed emoglobina). Questo risultato era emerso anche durante l'analisi esplorativa, nonostante  in quel caso le relazioni fossero valutate a coppie. Si era infatti notato che il secondo gruppo denotava una maggiore correlazione con le variabili rispetto al primo. Naturalmente non era scontato che l'analisi a coppie si riflettesse sul modello completo in quanto le correlazioni tra variabili esplicative e risposta non implicano correlazioni fra di esse ed i residui di un modello di regressione.
# Analizzando i risultati di questpo modello, risulta che la variabile weight non apporti una significativa informazione circa la distinzione fra gruppi per la variabile fosfatasi acida sierica. 

#Etichette suggerite dal modello: 

fit@cluster

#I grafici del modello:

plot(fit)

# Analizzando il grafico risulta importnante notare che le probabilità di appartenenza al gruppo sono altamente influenzate dalla variabilità dei coefficienti. In particolare, vi sono diverse u.s. con probailità di appartenenza (zj) al secondo gruppi vicine ad 1 (cioè situazioni di sicurezza statistica rispetto alla cluster). Ciò invece non accade per la prima componente, si notano moltissime u.s con probabilità di appartenenza al primo gruppo che sono alte ma non 'certe', e di conseguenza vi sono altrettante u.s. etichettate come molto poco probabili ma non nulle dalla seconda componente. Questo comportamento è spiegato dalla variabilità di quest'utlima, che avendo 'code' della distrubuzione molto pesanti tende ad avere valori di densità non nulle anche se ci si trova mlto lontani dal baricentro. Questa osservazione era stata del tutto simile quando si era analizzata la model based clustering con la variabile acid.

# Per valutare la distanza fra cluster del modello è stata calcolata la divergenza di Kullback-Leibler. Risulta però importante ricordare che tale misura non è standardizzata quindi va interpretata con cautela. Inoltre, la distanza calcolata dal comando KLdiv è empirica, cioè calcolata dalla formula sum(P* log(P/Q)), il valore tenderà a infito nonappena vi siano una coppia di valori del tipo (0.9999, 0.0001). Un valore molto grade di tale indice non significa che le distribuzioni siano ben distinte in assoluto, ma soltanto in quel punto. Vi sono quindi evidenti limiti sull'interpretabilità di questo risultato.

KLdiv(fit)

# Si evidenzia quanto detto precedentemente rispetto al grafico, cioè che la componente 1, per le caratteristiche descritte, ha un valore decisamente inferiore rispetto alla seconda dove sono presenti più osservazioni con probabilità vicine a zero.


bestfit <- stepFlexmix( Serum.prostatic.acid.phosphatase ~ Size.prim.tumor + Serum.haemoglob + Weight, data = pros_canc_quant, k = 1:4,  nrep = 30,verbose = TRUE, drop = F, unique = FALSE)

# Il modello con K = 2 sembra essere di fatto il migliore. E' interessante osservare come malgrado k = 3 abbia un BIC molto simile viene penalizzato dalla la creazione di un nuovo cluster probabilmente inutile e  sovrapposto ad un altro, rendendo l'Entropia molto elevata.

best_model <- getModel(bestfit, which = "ICL")  
MEM_labels <- clusters(best_model)
MEM_labels <- as.factor(MEM_labels)

cluster_true <- pros_canc$Stage
cluster_true <- ifelse(cluster_true==3, 1,2)
cluster_true <- as.factor(cluster_true)

confusionMatrix(cluster_true, MEM_labels)

# Si analizzi ora la confusion matrix. In particolare, il valore dell'indice kappa è circa 0.80, il che è alto, indicando una clustering molto soddisfacente. Rispetto alla model based clustering, che nel nostro caso evidenziava diversi problemi come scritto in precedenza (attibuzione di molti falsi negativi), la tecnica Expert Model sembra essere più bilancia, ottenendo falsi negativi e falsi positivi e egual misura, facendoci pensare che i problemi nopn siano più di tipo prettamente strutturale del modello ma intrinseci nella variabilità dei dati. Questo è evidenziato anche dal McNemar's Test, il quale p-value non permette di rifiutare l'ipotesi nulla: probabilità falsi negativi = probabilità falsi positivi (mentre nelle clustering fatte precedentemente il p-value era significativo indicando asimmetria negli errori).


#In seguito a tutte le analisi svolte, si può affermare che all'interno dei dati esiste una struttura latente, che il modello a due componenti descrive in modo ottimale, essendo in grado di spiegare l'eterogeneità presente nel dataset.

########## GRAZIE PER L'ATTENZIONE - Gruppo Fortinati ###############

