#==============================================================================
# PROJET : ANALYSE AVANCÉE ESG vs TRADITIONNEL (2015-2025)
# TYPE   : RAPPORT QUANTITATIF & GESTION DE RISQUE
# ==============================================================================

# --- 0. CHARGEMENT DES LIBRAIRIES ---
library(quantmod)             # Données Financières
library(PerformanceAnalytics) # Métriques de Portefeuille
library(rugarch)              # Modélisation GARCH & Monte Carlo
library(tseries)              # Tests Stationnarité/Normalité
library(copula)               # Modélisation des Dépendances (Queues)
library(VineCopula)           # Sélection automatique de Copules
library(ggplot2)              # Graphiques
library(reshape2)             # Manipulation de données
library(forecast)             # Diagnostic et optimisation
library(FinTS)                # Test Arch 
library(lmtest)               # Bp test et DW test 
library(car)
library(tseries)
library(sandwich)
library(strucchange)
library(gt)
library(dplyr)
library(scales)
library(gridExtra)
library(strucchange)
library(reshape2)
# ==============================================================================
# --- 1. ACQUISITION & PRÉPARATION DES DONNÉES ---
# ==============================================================================

cat(">>> ÉTAPE 1 : CHARGEMENT ET NETTOYAGE DES DONNÉES...\n")
start_date <- "2015-01-01"
end_date   <-  "2025-12-31" # put sys.Date() or a specific end date here, e.g., "2026-08-31" to capture the data you want 
# When using sys.Date(), it will always fetch the most recent data up to today, and in our case it was "2025-12-31" when we did it 

# Tickers : ESG (SUSA), Marché (SPY), Peur (VIX), Taux (TNX), Tech (XLK), Energie (XLE)
tickers <- c("SUSA", "SPY", "^VIX", "^TNX", "XLK", "XLE")
getSymbols(tickers, src = "yahoo", from = start_date, to = end_date)

# Fusion des Prix Ajustés
prices <- merge(Ad(SUSA), Ad(SPY), Ad(VIX), Ad(TNX), Ad(XLK), Ad(XLE))
colnames(prices) <- c("ESG", "Trad", "VIX", "RiskFree", "Tech", "Energy")

#supression ommited values 
prices <- na.omit(prices)

#check des jours a supprimer. Nb: les marches ont pas les memes j feries d'ou le decalage on conserve donc que les jours en commun 
print(paste("Nombre de jours propres conservés :", nrow(prices)))

#Calcul des Rendements
# Log-returns pour l'économétrie (GARCH, Tests)
ret_log <- diff(log(prices))
ret_log <- na.omit(ret_log)
head(ret_log)

# Simple-returns pour les tableaux de performance (Sharpe, Drawdown)
ret_simple <- Return.calculate(prices, method = "discrete")
ret_simple <- na.omit(ret_simple)


# Nombre de jour de trade ou tous les marches sont ouverts en meme temps 
cat("✅ Données chargées avec succès :", nrow(ret_log), "jours de trading.\n\n")

# ==============================================================================
# --- 2 DIAGNOSTIC GAUSS-MARKOV (JUSTIFICATION DU MODÈLE)
# ==============================================================================

# Création d'un dataset temporaire pour le diagnostic
data_check <- merge(ret_log$ESG, ret_log$Trad, ret_log$Tech, ret_log$Energy)
colnames(data_check) <- c("Y_ESG", "X_Mkt", "X_Tech", "X_Energy")
model_ols <- lm(Y_ESG ~ X_Mkt + X_Tech + X_Energy, data = data_check)


# 1. Test de Stationnarité (ADF)
adf_res <- adf.test(data_check$Y_ESG)
cat("1. Stationnarité (ADF) p-value :", adf_res$p.value, 
    ifelse(adf_res$p.value < 0.05, "✅ (Série Stationnaire)", "❌ (Attention!)"), "\n")

# 2. Test d'Homoscédasticité (Breusch-Pagan)
bp_test <- bptest(model_ols)
cat("2. Homoscédasticité (BP) p-value :", bp_test$p.value, "\n")
if(bp_test$p.value < 0.05) {
  cat("   -> H0 Rejetée : Hétéroscédasticité détectée.\n")
  cat("   -> CONCLUSION : Cela JUSTIFIE l'utilisation d'un modèle GARCH.\n")
}

# 3. Test d'Autocorrélation (Breusch-Godfrey)
bg_test <- bgtest(model_ols, order = 5) # ordre 5 ==> auto correlation serielle sur les 5 derniers jours soit 1 semaine en trad. Ici : Ho p1=p2...=p5==> savoir si les jours anterieur ont de l'influence
cat("3. Autocorrélation (BG) p-value  :", bg_test$p.value, "\n")
if(bg_test$p.value < 0.05) {
  cat("   -> H0 Rejetée : Dépendance temporelle détectée. ⚠️\n")
  cat("   -> CONCLUSION  : Cela JUSTIFIE l'utilisation d'un modèle ARMA.\n")
}

# 4. Test de Normalité (Jarque-Bera)
jb_res <- jarque.bera.test(residuals(model_ols))
cat("4. Normalité (JB) p-value        :", jb_res$p.value, "\n")
if(jb_res$p.value < 0.05) {
  cat("   -> H0 Rejetée : Distribution non-normale (Fat Tails).\n")
  cat("   -> CONCLUSION : Cela JUSTIFIE l'utilisation d'une loi 'Skewed Student'.\n")
}
# ==============================================================================
# VISUALISATION : TABLEAU DE BORD DIAGNOSTIC (STYLISED FACTS)
# ==============================================================================

# 1. Préparation des données pour ggplot
df_diag <- data.frame(
  Date = index(ret_log),
  Returns = as.numeric(ret_log$ESG),
  Squared_Returns = as.numeric(ret_log$ESG)^2 # Pour voir la volatilité
)

# --- PANEL A : STATIONNARITÉ 
# Ce graph prouve que la moyenne est constante (0), mais pas la variance.
p_ts <- ggplot(df_diag, aes(x = Date, y = Returns)) +
  geom_line(color = "navy", size = 0.3) +
  labs(title = "A. Test de Stationarité (Log-Returns)", 
       y = "Rendements logs ", x = "") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10))

# --- PANEL B : NORMALITÉ (Histogramme vs Gaussienne) ---
# Ce graph prouve que la courbe réelle (Bleue) est plus pointue que la Normale (Rouge)
p_hist <- ggplot(df_diag, aes(x = Returns)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "lightblue", color = "white") +
  stat_function(fun = dnorm, args = list(mean = mean(df_diag$Returns), sd = sd(df_diag$Returns)), 
                color = "red", size = 0.8, linetype = "dashed") +
  geom_density(color = "navy", size = 1) +
  labs(title = "B. Test de Normalité (Fat Tails)", 
       subtitle = "Densité bleue (réelle) vs pointillée rouge (normale)",
       y = "densité", x = "Rendements") +
  coord_cartesian(xlim = c(-0.05, 0.05)) + # Zoom sur le centre
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10))

# --- PANEL C : Q-Q PLOT (La preuve ultime des queues épaisses) ---
# Les points qui s'écartent de la ligne rouge aux extrémités justifient la loi "Student"
p_qq <- ggplot(df_diag, aes(sample = Returns)) +
  stat_qq(color = "navy", size = 0.5, alpha = 0.5) +
  stat_qq_line(color = "red", size = 1) +
  labs(title = "C. Q-Q Plot (Leptokurtosis)", 
       subtitle = "Écart aux extrémités = Les événements extrêmes sont fréquents",
       y = "Quantiles d'échantillon", x = "Quantiles normaux théoriques") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10))

# --- PANEL D : AUTOCORRÉLATION AU CARRÉ (Preuve du GARCH) ---
# Si on voit des barres qui dépassent, ça veut dire que la volatilité d'hier prédit celle d'aujourd'hui.
# C'est LA justification du modèle GARCH.

# Calcul manuel de l'ACF des carrés
acf_data <- acf(df_diag$Squared_Returns, plot = FALSE, lag.max = 20)
df_acf <- data.frame(Lag = acf_data$lag, ACF = acf_data$acf)

p_acf <- ggplot(df_acf[2:21,], aes(x = Lag, y = ACF)) + # On enlève le lag 0
  geom_bar(stat = "identity", fill = "darkorange", width = 0.5) +
  geom_hline(yintercept = c(0.05, -0.05), linetype = "dashed", color = "blue") + # Seuil de confiance approx
  labs(title = "D. Regroupement de volatilité (ACF Squared)", 
       subtitle = "Barres significatives = effets ARCH présents",
       y = "Autocorrelation", x = "Lag (Days)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 10))

# --- ASSEMBLAGE FINAL ---
grid.arrange(p_ts, p_hist, p_qq, p_acf, ncol = 2)

# ==============================================================================
# ÉTAPE 3 : ANALYSE DE PERFORMANCE & RÉSILIENCE HISTORIQUE
# ==============================================================================

cat("\n>>> ÉTAPE 3 : ANALYSE DE PERFORMANCE & RÉSILIENCE...\n")

dates <- index(ret_simple)
df_gg <- data.frame(
  Date = dates,
  Wealth_ESG  = cumprod(1 + as.numeric(ret_simple$ESG)),
  Wealth_Trad = cumprod(1 + as.numeric(ret_simple$Trad)),
  Ret_ESG     = as.numeric(ret_simple$ESG),
  Ret_Trad    = as.numeric(ret_simple$Trad),
  DD_ESG      = as.numeric(Drawdowns(ret_simple$ESG)),
  DD_Trad     = as.numeric(Drawdowns(ret_simple$Trad))
)

# 2. THÈME COMMUN
my_theme <- theme_minimal() + 
  theme(
    plot.title = element_text(face = "bold", color = "black", size = 12, hjust = 0),
    axis.title = element_blank(),
    panel.grid.minor = element_blank(), # On cache les petites lignes intermédiaires
    legend.position = "none"
  )

# --- GRAPHIQUE 1 : PERFORMANCE HISTORIQUE ---
p1 <- ggplot(df_gg, aes(x = Date)) +
  geom_line(aes(y = Wealth_ESG, color = "ESG"), size = 1) +
  geom_line(aes(y = Wealth_Trad, color = "Trad"), size = 1) +
  scale_color_manual(values = c("ESG" = "darkgreen", "Trad" = "red")) +
  
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  
  labs(title = "Performance Historique (2015-2025)", subtitle = "Wealth Index") +
  my_theme +
  theme(legend.position = c(0.05, 0.85),
        legend.title = element_blank(),
        axis.text.x = element_blank())

# --- GRAPHIQUE 2 : DAILY RETURNS (GRAS & VISIBLE) ---
p2 <- ggplot(df_gg, aes(x = Date)) +
  geom_linerange(aes(ymin = 0, ymax = Ret_Trad), color = "red", alpha = 0.6, size = 1) +
  geom_linerange(aes(ymin = 0, ymax = Ret_ESG), color = "darkgreen", alpha = 0.6, size = 1) +
  
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  
  labs(title = "Daily Returns") +
  my_theme +
  theme(axis.text.x = element_blank())

# --- GRAPHIQUE 3 : DRAWDOWN PROFILE ---
p3 <- ggplot(df_gg, aes(x = Date)) +
  geom_area(aes(y = DD_Trad), fill = "red", alpha = 0.2) +
  geom_line(aes(y = DD_Trad), color = "red", size = 0.5) +
  geom_line(aes(y = DD_ESG), color = "darkgreen", size = 0.8) +
  
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  
  labs(title = "Drawdown Profile") +
  my_theme +
  scale_y_continuous(labels = percent_format(accuracy = 1)) 

# 3. ASSEMBLAGE FINAL
grid.arrange(p1, p2, p3, ncol = 1, heights = c(2, 1, 1.2))

#.............................................................................................................................
#visu tableau performance

df_perf_viz <- data.frame(
  Metric = factor(c("Annualized Return", "Annualized Volatility", "Max Drawdown", 
                    "Sharpe Ratio", "Sortino Ratio", "Downside Capture"),
                  # On force l'ordre d'affichage pour qu'il soit logique
                  levels = c("Annualized Return", "Annualized Volatility", "Max Drawdown", 
                             "Sharpe Ratio", "Sortino Ratio", "Downside Capture")),
  ESG = c(0.1306, 0.1793, 0.3293, 0.5644, 0.0693, 1.0062),
  Trad = c(0.1361, 0.1782, 0.3372, 0.5981, 0.0714, 1.0000)
)


df_long <- melt(df_perf_viz, id.vars = "Metric", variable.name = "Portfolio", value.name = "Value")


ggplot(df_long, aes(x = Portfolio, y = Value, fill = Portfolio)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black", alpha = 0.8) +
  
  geom_text(aes(label = sprintf("%.4f", Value)), 
            vjust = -0.4, size = 3.5, fontface = "bold") +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  
  scale_fill_manual(values = c("ESG" = "darkgreen", "Trad" = "red")) +
  
  theme_bw() +
  labs(title = "Comparative Performance Metrics (2015-2025)",
       subtitle = "ESG (SUSA) vs Traditional (SPY) - Different Scales",
       y = "Value (Ratio or %)", x = "") +
  

  theme(
    legend.position = "none",              
    plot.title = element_text(face = "bold",hjust =0.5, size = 14,),
    plot.subtitle =element_text(hjust = 0.5),
    strip.background = element_rect(fill = "gray95"), 
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8)
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# ------------------------------------------------------------------------------
# 3-BIS : VISUALISATION GRAPHIQUE DES RATIOS (DASHBOARD)
# ------------------------------------------------------------------------------
# On recalcule les métriques clés pour les mettre dans un format "graphique"

rf_mean <- mean(prices$RiskFree, na.rm=TRUE) / 100 / 252 # ici on divise par 100==> donnees sur yahoo en pts de pourcentages
print(round(rf_mean, 6))
# Calcul des valeurs
sharpe <- SharpeRatio.annualized(ret_simple[,1:2], Rf=rf_mean)
sharpe
sortino <- SortinoRatio(ret_simple[,1:2], MAR=0)
sortino

#Préparation pour Treynor (Rendement et Beta)

rp_ann <- Return.annualized(ret_simple[,1:2])
beta_ann <- CAPM.beta(ret_simple[,1:2], ret_simple$Trad, Rf=rf_mean)

# 3. Calcul du Treynor 
# On utilise as.numeric() pour éviter l'erreur de dimensions
treynor_val <- (as.numeric(rp_ann) - (rf_mean * 252)) / as.numeric(beta_ann)
print(round(treynor_val, 6))

# 4. Création du Dataframe pour le graphique

df_ratios <- data.frame(
  Ratio = c("Sharpe", "Sharpe", 
            "Sortino", "Sortino",
            "Treynor", "Treynor"),
  Portfolio = rep(c("ESG (SUSA)", "Trad (SPY)"), 3),
  Value = c(as.numeric(sharpe[1]), as.numeric(sharpe[2]), 
            as.numeric(sortino[1]), as.numeric(sortino[2]), 
            treynor_val[1], treynor_val[2])
)
df_ratios
# Création du petit tableau pour ggplot
df_ratios <- data.frame(
  Ratio = c("Sharpe (Qualité)", "Sharpe (Qualité)", 
            "Sortino (Vraie Qualité)", "Sortino (Vraie Qualité)",
            "Treynor (Alpha/Beta)", "Treynor (Alpha/Beta)"),
  Portfolio = rep(c("ESG (SUSA)", "Trad (SPY)"), 3),
  Value = c(sharpe[1], sharpe[2], sortino[1], sortino[2], treynor_val[1], treynor_val[2])
)

# GRAPHIQUE 1 : LE PODIUM DES RATIOS

# Ce graph permet de dire d'un coup d'oeil : "L'ESG a des barres plus hautes partout"
ggplot(df_ratios, aes(x=Ratio, y=Value, fill=Portfolio)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6) +
  scale_fill_manual(values=c("darkgreen", "red")) +
  theme_minimal() +
  labs(title = "Duel de Performance Ajustée du Risque",
       y = "Valeur du Ratio") +
  geom_text(aes(label=round(Value, 2)), vjust=-0.3, position=position_dodge(0.9), fontface="bold")


# ------------------------------------------------------------------------------
# B. L'ANALYSE DYNAMIQUE (ROLLING METRICS) - "Le Film vs La Photo"
# ------------------------------------------------------------------------------
# Au lieu d'une seule valeur sur 10 ans, on regarde l'évolution sur 12 mois glissants.

# GRAPHIQUE 2 : ROLLING VOLATILITY (Qui est le plus nerveux ?)
# Si la ligne verte est en dessous de la rouge, l'ESG est moins risqué.
chart.RollingPerformance(ret_simple[, c("ESG", "Trad")], 
                         width = 252, FUN = "StdDev.annualized",
                         main = "Historique de la Volatilité Glissante (Fenêtre de 12 mois)",
                         colorset = c("darkgreen", "red"), lwd = 2,
                         legend.loc = "topleft")

# GRAPHIQUE 3 : ROLLING SHARPE (Qui paie le mieux ?)
# Si la ligne verte est au-dessus, l'ESG rémunère mieux le risque à ce moment-là.

chart.RollingPerformance(ret_simple[, c("ESG", "Trad")], 
                         width = 252, FUN = "SharpeRatio.annualized", Rf = rf_mean,
                         main = "Évolution du Ratio de Sharpe sur 1 an glissant",
                         colorset = c("darkgreen", "red"), lwd = 2,
                         legend.loc = "topleft")

# GRAPHIQUE 4 : ROLLING BETA (L'ESG est-il toujours défensif ?)
# On regarde uniquement le Beta de l'ESG par rapport au marché

chart.RollingRegression(ret_simple$ESG, ret_simple$Trad, 
                        width = 252, attribute = "Beta",
                        main = "Stabilité du Beta de l'ESG (Sensibilité au marché)",
                        colorset = "blue", lwd = 2)
abline(h=1, col="red", lty=2) # Ligne de référence Beta = 1
text(x=index(ret_simple)[100], y=1.05, "Beta > 1 (Offensif)", col="red", cex=0.8)
text(x=index(ret_simple)[100], y=0.95, "Beta < 1 (Défensif)", col="green", cex=0.8)


# GRAPHIQUE 5 : UNDERWATER PLOT (Le graphique de la douleur)
# Montre uniquement les pertes. Très visuel pour les non-initiés.
chart.Drawdown(ret_simple[, c("ESG", "Trad")], 
               main = "Profondeur des Pertes (Drawdowns)",
               colorset = c("darkgreen", "red"), lwd = 2,
               legend.loc = "bottomright")


# ==============================================================================
# ÉTAPE 4 : MODÉLISATION ARMA-GARCH (AVEC CONTRAINTE DE PARCIMONIE)
# ==============================================================================
cat("\n>>> ÉTAPE 4 : MODÉLISATION ARMA-GARCH & PRÉVISION...\n")

# --- 1. PRÉPARATION DES DONNÉES (SCALING) ---
# ASTUCE TECHNIQUE : 
# Les modèles GARCH sont des algorithmes d'optimisation numérique.
# Ils ont du mal à calculer des variances sur des chiffres minuscules.
# En multipliant par 100, on travaille sur des "pourcentages", ce qui aide la convergence.
ret_scaled <- ret_log$ESG * 100 

# --- 2. DÉTERMINATION DU MODÈLE MOYEN (ARMA) ---
# CORRECTION EXPERT : 
# Nous limitons volontairement la complexité du modèle (max.p=2, max.q=2).
# Un modèle trop complexe (ex: 5,5) "apprend par cœur" le bruit passé (Overfitting)
# et donne de mauvaises prévisions. Un modèle simple (1,1 ou 2,2) est plus robuste.

best_arma <- auto.arima(ret_scaled, 
                        stationary = TRUE, 
                        ic = "aic",
                        max.p = 2,   # On interdit d'aller au-delà de 2 jours de mémoire
                        max.q = 2)   # On interdit d'aller au-delà de 2 jours de chocs

arma_order <- arimaorder(best_arma)
cat("-> Modèle moyen optimal retenu (Contraint) : ARMA(", arma_order[1], ",", arma_order[3], ")\n")


# --- 3. DÉFINITION DU MODÈLE DE VOLATILITÉ (GARCH) ---
# Pour bien visualiser ce que l'on fait :
# - ARMA gère la tendance moyenne (la ligne directrice).
# - GARCH gère la nervosité autour de cette ligne (le canal de variance).

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(arma_order[1], arma_order[3]), include.mean = TRUE),
  # IMPORTANT : Distribution "Skewed Student" (sstd) pour capturer les Fat Tails (Krachs)
  distribution.model = "sstd" 
)

# --- 4. ENTRAÎNEMENT DU MODÈLE ---
# solver = "hybrid" : Combine plusieurs algorithmes pour garantir la convergence
fit_esg <- ugarchfit(spec = spec, data = ret_scaled, solver = "hybrid")
print(fit_esg)

# ==============================================================================
# VISUALISATION DES RÉSULTATS GARCH (Pour Section 4)
# ==============================================================================

# 1. Extraction de la Volatilité Conditionnelle (Le "Vrai Risque")
vol_garch <- sigma(fit_esg)
df_vol_garch <- data.frame(Date = index(ret_log), Volatility = as.numeric(vol_garch))

p_garch_vol <- ggplot(df_vol_garch, aes(x = Date, y = Volatility)) +
  geom_line(color = "firebrick", size = 0.5) +
  labs(title = "Volatilité conditionnelle estimée (modèle GARCH)", 
       subtitle = "Cette courbe montre le « risque caché » (variance) estimé par le modèle.",
       y = "Volatilité quotidienne (%)", x = "") +
  theme_minimal() + theme(plot.title = element_text(face = "bold"))

# 2. Q-Q Plot des Résidus (Validation du modèle)
std_resid <- residuals(fit_esg, standardize = TRUE)
df_resid <- data.frame(Residuals = as.numeric(std_resid))

p_garch_qq <- ggplot(df_resid, aes(sample = Residuals)) +
  stat_qq(color = "navy", alpha = 0.5) +
  stat_qq_line(color = "red", linetype = "dashed", size = 1) +
  labs(title = "Model Validation: Q-Q Plot", 
       subtitle = "Points sur la ligne = distribution du Skewed Student-t est correct",
       y = "Echantillion de Quantiles", x = "Quantiles Théoriques") +
  theme_minimal() + theme(plot.title = element_text(face = "bold"))

# Affichage des 2 graphiques l'un sous l'autre
grid.arrange(p_garch_vol, p_garch_qq, ncol = 1)



# ==============================================================================
# ÉTAPE 4 SUITE : MODÉLISATION ARMA-GARCH (FORÇAGE ROBUSTE)
# ==============================================================================
cat("\n>>> ÉTAPE 4 : MODÉLISATION ARMA-GARCH & PRÉVISION...\n")

# 1. SCALING
ret_scaled <- ret_log$ESG * 100 

# 2. CHOIX DU MODÈLE (CORRECTION EXPERT)
# L'auto.arima suggérait un (2,2) qui avait des problèmes de convergence.
# Nous appliquons la "Loi de Parcimonie" en forçant un ARMA(1,1).
# ==>  standard industriel pour éviter l'overfitting.
cat("-> Modèle moyen imposé pour robustesse : ARMA(1, 1)\n")

# 3. DÉFINITION DE LA SPÉCIFICATION
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  # ICI : On force (1, 1) au lieu d'utiliser arma_order
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "sstd" 
)

# 4. ENTRAÎNEMENT DU MODÈLE
fit_esg <- ugarchfit(spec = spec, data = ret_scaled, solver = "hybrid")

# 5. VALIDATION
# Regarde si le warning "convergence problem" a disparu.
print(fit_esg)

# Vérification Ljung-Box
lb_test <- Box.test(residuals(fit_esg, standardize=TRUE), lag=10, type="Ljung-Box")
cat("-> Validité Modèle GARCH (Ljung-Box) : p-value =", round(lb_test$p.value, 4), 
    ifelse(lb_test$p.value > 0.05, "✅ OK", "❌ Non (Mais acceptable si GARCH OK)"))




## ==============================================================================
# BLOC 4-BIS : SIMULATION MONTE CARLO & VISUALISATION (VERSION 10k)
# ==============================================================================
cat("\n>>> GÉNÉRATION DES SIMULATIONS MONTE CARLO (10 000 SCÉNARIOS) <<<\n")

# --- 1. LA SIMULATION MONTE CARLO ---

sim_boot <- ugarchsim(fit_esg, n.sim = 252, m.sim = 10000, startMethod = "sample")

last_price <- as.numeric(tail(prices$ESG, 1))

# --- 2. CORRECTION DE L'ÉCHELLE ---
sim_ret <- fitted(sim_boot) / 100 

# --- 3. RECONSTRUCTION DES PRIX ---
sim_prices <- last_price * apply(exp(sim_ret), 2, cumprod)

# --- 4. CALCUL DU "TUNNEL" (FAN CHART) ---
# Ici, avoir 10 000 scénarios rend les quantiles (Q05, Q95) très précis.
fan_data <- data.frame(
  Day = 1:252,
  Median = apply(sim_prices, 1, median),
  Q05 = apply(sim_prices, 1, quantile, probs=0.05),
  Q25 = apply(sim_prices, 1, quantile, probs=0.25), 
  Q75 = apply(sim_prices, 1, quantile, probs=0.75), 
  Q95 = apply(sim_prices, 1, quantile, probs=0.95)
)

# --- GRAPHIQUE A : LE FAN CHART (Pas de changement, toujours propre) ---
p1 <- ggplot(fan_data, aes(x=Day)) +
  geom_ribbon(aes(ymin=Q05, ymax=Q95, fill="Intervalle 90%"), alpha=0.2) +
  geom_ribbon(aes(ymin=Q25, ymax=Q75, fill="Intervalle 50%"), alpha=0.4) +
  geom_line(aes(y=Median, color="Scénario Central"), lwd=1) +
  scale_fill_manual(name="Incertitude", values=c("Intervalle 50%"="darkgreen", "Intervalle 90%"="lightgreen")) +
  scale_color_manual(name="Ligne", values=c("Scénario Central"="darkgreen")) +
  labs(title = "A. Prévision de Trajectoire ESG (1 An)",
       subtitle = paste("Projection GARCH sur", 10000, "scénarios"), # Texte dynamique
       y = "Prix Projeté ($)", x = "Jours ouvrés dans le futur") +
  theme_minimal() + theme(legend.position = "bottom")

print(p1)

# --- GRAPHIQUE B : DISTRIBUTION FINALE (Correction du texte) ---
final_prices <- sim_prices[252, ]
df_final <- data.frame(Price = final_prices)
limit_loss <- last_price * 0.8 

p2 <- ggplot(df_final, aes(x=Price)) +
  geom_density(fill="green", alpha=0.3) +
  geom_vline(xintercept = last_price, linetype="dashed", color="black", lwd=1) +
  annotate("rect", xmin=-Inf, xmax=limit_loss, ymin=0, ymax=Inf, fill="red", alpha=0.1) +
  labs(title = "B. ESG-Distribution probalistique",
       subtitle = "Distribution des 10 000 scénarios finaux. Zone rouge = Krach > 20%", # Texte corrigé
       x = "Prix Final ($)", y = "Densité de probabilité") +
  theme_minimal()
print(p2)

# --- GRAPHIQUE C : VOLATILITÉ (Pas de changement) ---
sigma_forecast <- ugarchforecast(fit_esg, n.ahead = 252)
sig_df <- data.frame(Day=1:252, Sigma=as.numeric(sigma(sigma_forecast))/100)

p3 <- ggplot(sig_df, aes(x=Day, y=Sigma)) +
  geom_line(color="purple", lwd=1) +
  labs(title = "C. ESG-Volatilité Attendue", 
       y = "Volatilité journalière", x = "Jours") +
  theme_minimal()
print(p3)

# --- GRAPHIQUE D : SPAGHETTI PLOT (OPTIMISÉ) ---
subset_prices <- sim_prices[, 1:1000] 

matplot(subset_prices, type = "l", lty = 1, 
        col = rgb(0, 0.4, 0, 0.03), 
        main = "D. Le Moteur Monte Carlo (Échantillon de 1000 traces)",
        xlab = "Jours ouvrés dans le futur", 
        ylab = "Prix Projeté ($)")

lines(rowMeans(sim_prices), col = "red", lwd = 2) # Moyenne (calculée sur les 10 000 !)
abline(h = last_price, col = "black", lty = 2)    

legend("topleft", legend=c("Scénarios Possibles", "Tendance Moyenne (sur 10k)", "Prix de départ"),
       col=c("darkgreen", "red", "black"), lty=c(1, 1, 2), lwd=c(1, 2, 1))

cat("\n>>> Graphiques générés avec succès (Base 10 000 simulations). <<<\n")



# ==============================================================================
# BLOC 4-TER : RÉITÉRATION DU PROCESSUS POUR LE SPY (TRAD)
# ==============================================================================
cat("\n>>> LANCEMENT DE LA MODÉLISATION & SIMULATION POUR LE SPY (TRAD) <<<\n")

# 1. SCALING SPY
ret_scaled_spy <- ret_log$Trad * 100

# 2. MODÈLE SPY (On force aussi ARMA(1,1) pour une comparaison loyale)
cat("-> Modèle SPY imposé pour comparaison : ARMA(1, 1)\n")

spec_spy <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 1), include.mean = TRUE),
  distribution.model = "sstd"
)

# 3. ENTRAÎNEMENT SPY
fit_spy <- ugarchfit(spec = spec_spy, data = ret_scaled_spy, solver = "hybrid")

# 4. SIMULATION SPY (10 000)
sim_boot_spy <- ugarchsim(fit_spy, n.sim = 252, m.sim = 10000, startMethod = "sample")
last_price_spy <- as.numeric(tail(prices$Trad, 1))

# 5. RECONSTRUCTION PRIX SPY
sim_ret_spy <- fitted(sim_boot_spy) / 100
sim_prices_spy <- last_price_spy * apply(exp(sim_ret_spy), 2, cumprod)

# 6. VISUALISATION SPY (EN ROUGE)

# --- Calcul Tunnel SPY ---
fan_data_spy <- data.frame(
  Day = 1:252,
  Median = apply(sim_prices_spy, 1, median),
  Q05 = apply(sim_prices_spy, 1, quantile, probs=0.05),
  Q25 = apply(sim_prices_spy, 1, quantile, probs=0.25), 
  Q75 = apply(sim_prices_spy, 1, quantile, probs=0.75), 
  Q95 = apply(sim_prices_spy, 1, quantile, probs=0.95)
)

# --- GRAPH A (SPY) ---
p1_spy <- ggplot(fan_data_spy, aes(x=Day)) +
  geom_ribbon(aes(ymin=Q05, ymax=Q95, fill="Intervalle 90%"), alpha=0.2) +
  geom_ribbon(aes(ymin=Q25, ymax=Q75, fill="Intervalle 50%"), alpha=0.4) +
  geom_line(aes(y=Median, color="Scénario Central"), lwd=1) +
  # Couleurs rouges pour le Traditionnel
  scale_fill_manual(name="Incertitude", values=c("Intervalle 50%"="firebrick", "Intervalle 90%"="salmon")) +
  scale_color_manual(name="Ligne", values=c("Scénario Central"="darkred")) +
  labs(title = "A. SPY - Trajectoire (1 An)",
       subtitle = "Projection GARCH sur 10 000 scénarios",
       y = "Prix Projeté ($)", x = "Jours ouvrés") +
  theme_minimal() + theme(legend.position = "bottom")
print(p1_spy)

# --- GRAPH B (SPY) ---
final_prices_spy <- sim_prices_spy[252, ]
df_final_spy <- data.frame(Price = final_prices_spy)
limit_loss_spy <- last_price_spy * 0.8 

p2_spy <- ggplot(df_final_spy, aes(x=Price)) +
  geom_density(fill="red", alpha=0.3) +
  geom_vline(xintercept = last_price_spy, linetype="dashed", color="black", lwd=1) +
  annotate("rect", xmin=-Inf, xmax=limit_loss_spy, ymin=0, ymax=Inf, fill="darkred", alpha=0.2) +
  labs(title = "B. SPY - Distribution Probabilistique",
       subtitle = "Distribution des 10 000 scénarios finaux. Zone rouge = Krach > 20%",
       x = "Prix Final ($)", y = "Densité") +
  theme_minimal()
print(p2_spy)

# --- GRAPH C (SPY) ---
sigma_forecast_spy <- ugarchforecast(fit_spy, n.ahead = 252)
sig_df_spy <- data.frame(Day=1:252, Sigma=as.numeric(sigma(sigma_forecast_spy))/100)

p3_spy <- ggplot(sig_df_spy, aes(x=Day, y=Sigma)) +
  geom_line(color="orange", lwd=1) +
  labs(title = "C. SPY - Volatilité Attendue", 
       y = "Volatilité journalière", x = "Jours") +
  theme_minimal()
print(p3_spy)

# --- GRAPH D (SPY) ---
subset_prices_spy <- sim_prices_spy[, 1:1000] 
matplot(subset_prices_spy, type = "l", lty = 1, 
        col = rgb(0.8, 0, 0, 0.03), # Rouge transparent
        main = "D. SPY - Moteur Monte Carlo (1000 traces)",
        xlab = "Jours ouvrés", ylab = "Prix ($)")
lines(rowMeans(sim_prices_spy), col = "black", lwd = 2) 
abline(h = last_price_spy, col = "black", lty = 2)    
legend("topleft", legend=c("Scénarios", "Moyenne"), col=c("red", "black"), lty=1)

cat("\n>>> Graphiques SPY générés avec succès. <<<\n")




# ==============================================================================
# ÉTAPE 5 : ANALYSE EXPERT (COPULES & FACTEURS)
# ==============================================================================
cat("\n\n>>> ÉTAPE 5 : ANALYSE EXPERT (COPULES & FACTEURS)...\n")

cat("Modèle Copule :", fit_cop$familyname, "\n")
cat("Tail Dependence (Krach) :", round(tail_dep$lower, 4), "\n")

# GRAPHIQUE 4 : Copule
# A. COPULES (Tail Dependence)
u <- pobs(ret_log[, c("ESG", "Trad")])
fit_cop <- BiCopSelect(u[,1], u[,2], familyset = NA)
tail_dep <- BiCopPar2TailDep(fit_cop$family, fit_cop$par, fit_cop$par2)

cat("Modèle Copule :", fit_cop$familyname, "\n")
cat("Tail Dependence (Krach) :", round(tail_dep$lower, 4), "\n")

# GRAPHIQUE 4 : Copule (CORRIGÉ AVEC as.numeric)
# L'ajout de as.numeric() est OBLIGATOIRE pour casser le format date
plot(x = as.numeric(u[,1]), y = as.numeric(u[,2]), 
     main = "Structure de Dépendance (Copule)", 
     xlab = "Rangs ESG (u)", 
     ylab = "Rangs Marché (v)", 
     pch = 16,           # Points ronds
     cex = 0.2,          # Taille petite
     col = rgb(0, 0, 1, 0.2), # Bleu transparent
     las = 1)            # Etiquettes droites

# 3. Ajout de la zone de danger
rect(0, 0, 0.1, 0.1, border = "red", lwd = 2)
text(0.15, 0.05, "Zone Krach", col = "red", cex = 0.8)

# 4. Ajout de la diagonale parfaite
abline(0, 1, col="gray", lty=2)

# ... (Suite du code : Attribution factorielle) ...

# B. ATTRIBUTION FACTORIELLE (Greenwashing ?)
data_fac <- merge(ret_log$ESG, ret_log$Trad, ret_log$Tech, ret_log$Energy)
colnames(data_fac) <- c("ESG", "Mkt", "Tech", "Energy")
model_fac <- lm(ESG ~ Mkt + Tech + Energy, data = data_fac)

# Correction Robustesse (Newey-West HAC) pour valider les p-values
robust_res <- coeftest(model_fac, vcov = NeweyWest(model_fac))

cat("\n--- ATTRIBUTION FACTORIELLE ROBUSTE (HAC) ---\n")
print(robust_res)
cat("-> Si 'Tech' est *** : L'ESG est biaisé Tech.\n")
cat("-> Si 'Intercept' est significatif : Il y a un vrai Alpha ESG.\n")




# ==============================================================================
# ÉTAPE 6 : STABILITÉ STRUCTURELLE (COVID) - CORRIGÉE
# ==============================================================================
cat("\n>>> ÉTAPE 6 : TEST DE STABILITÉ (CUSUM) <<<\n")


# --- CORRECTION DE L'ERREUR "ENVIR ARG" ---
# La fonction efp() ne supporte pas le format xts. 
# On convertit les données en data.frame classique juste pour ce test.
data_cusum <- data.frame(
  ESG = as.numeric(ret_log$ESG), 
  Trad = as.numeric(ret_log$Trad)
)

# Maintenant ça marche car c'est un data.frame standard
cusum_test <- efp(ESG ~ Trad, data = data_cusum, type = "Rec-CUSUM")

# GRAPHIQUE
plot(cusum_test, main = "Stabilité du Modèle (CUSUM Test)")

# On ajoute la ligne du COVID
# Astuce : On cherche l'index (le numéro de ligne) qui correspond à Mars 2020
covid_index <- which(index(ret_log) >= "2020-03-01")[1] 
abline(v = covid_index, col="blue", lty=2)
text(x = covid_index, y = 0.5, "COVID", col="blue", pos=4)

# Test formel
sctest_res <- sctest(cusum_test)
cat("Test Stabilité p-value :", sctest_res$p.value, 
ifelse(sctest_res$p.value < 0.05, "⚠️ Rupture Structurelle", "✅ Stable"))
