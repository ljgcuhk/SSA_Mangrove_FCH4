#This code is for SSA decompose and reconstruct multiscale FCH4 during DOY100-300 2017, HK mangrove
#For more details about SSA and its implement, readers are recommended for the monograph:
#Golyandina, N., Korobeynikov, A., & Zhigljavsky, A. (2018). Singular spectrum analysis with R. Berlin, Heidelberg: Spring-Verlag.

#install.packages("Rssa")
library(Rssa)

setwd("D:/Multiscale_Interaction_CH4/SSA/200d")

# Read fragmented original time series
# To reduce edge effects, artificial extension of time series is used (P+2P, i.e. 200d + 2*200d).
Flux <- read.csv(file = "Flux_200day_ext.csv")

# Read 30-min CH4 fluxes
F_CH4 <- Flux$CH4

# Decompose the original time series, with N/L = 2.5 (N = 28944, then window length = 11578).
# Default window length is 1/2 N.
# If SSA decomposition can not be implemented on fragmented time series, add "force.decompose = FALSE".
s_CH4 <- ssa(F_CH4, L = 11578, svd.method = "eigen", kind = "1d-ssa", force.decompose = FALSE)


# Gap fill the fragmented 30-min NEE fluxes using the iterative reconstruction method.
# igapfill() requires a SSA-decomposed structure as input.
# Here we only do gap-filling using the leading low-frequency component.
# It is time-consuming to complete the gap-filling
gf_CH4_1 <- igapfill(s_CH4, groups = list(c(1)),  base = "original", fill = NULL, trace = TRUE)

#SSA decompose the gap-filled time series
s_CH4_gf1 <- ssa(gf_CH4_1, L = 11578)

# Plot of eigenvalues
plot(s_CH4_gf1)
# Plot of eigenvectors 
# Here the 'idx' argument denotes the indices of vectors of interest 
# Because the time series was gap-filled with the leading components, and therefore the contribution of the first eigenvector is increased
plot(s_CH4_gf1, type = "vectors", idx=1:12)
# Calculate the w-correlation matrix using the first 30 components. 
# Here the 'groups' argument as usual denotes the grouping used. 
# This can be used to identify "signal"and "noises". 
w <- wcor(s_CH4_gf1, groups = as.list(1:30)) 
plot(w)
# Plot of elementary reconstructed series 
# Here the 'groups' argument specifies the grouping 
plot(s_CH4_gf1, type = "series", groups = as.list(1:6))


########################To determine thr########################################
# While the choice of the window length is well supported by the SSA theory, 
# the procedure for choosing the eigentriples for grouping is much less formal (Golyandina, 2019).

# There are two methods for grouping: one is based on eigentriples(Golyandina et al, 2018 monograph)
# Another is based on finding components with similar frequency characteristics(this study uses this one)

#To determine the argument "thr" for the function "grouping.auto":
#To understand what is a reasonable value of the threshold, we first choose an 
#arbitrary small threshold to draw the plot of component contributions in the chosen
#frequency ranges; we reorder the components by their contributions. 


s_thr_0 <- ssa(gf_CH4_1, L = 100)
s_thr_1 <- ssa(gf_CH4_1, L = 150)
s_thr_2 <- ssa(gf_CH4_1, L = 300)

#Frequency range corresponds to wavelet analysis (2n * measurement)


g0 <- grouping.auto(s_thr_0, base = "series", groups = 1:100, 
                    freq.bins = list(Hourly = c(1/4, 1/2)),
                    threshold = 0.1)

g1 <- grouping.auto(s_thr_1, base = "series", groups = 1:150,
                    freq.bins = list(Daily = c(1/(1.3*48), 1/8)),
                    threshold = 0.1)

g2 <- grouping.auto(s_thr_2, base = "series", groups = 1:300,
                    freq.bins = list(Multiday = c(1/(21.3*48), 1/(2.7*48))),
                    threshold = 0.1)

#combine
#s_thr_CH4 <- ssa(gf_CH4_1, L = 300)

#g3 <- grouping.auto(s_thr_CH4, base = "series", groups = 1:300, 
                          #freq.bins = list(Hourly = c(1/4, 1/2), Daily = c(1/(1.3*48), 1/8), 
                                           #Multiday = c(1/(21.3*48), 1/(2.7*48))),
                          #threshold = 0.1)

#plot(g3, order = TRUE, type = "b")
                  
                    
contrib0 <- attr(g0, "contributions")[, 1]
contrib1 <- attr(g1, "contributions")[, 1]
contrib2 <- attr(g2, "contributions")[, 1]

plot(g0, order = TRUE, type = "b")
plot(g1, order = TRUE, type = "b")
plot(g2, order = TRUE, type = "b")

#For example for hourly time scales
#Visualization shows that the threshold should be between the contributions of the 55th and 56th components.
#We choose the contribution of the 55th component, which is approximately equal to 0.15, as a new threshold

print(sort(contrib0, decreasing = TRUE)[55]) #0.15
print(sort(contrib1, decreasing = TRUE)[42]) #0.26
print(sort(contrib2, decreasing = TRUE)[4]) #0.57


################################RECONSTRUCTION###############################

#hourly
group.hourly_CH4 <- grouping.auto(s_CH4_gf1, base = "series", groups = 1:11578, 
                                  freq.bins = list(c(1/4, 1/2)),
                                  threshold = 0.15, grouping.method = "pgram")

r_hourly_CH4 <- reconstruct(s_CH4_gf1, groups = group.hourly_CH4)


#daily
group.daily_CH4 <- grouping.auto(s_CH4_gf1, base = "series", groups = 1:11578, 
                                 freq.bins = list(c(1/(1.3*48), 1/8)),
                                 threshold = 0.26, grouping.method = "pgram")

r_daily_CH4 <- reconstruct(s_CH4_gf1, groups = group.daily_CH4)


#Multiday
group.Multiday_CH4 <- grouping.auto(s_CH4_gf1, base = "series", groups = 1:28032, 
                                    freq.bins = list(c(1/(21.3*48), 1/(2.7*48))),
                                    threshold = 0.57, grouping.method = "pgram")

r_Multiday_CH4 <- reconstruct(s_CH4_gf1, groups = group.Multiday_CH4)

#write.csv(r_hourly_CH4$F1,"CH4_hourly_200d_thr0.14.csv")
#write.csv(r_daily_CH4$F1,"CH4_daily_200d_thr0.23.csv")
#write.csv(r_Multiday_CH4$F1,"CH4_Multiday_200d_thr0.52.csv")
#save.image(file = "SSA_Flux_200day_complete.Rdata")




