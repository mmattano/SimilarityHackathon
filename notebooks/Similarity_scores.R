library(MsCoreUtils)
library(tidyr)


lib_MS2 <- Spectra("C:/Users/rfm848/Desktop/EuBic_2023/ex_spectra_stds_NEG_scans_50.mgf", 
                   source = MsBackendMgf())

exp_MS2 <- Spectra("C:/Users/rfm848/Desktop/EuBic_2023/pseudo_query_neg_v1.mgf", 
                   source = MsBackendMgf())






table_combination <- data.frame(m = c(rep(0,4), rep(0.5,4), rep(1,4), rep(2,4)),
                                n= rep(c(0,0.5,1,2),4)
                                )


cosine <- list()

for (tolerance in c(0.01, 0.005)) {
  for (i in 1:nrow(table_combination)-1){
    
    m <- table_combination[i,"m"]
    n <- table_combination[i,"n"]
    
    scores_cosine <- compareSpectra(exp_MS2, 
                                    lib_MS2,
                                    m = m,
                                    n = n,
                                    FUN = ndotproduct, 
                                    tolerance = tolerance)
    scores_cosine <- scores_cosine %>% 
                     as_tibble %>% 
                     setNames(lib_MS2$FEATURE_ID) %>% 
                     mutate(id_query = exp_MS2$scanIndex)
    
    write.csv(scores_cosine, 
              row.names = FALSE,
              paste0("C:/Users/rfm848/Desktop/EuBic_2023/scores_cosine_", "m_", m, "_n_", n, "_tolerance_", tolerance, ".csv"))
    
    cosine[[paste(m, n, tolerance, sep = "_")]] <- scores_cosine
  }
  

}


nspecanglescore <- list()

for (tolerance in c(0.01, 0.005)) {
  for (i in 1:nrow(table_combination)){
    
    m <- table_combination[i,"m"]
    n <- table_combination[i,"n"]
    
    scores_nspecangle <- compareSpectra(exp_MS2, 
                                        lib_MS2,
                                        m = m,
                                        n = n,
                                        FUN = nspectraangle, 
                                        tolerance = tolerance)
    scores_nspecangle <- scores_nspecangle %>% 
                         as_tibble %>% 
                         setNames(lib_MS2$FEATURE_ID) %>% 
                         mutate(id_query = exp_MS2$scanIndex)
    
    write.csv(scores_nspecangle, 
              row.names = FALSE,
              paste0("C:/Users/rfm848/Desktop/EuBic_2023/scores_nspecangle_", "m_", m, "_n_", n, "_tolerance_", tolerance, ".csv"))
    
    nspecanglescore[[paste(m, n, tolerance, sep = "_")]] <- scores_nspecangle
  }
  
  
}



# 
# scores_nspecangle <- compareSpectra(exp_MS2, 
#                                     lib_MS2, 
#                                     FUN = nspectraangle, 
#                                     tolerance = 0.01)




navdistscore <- list()

for (tolerance in c(0.01, 0.005)) {
  for (i in 1:nrow(table_combination)){
    
    m <- table_combination[i,"m"]
    n <- table_combination[i,"n"]
    
    scores_navdis <- compareSpectra(exp_MS2, 
                                    lib_MS2,
                                    m = m,
                                    n = n,
                                    FUN = navdist, 
                                    tolerance = tolerance)
    scores_navdis <- scores_navdis %>% 
                     as_tibble %>% 
                     setNames(lib_MS2$FEATURE_ID) %>% 
                     mutate(id_query = exp_MS2$scanIndex)
    
    write.csv(scores_navdis, 
              row.names = FALSE,
              paste0("C:/Users/rfm848/Desktop/EuBic_2023/scores_navdis_", "m_", m, "_n_", n, "_tolerance_", tolerance, ".csv"))
    
    navdistscore[[paste(m, n, tolerance, sep = "_")]] <- scores_navdis
  }
  
  
}





# scores_navdis <- compareSpectra(exp_MS2, 
#                                 lib_MS2, 
#                                 FUN = navdist, 
#                                 tolerance = 0.01)



neuclideanscore <- list()

for (tolerance in c(0.01, 0.005)) {
  for (i in 1:nrow(table_combination)){
    
    m <- table_combination[i,"m"]
    n <- table_combination[i,"n"]
    
    scores_neuclidean <- compareSpectra(exp_MS2, 
                                        lib_MS2,
                                        m = m,
                                        n = n,
                                        FUN = neuclidean, 
                                        tolerance = tolerance)
    scores_neuclidean <- scores_neuclidean %>% 
                         as_tibble %>% 
                         setNames(lib_MS2$FEATURE_ID) %>% 
                         mutate(id_query = exp_MS2$scanIndex)
    
    write.csv(scores_neuclidean, 
              row.names = FALSE,
              paste0("C:/Users/rfm848/Desktop/EuBic_2023/scores_neuclidean_", "m_", m, "_n_", n, "_tolerance_", tolerance, ".csv"))
    
    neuclideanscore[[paste(m, n, tolerance, sep = "_")]] <- scores_neuclidean
  }

}

# scores_neuclidean <- compareSpectra(exp_MS2, 
#                                     lib_MS2, 
#                                     FUN = neuclidean, 
#                                     tolerance = 0.01)






