RNLachrom <- function(R1, R2=Rb, Rb, I, C, e,
                      interpolate = TRUE, nm = seq(300, 700, 1)) {
  
  #SEE WEBER CONTRAST IN SONKE JOHNSEN BOOK
  #https://academic.oup.com/beheco/article/29/2/273/4560781
  photo=1
  dependent <- FALSE
  nphoto = ncol(C) - 1
  photo1<-photo
  contrast <- "Weber"
  noise <- TRUE
  v<-NA
  n<-NA
  model = "log"
  
  #warnings
  ifelse(photo1 != nphoto, yes = warning("Argument 'C' has more than one sensitivity curve. Using only the first one.", 
                                         call. = FALSE), no = "")
  ifelse(photo1 != length(e) && noise == T, 
         yes = warning("Argument 'e' has more than one value. Using only the first one.", 
                       call. = FALSE), no = "")
  ifelse(photo1 != length(n) && noise == F, 
         yes = warning("Argument 'n' has a number of parameters different than 'photo'.", 
                       call. = FALSE), no = "")
  ifelse(any(ncol(I) > 2), yes = warning("'I' argument with more than two columns. Only the first two will be used.", 
                                         call. = FALSE), no = "")
  ifelse(any(ncol(Rb) > 2), yes = warning("'Rb' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  ifelse(any(ncol(R2) > 2), yes = warning("'R2' argument with more than two columns. Only the first two will be used.", 
                                          call. = FALSE), no = "")
  maxR <- apply(data.frame(R1[, 2:ncol(R1)]), 2, max)
  ifelse(any(maxR <= 1) && max(Rb[, 2]) > 1 || any(maxR > 1) && max(Rb[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R1' and 'Rb' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                       call. = FALSE), no = "")
  ifelse(any(maxR <= 1) && max(R2[, 2]) > 1 || any(maxR > 1) && max(R2[, 2]) <= 1,
         yes = warning("There seems to be a problem with input files. 'R1' and 'R2' must be in the same scale. Both must be either in percentage (0-100%) or proportion (0-1).", 
                       call. = FALSE), no = "")
  if (noise==FALSE&length(n)>1) {
    message("More than one estimated noise. The model assumes that noise refers to the most common receptor.")
  }

  internal <- function(model, photo, R1, R2=Rb, Rb, I, C,
                       e, interpolate, nm) {
    
    #relative photon catches
    S1<-vector(length=photo1)
    S2<-vector(length=photo1)
    for (i in 1:photo1) {
      S1[[i]] <- Qr(I = I, R = R1, Rb = Rb, C = C[, c(1, i+1)], interpolate = interpolate, nm = nm)
      S2[[i]] <- Qr(I = I, R = R2, Rb = Rb, C = C[, c(1, i+1)], interpolate = interpolate, nm = nm)
    }
    
    S1.Qr<-S1
    S2.Qr<-S2
    
    if (model == "log") {
      S1 <- log(S1)
      S2 <- log(S2)
    }
    
    #noise
    if (dependent == FALSE) {
        noise_values<-noise_e(noise = noise, e = e, v = v, n = n)
        noise_values<-noise_values[[1]]
    }
    
    #deltafi=ln[qi(spec1)/qi(spec2)] = ln(1+deltaqi/qi(spec2)] = deltaqi/qi(spec2)

    if (model == "log") {
      if (contrast=="Weber") {
        delta_f<- S1 - S2
      }
      if (contrast=="Michelson") {
        stop("Michelson contrast not implemented")
      }
    }
    
    if (model == "linear") {
      if (contrast=="Weber") {
        delta_f<- (S1 - S2)/S2
      }
      if (contrast=="Michelson") {
        stop("Michelson contrast not implemented")
      }
    }
    
    #deltaS = |deltafi/w|
    
    if (contrast=="Weber") {
      deltaS<-abs(delta_f)/noise_values
    }
    
    if (contrast=="Michelson") {
      stop("Michelson contrast not implemented")
      #deltaS<-abs(delta_f)/noise_values
    }
    
    r<-c(noise_values, S1.Qr, S2.Qr, S1, S2, deltaS)
    
    return(r)
  }
  
  #apply function for several spectra
  n.spectra <- ncol(R1) - 1
  R1.list <- vector(length = n.spectra, "list")
  for (i in 1:n.spectra) {
    R1.list[[i]] <- data.frame(R1[, 1], R1[, i + 1])
  }
  r <- sapply(X = R1.list, FUN = internal, model = model, photo = photo, 
              I = I, R2=R2, Rb = Rb, C = C, e = e, interpolate = interpolate, nm = nm)
  r <- as.data.frame(t(r))
  
  e_names<-"e1"
  Qr_1names<-"Qr1_R1"
  Qr_2names<-"Qr1_R2"
  E_1names<-"E1_R1"
  E_2names<-"E1_R2"

  colnames(r) <- c(e_names, Qr_1names, Qr_2names, E_1names, E_2names, "deltaS")
  rownames(r) <- names(R1)[2:ncol(R1)]
  
  test1<- r[,c(E_1names, E_2names)]
  ifelse(any(test1 < 0), yes = warning("Photoreceptor output < 0.", 
                                       call. = FALSE), no = "")
  
  class(r)<-c("colourvision", "data.frame")
  attr(r, "model_name") <- "Receptor noise limited model - Achromatic"
  attr(r, "model_input_output") <- model
  attr(r, "n_photor_types") <- photo1
  attr(r, "R2") <- R2
  attr(r, "Rb") <- Rb
  attr(r, "I") <- I
  attr(r, "C") <- C
  attr(r, "noise calculated") <- noise
  attr(r, "v") <- v
  attr(r, "n") <- n
  attr(r, "Interpolate") <- interpolate
  attr(r, "nm") <- nm
  attr(r, "contrast") <- contrast
  
  return(r)
  
}
