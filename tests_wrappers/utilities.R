library(sn)  
library(purrr)
library(pracma)
library(MASS)

# Plot functions  
plot_field2D = function(x, y, z, title = "", levels_contour_by = 0.2, 
                        path_img, width, height, xtick, ytick, colkey = TRUE, 
                        mar = c(5.1, 4.1, 4.1, 2.1), number.figures = 1, mfrow = c(1,1), 
                        font.size = 3){
  
  cols <- hcl.colors(100, "YlOrRd", rev = FALSE)
  
  # Show 
  par(mar = mar, mfrow = mfrow)
  for(i in 1:number.figures){
    if(number.figures > 1){
      contour.list = list(levels = round(range(z[[i]])*seq(0,1,by=levels_contour_by),3), 
                          drawlabels = TRUE, labcex = font.size)
      image2D(x = x, y = y, z = z[[i]],
              xlab = "", ylab = "", xaxt="n", yaxt="n",  
              col = cols, colkey = colkey,  
              contour = contour.list)
      title(main = title[[i]])
    }
    
    else{
      contour.list = list(levels = round(range(z)*seq(0,1,by=levels_contour_by),3), 
                          drawlabels = TRUE, labcex = font.size)
      image2D(x = x, y = y, z = z,
              xlab = "", ylab = "", xaxt="n", yaxt="n",  
              col = cols, colkey = colkey,  
              contour = contour.list)
      title(main = title)
    }
    
  }
  if(!is.null(xtick))
    axis(side=1, at=xtick)
  if(!is.null(xtick))
    axis(side=2, at=ytick)
  
  # Save 
  pdf(path_img, width = width, height = height)
  par(mar = mar, mfrow = mfrow)
  for(i in 1:number.figures){
    if(number.figures > 1){
      contour.list = list(levels = round(range(z[[i]])*seq(0,1,by=levels_contour_by),3), 
                          drawlabels = TRUE, labcex = font.size)
      image2D(x = x, y = y, z = z[[i]],
              xlab = "", ylab = "", xaxt="n", yaxt="n",  
              col = cols, colkey = colkey,  
              contour = contour.list)
      title(main = title[[i]])
    }
    else{
      contour.list = list(levels = round(range(z)*seq(0,1,by=levels_contour_by),3), 
                          drawlabels = TRUE, labcex = font.size)
      image2D(x = x, y = y, z = z,
              xlab = "", ylab = "", xaxt="n", yaxt="n",  
              col = cols, colkey = colkey,  
              contour = contour.list)
      title(main = title)
    }
  }
  
  if(!is.null(xtick))
    axis(side=1, at=xtick)
  if(!is.null(xtick))
    axis(side=2, at=ytick)
  
  dev.off()
  
  
}

plot_boxplot_RvsCpp = function(matrix_R, matrix_Cpp, 
                               seq_x, title = "", font.labels = 1, xlab = "",
                               ylim = NULL, names = NULL, true = NULL,
                               shift = FALSE){
  
  myred <- rgb(255, 0, 0, max = 255, alpha = 125)
  myblue <- rgb(0, 0, 255, max = 255, alpha = 125)
  myred_dark <- rgb(255, 0, 0, max = 255, alpha = 255)
  myblue_dark <- rgb(0, 0, 255, max = 255, alpha = 255)
  mypink <- rgb(255, 102, 204, max = 255, alpha = 255)
  mylightblue <- rgb(102, 204, 255, max = 255, alpha = 255)
  mygrey_dark = rgb(105, 105, 105, max = 255, alpha = 255)
  
  dist_axis = 1
  dist_title = 1
  
  if(shift) {
    width_box = 0.3
    shift_box = 0.2
  }
  if(!shift){
    width_box = 0.6
    shift_box = 0.0
  }
  
  
  par(mar=c(6,4,4,4))
  boxplot(matrix_R, cex.axis = font.labels, 
          ylim = ylim, names = names, 
          col = myblue,
          medcol = myblue_dark,
          boxcol = myblue_dark,
          whiskcol = myblue_dark,
          staplecol = myblue_dark,
          outcol = myblue,
          outpch = 19, 
          at = seq(1,length(seq_x))-shift_box,
          boxwex = width_box,
          xaxt="n")
  axis(1, axTicks(1), labels=F, cex.axis = font.labels)
  mtext(names, 
        1, dist_axis, at=axTicks(1), cex = font.labels)
  mtext(xlab,  cex = font.labels, side = 1, line = dist_title)
  
  if(title != "")
    title(title, cex = 2)
  if(!is.null(true))
    abline(h=true, col = "darkgreen", lwd=2.5, lty=1)
  
  boxplot(matrix_Cpp, cex.axis = font.labels,
          ylim = ylim, names = names, add = TRUE,
          col = myred,
          medcol = myred_dark,
          boxcol = myred_dark,
          whiskcol = myred_dark,
          staplecol = myred_dark,
          outcol = myred,
          outpch = 19, 
          at = seq(1,length(seq_x))+shift_box,
          boxwex = width_box,
          xaxt="n")
  axis(1, axTicks(1), labels=F) 
  mtext(names, 
        1, dist_axis, at=axTicks(1), cex=font.labels)
  
  
  if(title != "")
    title(title, cex = 2)

  n_string = "n"
  if(xlab == "log10(Nodes)"){
    n_string = "N"
  }
  
}


# Metrics 
RMSE = function(x, y){
  return(sqrt(mean((x-y)^2)))
}

SpRMSE = function(fitted_vec, true_vec){
  SpRMSE_vec = c()
  n = dim(fitted_vec)[1]
  M = dim(fitted_vec)[2]
  for(i in 1:n){
    sum = 0
    for(m in 1:M){
        fn_hat = fitted_vec[, m]
        fn_true = true_vec[, m]
        sum = sum + (fn_hat[i] - fn_true[i])^2
      }
      SpRMSE_vec = c(SpRMSE_vec, sqrt(sum/M))
  }
  return(SpRMSE_vec)
}

# C-shaped functions 
area.horseshoe = function(r, w, l){
  options(digits=16)
  return(3/2*pi*w^2 + pi/2*r^2 - 2*pi*w*r + 2*l*(w-r))
}

a_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    pi/4 + p[1]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -pi/4 -p[1]
    }else{
      if(p[1] < 0){
        -0.5*atan(p[2]/p[1])
      }
    }
  }
}

d_fun <- function(p){
  
  if(p[1]>= 0 && p[2] > 0){
    -0.5 + p[2]
  }else{
    if(p[1]>= 0 && p[2] <= 0){
      -0.5 - p[2]
    }else{
      if(p[1] < 0){
        sqrt(p[1]^2 + p[2]^2) - 0.5
      }
    }
  }
}

z <- function(p){
  a_fun(p) + d_fun(p)^2
}

z_xy <- function(x,y){z(c(x,y))}

generate_mesh_Cnetwork = function(eps = 0.5){
  # eps: vertical dimension 
  x = c(0.,1)
  y = c(0.,eps) 
  vertices.coarse = cbind(expand.grid(x,y)[,1], expand.grid(x,y)[,2])
  edges.coarse = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)
  mode(edges.coarse) = "integer"
  
  mesh.coarse = create.mesh.1.5D(vertices.coarse, edges.coarse)     
  mesh = refine.mesh.1.5D(mesh.coarse, delta=0.0125) 
  
  nodes = mesh$nodes
  edges = mesh$edges
  mode(edges) = "integer"
  neigh = mesh$neigh
  boundary = as.integer(mesh$nodesmarkers)
  
  # Adjust neighbors
  nedges = nrow(mesh$nodes)
  neigh = Matrix(0, nrow = nedges, ncol=nedges)
  for(i in 1:nrow(mesh$edges)){ neigh[ mesh$edges[i,1], mesh$edges[i,2] ] = 1; neigh[ mesh$edges[i,2], mesh$edges[i,1] ] = 1}
  neigh <- as(neigh, "TsparseMatrix")
  neigh.TMP = matrix(0, nrow=length(neigh@i), ncol=3)
  for(k in 1:length(neigh@i)){
    neigh.TMP[k, ] = c( neigh@i[k], neigh@j[k], neigh@x[k] )
  }
  neigh.TMP[1:nrow(neigh.TMP), 1:2] = neigh.TMP[1:nrow(neigh.TMP), 1:2] + matrix(1, nrow=nrow(neigh.TMP), ncol=2)
  neigh = neigh.TMP
  storage.mode(neigh) <- "integer"
  
  return(list(nodes=nodes, edges=edges, neigh=neigh, boundary=boundary))
}

# Data generation functions

data_creation_hetero_fun = function(n, mean.function, std.function, beta_true = NULL,
                                    locations, nodes, covariates = NULL, seed){
  # heteroscedastic data starting from true mean & std fields functions
  #         y = rnorm(n, mean = field_mean, sd = field_std)
  
  if( min(min(std.function(locations))) <= 0 ){
    cat("\n Error: negative standard deviation \n")
    return(NULL)
  }
  
  set.seed(seed)
  data = rnorm(n, mean = mean.function(locations), sd = std.function(locations))
  if(!is.null(covariates)){
    q = dim(covariates)[2]
    if(is.null(beta_true)){
      beta_true = seq(1,q)
    }
    data <- data + covariates%*%beta_true
  }
  
  # Field true
  fn_true <-  mean.function(locations) + qnorm(alpha) * std.function(locations)
  f_true <-  mean.function(nodes) + qnorm(alpha) * std.function(nodes)
  
  return(list(data = data, beta_true = beta_true, fn_true = fn_true, f_true = f_true))
}

data_creation_hetero_grf = function(nxx, nodes, locations, covariates = NULL, beta_true = NULL, seed, 
                                    tau2_mean = 1, nu_mean = 2, rho_mean = 0.3, tau2_sd = 0.3, nu_sd = 2, 
                                    rho_sd = 0.6, anisotropy = NULL){
  
  # heteroscedastic data starting from true mean & std fields Gaussian Random Field
  #         y = rnorm(n, mean = field_mean, sd = field_std)
  
  library(geoR)
  
  seed.fixed = 10 # a fixed seed for the grf (we do not want that mu,sigma change with 
  # the simulation seed)
  
  n = dim(locations)[1]
  nnodes = dim(nodes)[1]
  
  if(!is.null(anisotropy)){    # anisotropic case
    aniso.pars = c(anisotropy$angle, anisotropy$intensity)
  } else{                      # isotropic case
    aniso.pars = c(0, 1) 
  }
  
  # Mean field
  rho_mean = rho_mean / (2*sqrt(nu_mean))
  cov.pars_mean = c(tau2_mean, rho_mean)
  
  set.seed(seed.fixed)
  mean_loc = grf(n, grid = locations, sqrt(n), sqrt(n), xlims = c(0, 1), ylims = c(0, 1),
                 nsim = 1, cov.model = "matern",
                 cov.pars = cov.pars_mean , 
                 kappa = nu_mean, 
                 aniso.pars = aniso.pars,
                 messages = FALSE)$data   # nugget = 0, lambda = 1, aniso.pars,
  set.seed(seed.fixed)
  mean_nodes = grf(nnodes, grid = nodes, sqrt(nnodes), sqrt(nnodes), xlims = c(0, 1), ylims = c(0, 1),
                   nsim = 1, cov.model = "matern",
                   cov.pars = cov.pars_mean , 
                   kappa = nu_mean, 
                   aniso.pars = aniso.pars,
                   messages = FALSE)$data   # nugget = 0, lambda = 1, aniso.pars,
  
  
  # Std field
  rho_sd = rho_sd / (2*sqrt(nu_sd))
  cov.pars_sd = c(tau2_sd, rho_sd)
  
  set.seed(seed.fixed)
  sd_loc = grf(n, grid = locations, sqrt(n), sqrt(n), xlims = c(0, 1), ylims = c(0, 1),
               nsim = 1, cov.model = "matern",
               cov.pars = cov.pars_sd , 
               kappa = nu_sd, 
               aniso.pars = aniso.pars, 
               messages = FALSE)$data  # nugget = 0, lambda = 1, aniso.pars,
  set.seed(seed.fixed)
  sd_nodes = grf(nnodes, grid = nodes, sqrt(nnodes), sqrt(nnodes), xlims = c(0, 1), ylims = c(0, 1),
                 nsim = 1, cov.model = "matern",
                 cov.pars = cov.pars_sd , 
                 kappa = nu_sd, 
                 aniso.pars = aniso.pars, 
                 messages = FALSE)$data  # nugget = 0, lambda = 1, aniso.pars,
  
  if(min(sd_loc) <= 0){
    cat("\n Adjusting negative values of the std gaussian field... \n")
    sd_loc[which(sd_loc <= 0)] = sd_loc[which(sd_loc <= 0)] - 1.2*sd_loc[which(sd_loc <= 0)]
  }
  if(min(sd_nodes) <= 0){  # if there are negative values 
    cat("\n Adjusting negative values of the std gaussian field... \n")
    sd_nodes[which(sd_nodes <= 0)] = sd_nodes[which(sd_nodes <= 0)] - 1.2*sd_nodes[which(sd_nodes <= 0)]
  }
  
  set.seed(seed)   # this seed is simulation dependent 
  data = rnorm(n, mean = mean_loc, sd = sd_loc)
  if(!is.null(covariates)){
    q = dim(covariates)[2]
    if(is.null(beta_true)){
      beta_true = seq(1,q)
    }
    data = data + covariates%*%beta_true  
  }
  
  # Field true
  fn_true <- mean_loc + qnorm(alpha) * sd_loc
  f_true <- mean_nodes + qnorm(alpha) * sd_nodes
  
  return(list(data = data, beta_true = beta_true, 
              mean.field_nodes = mean_nodes, std.field_nodes = sd_nodes,
              mean.field_loc = mean_loc, std.field_loc = sd_loc, 
              fn_true = fn_true, f_true = f_true))
} 

data_creation_skewed = function(data_generation, alpha, locations, nodes, 
                                beta_true = NULL, covariates = NULL, seed){
  
  # skewed data generated as a deterministic function of the position + skewed-gaussian noise 
  #         y = data_generation(x,y) + rsn 
  
  n = dim(locations)[1]
  y_true <- data_generation(locations[,1], locations[,2]) 
  if(!is.null(covariates)) {
    q = dim(covariates)[2]
    if(is.null(beta_true)){
      beta_true = seq(1,q)
    }
    y_true = y_true + covariates%*%beta_true  
  }
  
  # Add (skewed zero-mean) noise 
  library(sn)  # Skewed normal Distribution
  xi_ = 4 
  omega_ = 0.05*(max(y_true)-min(y_true))
  alpha_noise = 5
  delta_ = alpha_noise / (sqrt(1+alpha_noise^2))
  scale_noise = 1
  
  set.seed(seed)
  skewed_t <- scale_noise*rsn(n, xi = xi_, omega = omega_, alpha = alpha_noise) 
  
  noise_true_mean = scale_noise*(xi_ + omega_*delta_*sqrt(2/pi))
  noise_true_sd = scale_noise*sqrt(omega_*(1-2*delta_^2/pi))
  
  skewed_t <- skewed_t - noise_true_mean   # zero mean noise   
  noise_true_quantile = scale_noise*qsn(alpha, xi = xi_, omega = omega_, alpha = alpha_noise) - noise_true_mean
  
  # Generate simulated data (pointwise at nodes) 
  data <- y_true + skewed_t
  
  # True field
  fn_true <- noise_true_quantile + data_generation(locations[,1], locations[,2]) 
  f_true <- noise_true_quantile + data_generation(nodes[,1], nodes[,2]) 
  
  return(list(data = data, beta_true = beta_true, fn_true = fn_true, f_true = f_true, 
              noise_true_mean = 0, noise_true_sd = noise_true_sd))
}  

data_creation_areal = function(mesh, incidence_matrix, alpha, beta_true, m){
  n = dim(incidence_matrix)[1]
  integration_nodes = data.frame(matrix(nrow = n, ncol = 11))
  names(integration_nodes) = c("T1","T2","N1","N2","N3","N4", "xmin", "xmax", "ymin", "ymax", "label")
  
  tri = mesh$triangles
  nodi_plot = NULL
  for(i in 1:n)
  {
    tri_used = which(incidence_matrix[i,] == 1)
    integration_nodes$T1[i] = tri_used[1]
    integration_nodes$T2[i] = tri_used[2]
    nodes_used = unique(c(tri[tri_used[1],],tri[tri_used[2],]))
    integration_nodes$N1[i] = nodes_used[1]
    integration_nodes$N2[i] = nodes_used[2]
    integration_nodes$N3[i] = nodes_used[3]
    integration_nodes$N4[i] = nodes_used[4]
    integration_nodes$label[i] = i
    xvec = c(nodes[integration_nodes$N1[i],1],nodes[integration_nodes$N2[i],1],nodes[integration_nodes$N3[i],1],nodes[integration_nodes$N4[i],1])
    yvec = c(nodes[integration_nodes$N1[i],2],nodes[integration_nodes$N2[i],2],nodes[integration_nodes$N3[i],2],nodes[integration_nodes$N4[i],2])
    integration_nodes$xmin[i] = min(xvec)
    integration_nodes$xmax[i] = max(xvec)
    integration_nodes$ymin[i] = min(yvec)
    integration_nodes$ymax[i] = max(yvec)
    
    nodi_tmp = rbind(mesh$nodes[integration_nodes$N1[i],],mesh$nodes[integration_nodes$N2[i],],mesh$nodes[integration_nodes$N3[i],],mesh$nodes[integration_nodes$N4[i],])
    nodi_tmp = cbind(nodi_tmp,rep(i,4))
    nodi_plot = rbind(nodi_plot,nodi_tmp)
  }
  
  area_domains = rep(0.04, n)
  
  sol_integrand = matrix(nrow = n , ncol = 1)
  for(i in 1:n){
    v = c(z(mesh$nodes[integration_nodes$N1[i],]), z(mesh$nodes[integration_nodes$N2[i],]), z(mesh$nodes[integration_nodes$N3[i],]), z(mesh$nodes[integration_nodes$N4[i],]))

    res = integral2(fun = z_xy, xmin = integration_nodes$xmin[i], xmax = integration_nodes$xmax[i] ,
                    ymin = integration_nodes$ymin[i], ymax = integration_nodes$ymax[i], 
                    vectorized = FALSE)
    sol_integrand[i] = (res$Q)/area_domains[i]

    
  }
  
  nnodes = dim(mesh$nodes)[1]
  sol_exact <- numeric(nnodes)
  
  set.seed(3)
  covariates=matrix(0,nrow= n ,ncol=1)
  covariates[,1]=rbeta(n ,shape1=2,shape2=2)  # sampling covariates from beta distr
  
  y_true <- sol_integrand + covariates%*%beta_true
  
  # Add (skewed zero-mean) noise 
  xi_ = 4 
  omega_ = 0.1*(max(y_true)-min(y_true))
  alpha_noise = 5
  delta_ = alpha_noise / (sqrt(1+alpha_noise^2))
  scale_noise = 1
  
  noise_true_mean = scale_noise*(xi_ + omega_*delta_*sqrt(2/pi))
  noise_true_sd = scale_noise*sqrt(omega_*(1-2*delta_^2/pi))
  
  for(i in 1:nnodes){
    sol_exact[i] <- z(mesh$nodes[i,])
  }
  
  set.seed(21*m)
  skewed_t <- scale_noise*rsn(n, xi = xi_, omega = omega_, alpha = alpha_noise) 
  
  skewed_t <- skewed_t - noise_true_mean   # zero mean noise   
  noise_true_quantile = scale_noise*qsn(alpha, xi = xi_, omega = omega_, alpha = alpha_noise) - noise_true_mean
  
  # Generate simulated data (pointwise at nodes) 
  data <- y_true + skewed_t
  
  # True field
  f_true <- noise_true_quantile + sol_exact 
  
  fn_true <- noise_true_quantile + sol_integrand 
  
  
  return(list(data = data, beta_true = beta_true, covariates = covariates, fn_true = fn_true, f_true = f_true, 
              noise_true_sd = noise_true_sd))
}  

# Others 
almost_equal = function(double1, double2){
  return( abs(double1 - double2) < 1e-11 )
}

