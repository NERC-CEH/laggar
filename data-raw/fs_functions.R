# exploration plots and tables ####
model_method_plot <- function(locs, mnem, poly){
  dat <- filter(locs, Mnemonic == mnem)
  dat_S <- filter(dat, DBH <= 20)
  dat_L <- filter(dat, DBH >= 100)

  dat_S1 <- dat_S[128,]
  dat_S1_buff1 <- st_buffer(dat_S1, 25) %>%
    st_intersection(poly)
  dat_S1_buff2 <- st_buffer(dat_S1, 50) %>%
    st_intersection(poly)
  dat_S1_buff3 <- st_buffer(dat_S1, 75) %>%
    st_intersection(poly)
  dat_S1_buff4 <- st_buffer(dat_S1, 100) %>%
    st_intersection(poly)
  dat_S1_buff5 <- st_buffer(dat_S1, 125) %>%
    st_intersection(poly)
  dat_S1_buff6 <- st_buffer(dat_S1, 150) %>%
    st_intersection(poly)
  dat_S1_buff7 <- st_buffer(dat_S1, 175) %>%
    st_intersection(poly)
  dat_S1_buff8 <- st_buffer(dat_S1, 200) %>%
    st_intersection(poly)

  pal <- viridis::viridis(10)[3:10]

  s1 <-ggplot() +
      theme_void() +
      geom_sf(data = poly) +
      geom_sf(data = dat_S1_buff8, fill = pal[1]) +
      geom_sf(data = dat_S1_buff7, fill = pal[2]) +
      geom_sf(data = dat_S1_buff6, fill = pal[3]) +
      geom_sf(data = dat_S1_buff5, fill = pal[4]) +
      geom_sf(data = dat_S1_buff4, fill = pal[5]) +
      geom_sf(data = dat_S1_buff3, fill = pal[6]) +
      geom_sf(data = dat_S1_buff2, fill = pal[7]) +
      geom_sf(data = dat_S1_buff1, fill = pal[8]) +
      geom_sf(data = dat_L) +
      geom_sf(data = dat_S1, colour = "red")

  # sum y in buffers
  sumy_buff <- lapply(seq(25,200,25), function(x) sum_within_dist(dat_S1, dat_L, poly, x))
  sumy_buff_df <- data.frame(num_points = diff(c(0, sapply(sumy_buff, "[[", 1))),
                             buffer_area = 1e-4*diff(c(0,sapply(sumy_buff, "[[", 2))),
                             distance = seq(25,200,25))
  p1 <- ggplot(sumy_buff_df, aes(x = distance, y = num_points)) +
    geom_bar(stat = "identity", fill = rev(pal)) +
    labs(x = "Distance (m)", y = "Number of adult trees") +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic()
  p2 <- ggplot(sumy_buff_df, aes(x = distance, y = num_points/buffer_area)) +
    geom_bar(stat = "identity", fill = rev(pal)) +
    labs(x = "Distance (m)", y = "Number of adult trees per hectare") +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic()
  allplot <- s1/(p1+p2) + plot_layout(heights = c(1.4,1))
  return(allplot)
}


# Modelling setup ####
sum_within_dist <- function(x, y, poly, dist){
  buff <- st_buffer(x, dist)
  y_in_buff <- st_contains(buff, y)
  sum_y_in_buff <- sapply(y_in_buff, length)
  suppressWarnings(area_in_buff <- st_area(st_intersection(buff, poly)))
  list(sum_y_in_buff,area_in_buff)
}

get_pts_buffer <- function(x, y, poly, mindist, maxdist, incdist){
  dist <- seq(mindist,maxdist,incdist)
  sumdist <- lapply(dist, sum_within_dist, x = x, y = y, poly = poly)
  num_points <- sapply(sumdist, "[[", 1)
  buffer_area <- sapply(sumdist, "[[", 2)

  # convert from cumulative
  num_points <- t(apply(num_points, 1, function(x) diff(c(0,x))))
  buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
  buffer_area <- 1e-4*buffer_area #convert to hectares

  return(list(num_points = num_points,
              buffer_area = buffer_area))
}

create_poly <- function(maxX = 1000, maxY = 500,
                        edgedist = 0){
  data.frame(
    lon = c(edgedist, maxX - edgedist, maxX- edgedist, edgedist),
    lat = c(edgedist, edgedist, maxY- edgedist, maxY- edgedist),
    var = c(1, 1, 1, 1)
  ) %>%
    sfheaders::sf_polygon(x = "lon", y = "lat"
    )
}

create_grid <- function(maxX = 1000, maxY = 500,
                        distance_fromedge = 20,
                        increment = 10){
  xseq <- seq(distance_fromedge,maxX - distance_fromedge,increment)
  yseq <- seq(distance_fromedge,maxY - distance_fromedge,increment)
  grid_pts <- data.frame(PX = rep(xseq, length(yseq)),
                         PY = rep(yseq, each = length(xseq))) %>%
    st_as_sf(coords = c("PX","PY"))
  return(grid_pts)
}


get_grid_data <- function(loc_data, grid_pts, mnem,
                          dist_seq,
                          poly, juv_cutoff, adult_cutoff){
  dat <- filter(loc_data, Mnemonic == mnem)
  dat_S <- filter(dat, DBH <= juv_cutoff)
  countrees <- sapply(st_is_within_distance(grid_pts, dat_S, 5), length)

  dat_L <- filter(dat, DBH > adult_cutoff)
  sumy_buff <- lapply(dist_seq, function(x) sum_within_dist(grid_pts, dat_L, poly, x))
  num_points <- sapply(sumy_buff, "[[", 1)
  buffer_area <- sapply(sumy_buff, "[[", 2)

  # convert from cumulative
  num_points <- t(apply(num_points, 1, function(x) diff(c(0,x))))
  buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
  buffer_area <- 1e-4*buffer_area #convert to hectares

  trees_pha <- num_points/buffer_area

  return(trees_pha)
}


get_grid_data_bydbhclass <- function(loc_data, grid_pts, mnem,
                                     dist_seq, dbh_class_lower_bounds,
                                     poly){
  dat <- filter(loc_data, Mnemonic == mnem)
  dat_S <- filter(dat, DBH <= 20)
  countrees <- sapply(st_is_within_distance(grid_pts, dat_S, 5), length)


  sumy_grid <- mapply(function(lower,upper) {
    dat_L <- filter(dat, DBH > lower & DBH <= upper)
    sumy_buff <- lapply(dist_seq, function(x) sum_within_dist(grid_pts, dat_L, poly, x))
    num_points <- sapply(sumy_buff, "[[", 1)
    buffer_area <- sapply(sumy_buff, "[[", 2)

    # convert from cumulative
    num_points <- t(apply(num_points, 1, function(x) diff(c(0,x))))
    buffer_area <- t(apply(buffer_area, 1, function(x) diff(c(0,x))))
    buffer_area <- 1e-4*buffer_area #convert to hectares

    trees_pha <- num_points/buffer_area

    return(trees_pha)
  }, lower = dbh_class_lower_bounds,
  upper = c(dbh_class_lower_bounds[-1],1e6), SIMPLIFY = FALSE)


  trees_pha_wide <- do.call(cbind, sumy_grid)
  return(trees_pha_wide)
}

grid_dat_ct <- function(mnem, loc_data, grid_pts){
  dat <- filter(loc_data, Mnemonic == mnem)
  dat_S <- filter(dat, DBH <= 20)
  y <- sapply(st_is_within_distance(grid_pts, dat_S, 5), length)
  return(y)
}

# model runs  ####
grid_mods_bydbhclass <- function(loc_data, grid_pts, mnem,
                                 dist_seq, dbh_class_lower_bounds,
                                 poly, data_list, dbhmax, ...){
  mnem <- mnem[dbhmax[mnem] > min(dbh_class_lower_bounds)]
  if(any(!mnem %in% names(data_list)))
    stop("All species must be in data list")
  mod_list <- lapply(mnem, function(mnem){
    y <- grid_dat_ct(mnem, loc_data, grid_pts)

    x <- data_list[[mnem]]

    treedist <- matrix(rep(dist_seq, length(dbh_class_lower_bounds)),
                       nrow=nrow(x), ncol = ncol(x),
                       byrow = TRUE)
    treeDBH <- matrix(rep(dbh_class_lower_bounds, each = length(dist_seq)),
                      nrow=nrow(x), ncol = ncol(x),
                      byrow = TRUE)
    data_list <- list(y = y, x = x, treedist = treedist, treeDBH = treeDBH)


    gamfit <- gam(y ~ te(treedist, treeDBH, by = x, k =c(15,6)), data = data_list,
                  family = "nb", method = "REML", ...)
    return(gamfit)
  })
  names(mod_list) <- mnem
  return(mod_list)
}

grid_mods_noDBH <- function(loc_data, grid_pts, mnem,
                            dist_seq, adult_dbhcutoff,
                            poly, data_list, dbhmax, ...){
  mnem <- mnem[dbhmax[mnem] > adult_dbhcutoff]
  if(any(!mnem %in% names(data_list)))
    stop("All species must be in data list")
  mod_list <- lapply(mnem, function(mnem){
    y <- grid_dat_ct(mnem, loc_data, grid_pts)

    x <- data_list[[mnem]]

    treedist <- matrix(dist_seq,
                       nrow=nrow(x), ncol = ncol(x),
                       byrow = TRUE)
    data_list <- list(y = y, x = x, treedist = treedist)


    gamfit <- gam(y ~ te(treedist, by = x, bs = "tp", k = 15),
                  data = data_list,
                  family = "nb", method = "REML", ...)
    return(gamfit)
  })
  names(mod_list) <- mnem
  return(mod_list)
}


grid_mods_mv <- function(loc_data, grid_pts, mnem,
                         dist_seq, dbh_class_lower_bounds,
                         poly, data_list, dbhmax, ...){
  # grid_data_bydbhclass <- lapply(mnem,get_grid_data_bydbhclass,
  #                                loc_data = loc_data, grid_pts = grid_pts,
  #                                dist_seq = dist_seq,
  #                                dbh_class_lower_bounds= dbh_class_lower_bounds,
  #                                poly = poly)
  mnem <- mnem[dbhmax[mnem] > min(dbh_class_lower_bounds)]
  if(any(!mnem %in% names(data_list)))
    stop("All species must be in data list")
  grid_data_bydbhclass <- data_list[mnem]
  ngrpts <- nrow(grid_pts)
  grid_data_x <- do.call(rbind, grid_data_bydbhclass)
  grid_data_y <- unlist(lapply(mnem, grid_dat_ct,
                               loc_data = loc_data, grid_pts = grid_pts))

  treedist <- matrix(rep(dist_seq, length(dbh_class_lower_bounds)),
                     nrow=nrow(grid_data_x), ncol = ncol(grid_data_x),
                     byrow = TRUE)
  treeDBH <- matrix(rep(dbh_class_lower_bounds, each = length(dist_seq)),
                    nrow=nrow(grid_data_x), ncol = ncol(grid_data_x),
                    byrow = TRUE)
  treeSpec <- matrix(rep(rep(as.integer(as.factor(mnem)),
                             each = ngrpts), ncol(grid_data_x)),
                     nrow=nrow(grid_data_x), ncol = ncol(grid_data_x))
  data_list <- list(y = grid_data_y, x = grid_data_x,
                    treedist = treedist, treeDBH = treeDBH, treeSpec = treeSpec)

  grid_coords <- st_coordinates(grid_pts)

  gresp <- cbind(y = grid_data_y,
                 spec = as.integer(as.factor(rep(mnem, each = ngrpts))))

  data_list <- list(y = gresp, x = grid_data_x,
                    species = as.factor(rep(mnem, each = ngrpts)),
                    PX = rep(grid_coords[,1],length(mnem)),
                    PY = rep(grid_coords[,2],length(mnem)),
                    treedist = treedist, treeDBH = treeDBH, treeSpec = treeSpec)


  mvgamfit <- bam(y ~ species + s(PX, PY) +
                    te(treedist, treeDBH,
                       bs = c("tp","tp"), by = x, k =c(15,6)) +
                    te(treedist, treeDBH, treeSpec,
                       bs = c("tp","tp","re"), by = x, k =c(15,6,5)),
                  data = data_list,
                  family = gfam(replicate(length(mnem),"nb",
                                          simplify = FALSE)),
                  method = "REML")

  return(mvgamfit)

}


grid_mods_mv_noDBH <- function(loc_data, grid_pts, mnem,
                               dist_seq, data_list, dbhmax,
                               adult_dbhcutoff, ...){
  # grid_data_bydbhclass <- lapply(mnem,get_grid_data_bydbhclass,
  #                                loc_data = loc_data, grid_pts = grid_pts,
  #                                dist_seq = dist_seq,
  #                                dbh_class_lower_bounds= dbh_class_lower_bounds,
  #                                poly = poly)
  mnem <- mnem[dbhmax[mnem] > adult_dbhcutoff]
  if(any(!mnem %in% names(data_list)))
    stop("All species must be in data list")
  grid_data <- data_list[mnem]
  ngrpts <- nrow(grid_pts)
  grid_data_x <- do.call(rbind, grid_data)
  grid_data_y <- unlist(lapply(mnem, grid_dat_ct,
                               loc_data = loc_data, grid_pts = grid_pts))

  treedist <- matrix(dist_seq,
                     nrow=nrow(grid_data_x), ncol = ncol(grid_data_x),
                     byrow = TRUE)
  treeSpec <- matrix(rep(rep(as.integer(as.factor(mnem)),
                             each = ngrpts), ncol(grid_data_x)),
                     nrow=nrow(grid_data_x), ncol = ncol(grid_data_x))

  grid_coords <- st_coordinates(grid_pts)

  gresp <- cbind(y = grid_data_y,
                 spec = as.integer(as.factor(rep(mnem, each = ngrpts))))

  data_list <- list(y = gresp, x = grid_data_x,
                    species = as.factor(rep(mnem, each = ngrpts)),
                    PX = rep(grid_coords[,1],length(mnem)),
                    PY = rep(grid_coords[,2],length(mnem)),
                    treedist = treedist, treeSpec = treeSpec)


  mvgamfit <- bam(y ~ species + s(PX, PY) +
                    te(treedist,
                       bs = "tp", by = x, k =15) +
                    te(treedist, treeSpec,
                       bs = c("tp","re"), by = x, k =c(8,5)),
                  data = data_list,
                  family = gfam(replicate(length(mnem),"nb",
                                          simplify = FALSE)),
                  method = "REML")

  return(mvgamfit)

}


# Model checks and plots ####
mods_residcheck <- function(modlist){
  simresid_list <- lapply(modlist, simulateResiduals)
  t1 <- lapply(simresid_list, testZeroInflation, plot = FALSE)
  t2 <- lapply(simresid_list, testUniformity, plot = FALSE)
  t3 <- lapply(simresid_list, testQuantiles, plot = FALSE)
  return(list(ZeroInflation = t1,
              Uniformity = t2,
              Quantiles = t3))
}

resid_signif_hist <- function(object){
  dat <- data.frame(mnem = names(object$ZeroInflation),
                    ZeroInflation = unlist(sapply(object$ZeroInflation,
                                            "[","p.value")),
                    Uniformity = unlist(sapply(object$Uniformity,
                                            "[","p.value")),
                    Quantiles = unlist(sapply(object$Quantiles,
                                            "[","p.value"))
  ) %>%
    tidyr::pivot_longer(ZeroInflation:Quantiles, names_to = "test",
                        values_to = "p.value") %>%
    mutate(alpha05 = ifelse(p.value < 0.05, "<0.05",">0.05"))
  pl <- ggplot(dat, aes(x = p.value, fill = alpha05)) +
    geom_histogram(binwidth = 0.05) +
    scale_fill_manual(values = c("#D55E00","grey35")) +
    facet_wrap(~test) +
    theme(legend.position = "none")
  return(pl)

}

reduce_ggplot_size <- function(object){
  object <-  ggplot2::ggplotGrob(object)
  object <- ggpubr::as_ggplot(object)
  return(object)
}

reduce_mgcViz_size <- function(object){
  object <-  ggplot2::ggplotGrob(object$ggObj)
  object <- ggpubr::as_ggplot(object)
  return(object)
}

mods_plots <- function(mnem, spec_file, modlist){
  pl_list_nb_viz <- lapply(mnem, function(i){
    pl1 <- plot(sm(modlist[[i]],1)) +
      l_fitRaster() +
      l_fitContour() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Distance (m)", y = "DBH (mm)",
           title = paste0(paste(filter(spec_file, Mnemonic == i)[,c("Genus", "Species")],
                                collapse = " "), " (",i,")")) +
      theme(legend.position = "none")
    pl1 <- reduce_mgcViz_size(pl1)
    return(pl1)
  })
  return(pl_list_nb_viz)
}
mods_plots_1d <- function(mnem, spec_file, modlist){
  pl_list_nb_viz <- lapply(mnem, function(i){
    pl1 <- plot(sm(modlist[[i]],1)) +
      l_fitLine() + l_ciLine() +
      # scale_x_continuous(expand = c(0,0)) +
      # scale_y_continuous(expand = c(0,0)) +
      labs(x = "Distance (m)", #y = "DBH (mm)",
           title = paste0(paste(filter(spec_file, Mnemonic == i)[,c("Genus", "Species")],
                                collapse = " "), " (",i,")")) +
      theme(legend.position = "none")
    pl1 <- reduce_mgcViz_size(pl1)
    return(pl1)
  })
  return(pl_list_nb_viz)
}

mv_mods_plots <- function(object, mnem, dist_seq, dbh_seq, maxdbh){
  mnem <- mnem[maxdbh[mnem] > min(dbh_seq)]
  maxdbh <- maxdbh[mnem]
  pred_df <- data.frame(species = rep(mnem, each = length(dist_seq)*length(dbh_seq)),
                        treedist = rep(rep(dist_seq, length(dbh_seq)),length(mnem)),
                        treeDBH = rep(rep(dbh_seq,each = length(dist_seq)), length(mnem)),
                        PX = 500, PY = 250, x = 1,
                        maxdbh = rep(maxdbh, each = length(dist_seq)*length(dbh_seq)))
  pred_df$treeSpec <- as.integer(as.factor(pred_df$species))

  mod_predict <- predict(object, newdata = pred_df,
                         type = "iterms")
  pred_df <- mutate(pred_df,
                    pp = mod_predict[,3],
                    ps = mod_predict[,4]) %>%
    mutate(total = pp + ps,
           overMaxDBH = ifelse(treeDBH > maxdbh, 0.5, 1))

  pl_list <- lapply(mnem, function(mnem) {
    pl1 <- pred_df %>% filter(species == mnem) %>%
      ggplot(aes(x = treedist, y = treeDBH)) +
      geom_tile(aes(fill = total, alpha = overMaxDBH)) +
      geom_contour(aes(z = total), col = "black", bins = 10) +
      scale_fill_viridis_c(begin = 0.2) +
      ggtitle(mnem) +
      scale_alpha_identity() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))+
      theme(legend.position = "none")
    pl1 <- reduce_ggplot_size(pl1)
    return(pl1)
  })
  return(pl_list)

}

mv_mods_plots_noDBH <- function(object){
  mnem <- as.character(unique(object$model$species))
  dist_seq <- c(unique(object$model$treedist))
  pred_df <- data.frame(species = rep(mnem, each = length(dist_seq)),
                        treedist = rep(dist_seq,length(mnem)),
                        PX = 500, PY = 250, x = 1)
  pred_df$treeSpec <- as.integer(as.factor(pred_df$species))

  mod_predict <- predict(object, newdata = pred_df,
                         se.fit = TRUE, type = "iterms")
  pred_df <- mutate(pred_df,
                    pp = mod_predict$fit[,3],
                    ps = mod_predict$fit[,4],
                    ps_se = mod_predict$se.fit[,4]) %>%
    mutate(upper = ps + ps_se, lower = ps-ps_se)

  pl_list <- lapply(mnem, function(mnem) {
    pl1 <- pred_df %>% filter(species == mnem) %>%
      ggplot(aes(x = treedist, y = ps)) +
      geom_ribbon(aes(ymin = lower, ymax  = upper), fill = "grey")+
      geom_line() +
      ggtitle(mnem)
    pl1 <- reduce_ggplot_size(pl1)
    return(pl1)
  })
  return(pl_list)

}
