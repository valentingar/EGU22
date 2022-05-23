################################
## redoing parameter sweep    ##
## for ZALF production model  ##
## Jochheim et al. 2022       ##
################################
# last update 2022-04-26

# dependencies ------------
library(ConFluxPro)
library(dplyr)
library(DBI)
library(tictoc)
library(ggplot2)
library(ggpubr)
library(future)
library(progressr)

# load data 
zalf_db <- dbConnect(RSQLite::SQLite(), "../Daten/DB/zalf_db")

gasdata <- zalf_db %>% 
  tbl("gasdata") %>% 
  collect() %>%
  rename(x_ppm = CO2,
         site = Plot) %>%
  #mutate(datetime = lubridate::as_datetime(datetime)) %>%
  mutate(datetime = as.POSIXct(datetime, origin = "1970-01-01")) %>%
  mutate(count = ifelse(rep == "surface",3,1)) %>%
  tidyr::uncount(count,.id = "new_rep") %>%
  mutate(rep = ifelse(rep == "surface", as.numeric(new_rep), rep)) %>%
  mutate(gas = "CO2") %>%
  filter(!is.na(datetime),
         !is.na(x_ppm)) %>%
  cfp_gasdata(id_cols = c("datetime", "rep", "site"))

soilphys <-
zalf_db %>% 
  tbl("soilphys_complete") %>%
  collect() %>%
  rename(site = Plot, 
         c_air = rho_air,
         t = Temp) %>%
  mutate(datetime = as.POSIXct(Date, origin = "1970-01-01")) %>%
  filter(!is.na(datetime)) %>%
  mutate(TPS = ifelse(upper == -5 &site == "Kiefer", 0.679, TPS),
         a = ifelse(upper == -5 & site == "Kiefer", 0.745, a),
         b = ifelse(upper == -5 & site == "Kiefer", 1.86, b)) %>%
  complete_soilphys(overwrite = TRUE) %>%
  cfp_soilphys(id_cols = c("datetime", "site", "rep"))

layers_map <- 
  data.frame(site = rep(c("Buche","Kiefer"),each = 4),
             layer = rep(c("HU","M1","M2","M3"),2),
             upper = c(5,0,-8,-30,6,0,-10,-30),
             lower = c(0,-8,-30,-100,0,-10,-30,-100)
  ) %>%
  left_join(data.frame(rep = c("1","2","3")), by = character()) %>%
  cfp_layers_map(id_cols = c("site", "rep"),
                 gas = "CO2", 
                 lowlim = 0, 
                 highlim = 1000, 
                 layer_couple = 0)


chamber <- 
  zalf_db %>%
  tbl("chamber") %>%
  collect() %>%
  mutate(datetime = as.POSIXct(Date,origin = "1970-01-01"))%>% 
  tidyr::extract(col = "gr_plot",
          into = "rep",regex = ".*([0-9])+$",remove = F) %>%
  select(-c(Date, gr_plot)) %>%
  rename(site = Plot)
  

zalf <- 
  cfp_dat(gasdata, 
          soilphys, 
          layers_map)

rm(gasdata)
rm(soilphys)
gc()

future::plan(future::multisession)
tic()
progressr::with_progress(
PROFLUX <- 
  zalf %>%
  pro_flux()
)
toc()
future::plan(future::sequential)

rm(zalf)
gc()

set.seed(630933)
profs_subsample <-
  PROFLUX$profiles %>% 
  group_by(site, rep) %>%
  slice_sample(n = 100) %>%
  pull(prof_id)

PF_subsample <- 
  PROFLUX %>%
  filter(prof_id %in% profs_subsample)

rmap_subsample <-
run_map(PROFLUX,
        params = list("TPS" = c(0.9,1.1),
                   "a" = c(0.9,1.1),
                   "b" = c(0.9,1.1),
                   "topheight" = c(-2,2)),
        type = c("factor", "factor", "factor", "addition"),
        method = "random",
        layers_different = TRUE,
        n_runs = 100
        )

rmap_sobol <- 
  sobol_run_map(rmap_subsample)

i <- 0

future::plan(future::multisession)
progressr::with_progress(
PF_sobol <-
  alternate(PF_subsample,
            f = function(x) {i <<- i+1;message(i);complete_soilphys(x, overwrite = TRUE)},
            run_map = rmap_sobol,
            return_raw = TRUE)
)

future::plan(future::sequential)

EF_sobol <-
  progressr::with_progress(
  cfp_altapply(PF_sobol, efflux)
  )


NRMSE_sobol <-
  progressr::with_progress(
    cfp_altapply(PF_sobol, function(x) error_concentration(x, normer = "sd", param_cols = c("rep", "site")))
  )


sobol_indices <- 
  sobol_calc_indices(EF_sobol, "efflux", id_cols = c("site", "rep"), run_map = rmap_sobol) %>%
  bind_rows(
    sobol_calc_indices(NRMSE_sobol, "NRMSE", id_cols = c("site", "rep"), run_map = rmap_sobol)
  )

sobol_plot <- 
sobol_indices %>%
  tidyr::pivot_longer(c("Si", "ST"),
                      names_to = "index",
                      values_to = "sensitivity") %>%
  ggplot(aes(x=paste0(pmap, param), y = sensitivity, fill = index))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_grid(cols = vars(rep),
             rows = vars(site, effect_param))+
  theme_bw()

dir.create("figs")
ggsave("sobol_indices.pdf", 
       path = "figs/", 
       plot = sobol_plot,
       height = 8, 
       width = 14, 
       units = "cm",
       scale = 3)

dir.create("output")
saveRDS(PF_sobol, "output/PF_sobol.rds")
rm(PF_sobol)
gc()

# parameter sweep ------------------------
#dates with chamber measurements.
chamber_dates <-
chamber %>%
  mutate(Date = lubridate::date(datetime)) %>%
  select(site, rep, Date) %>%
  distinct()

chamber_profs <- 
  PROFLUX$profiles %>%
  mutate(Date = lubridate::date(datetime)) %>%
  right_join(chamber_dates) %>%
  pull(prof_id)

sweep_profs <- 
  unique(c(chamber_profs, profs_subsample))

PF_sweep <- 
  PROFLUX %>% 
  filter(prof_id %in% sweep_profs)

# finding parameters with at least 15 % contribution
# to either NRMSE or efflux
param_selection <-
sobol_indices %>%
  filter(ST > 0.15) %>%
  select(site, rep, pmap, param, param_id) %>%
  distinct()

# fitting 4 most influential parameters per site/rep
param_selection <-
sobol_indices %>%
  group_by(site, rep,param, pmap) %>%
  summarise(ST = max(Si)) %>%
  group_by(site, rep) %>%
  slice_max(ST, n = 4)



rmap_sweep <-
rmap_subsample %>%
  left_join(param_selection %>% select(-any_of("param_id")) %>% mutate(selected = TRUE)) %>%
  mutate(value = ifelse(is.na(selected),ifelse(type == "factor",1,0),value)) %>%
  select(-any_of(c("param_id", "selected"))) %>%
  filter(run_id < 50)
attr(rmap_sweep, "n_runs") <- 50


future::plan(future::multisession)

progressr::with_progress(
PF_alt_sweep <- 
  alternate(PF_sweep,
            function(x) complete_soilphys(x, overwrite = TRUE),
            run_map = rmap_sweep,
            return_raw = TRUE)
)

future::plan(future::sequential)

EF_sweep <- 
  with_progress(
  cfp_altapply(PF_alt_sweep,
               efflux)
  )

NRMSE_sweep <- 
  with_progress(
    cfp_altapply(PF_alt_sweep,
                 function(x) x %>% 
                   filter(prof_id %in% profs_subsample) %>% 
                   error_concentration(normer = "sd", param_cols = c("site", "rep")))
  )

chamber_day <-
  chamber %>%
  mutate(date = lubridate::date(datetime)) %>%
  group_by(site, rep, date) %>%
  summarise(flux_ch = mean(flux_ch, na.rm = TRUE))

NRMSE_ef_sweep <-
EF_sweep %>%
  mutate(date = lubridate::date(datetime)) %>%
  group_by(site, rep, date, run_id) %>%
  summarise(efflux = mean(efflux, na.rm = TRUE)) %>%
  left_join(chamber_day) %>%
  group_by(site, rep, run_id) %>%
  summarise(NRMSE_ef = ConFluxPro:::nrmse(efflux, flux_ch, "sd"))

NRMSE_efsingle_sweep <-
  EF_sweep %>%
  left_join(chamber) %>%
  group_by(site, rep, run_id) %>%
  summarise(NRMSE_ef = ConFluxPro:::nrmse(efflux, flux_ch, "sd"))


sweep_scatter <-
NRMSE_ef_sweep %>% 
  left_join(NRMSE_sweep) %>%
  ggplot(aes(x=NRMSE_ef, y = NRMSE))+
  geom_point()+
  facet_grid(cols = vars(rep),
             rows = vars(site))

ggsave("sweep_scatter.pdf", 
       path = "figs/", 
       plot = sweep_scatter,
       height = 8, 
       width = 14, 
       units = "cm",
       scale = 3)

sweep_per_param <-
NRMSE_efsingle_sweep %>%
  rename(NRMSE = NRMSE_ef) %>%
  mutate(effect_param = "chamber") %>%
  bind_rows(NRMSE_sweep %>% mutate(effect_param = "profile")) %>%
  left_join(rmap_sweep %>% right_join(param_selection)) %>%
  ggplot(aes(x=value, y=NRMSE, col = factor(rep), shape = site))+
  geom_point()+
  facet_grid(rows = vars(site, rep, effect_param), cols = vars(pmap, param),
             scales = "free")+
  theme_bw()

ggsave("sweep_per_param.pdf", 
       path = "figs/", 
       plot = sweep_per_param,
       height = 8, 
       width = 14, 
       units = "cm",
       scale = 3)
