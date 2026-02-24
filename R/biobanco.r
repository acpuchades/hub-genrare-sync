library(readxl)
library(dplyr)
library(here)
library(janitor)
library(lubridate)
library(tidyr)


data_dir <- here("data", "20260221")
biobanco_muestras_ela_path <- file.path(data_dir, "Muestras ELA.xlsx")
biobanco_muestras_pals_path <- file.path(data_dir, "SB24-11 Registre activitat mostres.xlsx")

biobanco_muestras_ela <- readxl::read_excel(biobanco_muestras_ela_path, skip=1) |>
  janitor::clean_names() |>
  mutate(
    across(nhc, as.integer),
    across(starts_with("fecha_"), lubridate::dmy)
  ) |>
  drop_na(nhc)

biobanco_muestras_pals <- read_excel(biobanco_muestras_pals_path, skip = 2) |>
  clean_names() |>
  rename(
    volum_ul = volum_ml,
    ubicacio_1 = ubicacio_14,
    ubicacio_2 = x15,
    ubicacio_3 = x16
  ) |>
  fill(
    codi_origen, nhc_dni, data_obtencio, codi_biobanc, muestra,
    data_entrada_bb, hora_entrada, estat_mostra,
    hora_inici_processament, tipus_de_mostra, hora_congelacio,
    conservacio, data_sortida, ubicacio_1
  ) |>
  drop_na(codi_mostra) |>
  group_by(codi_biobanc, muestra) |>
  fill(observacions) |>
  ungroup() |>
  filter(estat_mostra != "NO ARRIBA") |>
  mutate(
    across(starts_with("data_") & -data_sortida,
           ~convert_to_date(.x, character_fun = dmy, string_conversion_failure = "warning")),
    data_entrada_bb = data_entrada_bb + as.numeric(hora_entrada) * ddays(1),
    data_congelacio = data_entrada_bb + as.numeric(hora_congelacio) * ddays(1),
    ) |>
  select(-hora_entrada, -hora_congelacio) |>
  mutate(
    index_mostra = (
      (data_obtencio - min(data_obtencio, na.rm = TRUE)) / dmonths(6)
    ) |> round() |> pmax(0),
    .by = nhc_dni
  )
