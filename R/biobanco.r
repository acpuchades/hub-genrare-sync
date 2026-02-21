library(dplyr)
library(tidyr)

data_dir <- here::here("data", "20260221")
biobanco_muestras_path <- file.path(data_dir, "Biobanco.xlsx")

biobanco_muestras <- readxl::read_excel(biobanco_muestras_path, skip=1) |>
  janitor::clean_names() |>
  mutate(
    across(nhc, as.integer),
    across(starts_with("fecha_"), lubridate::dmy)
  ) |>
  drop_na(nhc)
