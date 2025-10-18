library(dplyr)
library(tidyr)

biobanco_muestras_path <- here::here("data", "Biobanco 20250701.xlsx")

biobanco_muestras <- readxl::read_excel(biobanco_muestras_path, skip=1) |>
  janitor::clean_names() |>
  mutate(
    across(nhc, as.integer),
    across(starts_with("fecha_"), lubridate::dmy)
  ) |>
  drop_na(nhc)
