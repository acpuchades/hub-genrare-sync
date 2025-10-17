library(dplyr)

consents_path <- here::here('data', 'GENRARE 20251017.xlsx')

consents <- readxl::read_excel(consents_path) |>
  janitor::clean_names() |>
  mutate(across(c(comparte_con_grupos:consiente_contacto, interesado_ensayos), is.na))
