library(here)

consents_path <- here('data', 'GENRARE 2026-01-09.xlsx')

consents <- readxl::read_excel(consents_path) |>
  janitor::clean_names() |>
  dplyr::mutate(
    dplyr::across(c(comparte_con_grupos:consiente_contacto, interesado_ensayos), is.na)
  )
