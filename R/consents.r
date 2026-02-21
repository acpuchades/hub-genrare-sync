data_dir <- here::here('data', '20260221')
consents_path <- file.path(data_dir, 'Pacientes GENRARE.xlsx')

consents <- readxl::read_excel(consents_path) |>
  janitor::clean_names() |>
  dplyr::mutate(
    dplyr::across(c(comparte_con_grupos:consiente_contacto, interesado_ensayos), is.na)
  )
