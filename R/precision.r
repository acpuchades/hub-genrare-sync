pals_info_path <- here::here("data", "PRECISIONALS 20251013.xlsx")

pals_patients <- readxl::read_excel(pals_info_path, sheet = "Pacientes") |>
  janitor::clean_names()

pals_ecas <- readxl::read_excel(pals_info_path, sheet = "ECAS") |>
  janitor::clean_names()

pals_eq5 <- readxl::read_excel(pals_info_path, sheet = "EQ-5D-5L") |>
  janitor::clean_names()
