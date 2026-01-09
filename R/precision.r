library(dplyr)
library(here)

pals_info_path <- here("data", "Pacientes P-ALS 2026-01-09.xlsx")

pals_patients <- readxl::read_excel(
  pals_info_path, sheet = "Pacientes", na = c("", "N/A")
) |>
  janitor::clean_names()

pals_ecas <- readxl::read_excel(pals_info_path, sheet = "ECAS") |>
  janitor::clean_names() |>
  mutate(
    across(nombrar:resultado_ec & -starts_with("resultado_"),
           ~.x |> as.character() |> na_if("-") |> as.numeric())
  ) |>
  arrange(pals_id, fecha)

pals_alsftdq <- readxl::read_excel(pals_info_path, sheet = "ALS-FTD-Q") |>
  janitor::clean_names() |>
  arrange(pals_id, fecha)

pals_eq5 <- readxl::read_excel(pals_info_path, sheet = "EQ-5D-5L") |>
  janitor::clean_names() |>
  arrange(pals_id, fecha)
