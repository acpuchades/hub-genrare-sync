library(dplyr)
library(tibble)
library(tidyr)

data_dir <- here::here('data', '20260221')
ufela_db_path <- file.path(data_dir, 'formulario.sqlite')
ufela_db = DBI::dbConnect(RSQLite::SQLite(), ufela_db_path)

ufela_pacientes = DBI::dbGetQuery(ufela_db, "SELECT * FROM pacientes") |>
  rows_update(tibble(pid = "43c5588a-fea1-11ef-823b-3b027ebef8ae", fecha_nacimiento = "04-02-1984"), by = "pid") |>
  rows_update(tibble(pid = "c2ab61d6-fddd-11ef-beed-edbe42534d56", fecha_nacimiento = "29-11-1946"), by = "pid") |>
  rows_update(tibble(pid = "dea43e4a-596c-11eb-9950-239b631e5519", exitus = "Sí", fecha_exitus = "29-07-2024"), by = "pid") |>
  mutate(
    across(nhc, as.integer),
    across(starts_with("fecha_"), lubridate::dmy),
    across(ends_with("_datetime"), lubridate::ymd_hms),
    across(exitus, ~case_match(.x, "Sí" ~ TRUE, "No" ~ FALSE)),
  )

ufela_clinica = DBI::dbGetQuery(ufela_db, "SELECT * FROM datos_clinicos") |>
  rows_update(tibble(pid = "9342fe7c-d949-11e9-842a-ebf9c1d8fdac", fecha_visita_datos_clinicos = "26-01-2016"), by = "pid") |>
  rows_update(tibble(pid = "cd45be9a-fddb-11ef-beed-edbe42534d56", fecha_inicio_riluzol = "11-04-2023"), by = "pid") |>
  rows_update(tibble(pid = "ebf7ca94-7397-11ec-803e-e338475c84d8", fecha_diagnostico_ELA = "03-01-2022"), by = "pid") |>
  mutate(
    across(starts_with("fecha_"), lubridate::dmy),
    across(c(starts_with("historia_familiar_"), matches("estudio_genetico_(c9|sod1)"), riluzol), ~case_match(.x, "Sí" ~ TRUE, "No" ~ FALSE)),
    across(c(starts_with("resultado_estudio_"), fumador), ~na_if(.x, "NS/NC")),
    estudio_genetico_atxn2 = str_detect(estudio_genetico_otro, "ATA?XN2"),
    estudio_genetico_ar = str_detect(estudio_genetico_otro, "KENNEDY|AR"),
    estudio_genetico_tardbp = str_detect(estudio_genetico_otro, "TARDBP"),
    estudio_genetico_fus = str_detect(estudio_genetico_otro, "FUS"),
    estudio_genetico_unc13a = str_detect(estudio_genetico_otro, "UNCA?13A"),
    numero_repeticiones_atxn2_a1 = estudio_genetico_otro |>
      str_extract("ATA?XN2 *([^@]+)", group=1) |>
      str_extract("([0-9]+) */ *([0-9]+)", group=1) |>
      as.integer(),
    numero_repeticiones_atxn2_a2 = estudio_genetico_otro |>
      str_extract("ATA?XN2 *([^@]+)", group=1) |>
      str_extract("([0-9]+) */ *([0-9]+)", group=2) |>
      as.integer(),
    resultado_estudio_atxn2 = case_when(
      is.na(estudio_genetico_atxn2) ~ NA,
      numero_repeticiones_atxn2_a1 > 26 ~ "Alterado",
      numero_repeticiones_atxn2_a2 > 26 ~ "Alterado",
      numero_repeticiones_atxn2_a1 <= 26 &
        numero_repeticiones_atxn2_a2 <= 26 ~ "Normal",
      TRUE ~ if_else(estudio_genetico_otro |>
        str_extract("ATA?XN2 ([^@]+)", group = 1) |>
        str_detect("INTERMEDIO?"),
        "Alterado", "Normal"
      )
    ),
    resultado_estudio_tardbp = case_when(
      is.na(estudio_genetico_tardbp) ~ NA,
      TRUE ~ if_else(estudio_genetico_otro |>
        str_extract("TARDBP ([^@]+)", group = 1) |>
        str_detect("ALTERA(T|DO)|HOMOCIGOTO|PATOGENICO?|POSITI(U|VO)"),
        "Alterado", "Normal"
      )
    ),
    resultado_estudio_fus = case_when(
      is.na(estudio_genetico_fus) ~ NA,
      TRUE ~ if_else(estudio_genetico_otro |>
        str_extract("FUS ([^@]+)", group = 1) |>
        str_detect("ALTERA(T|DO)|HOMOCIGOTO|PATOGENICO?|POSITI(U|VO)"),
        "Alterado", "Normal"
      )
    ),
    resultado_estudio_unc13a = case_when(
      is.na(estudio_genetico_unc13a) ~ NA,
      TRUE ~ if_else(estudio_genetico_otro |>
        str_extract("UNCA?13A ([^@]+)", group = 1) |>
        str_detect("ALTERA(T|DO)|HOMOCIGOTO|PATOGENICO?|POSITI(U|VO)"),
        "Alterado", "Normal"
      )
    ),
  )

ufela_alsfrs = DBI::dbGetQuery(ufela_db, "SELECT * FROM esc_val_ela") |>
  rename(fecha_visita = fecha_visita_esc_val_ela) |>
  mutate(
    across(starts_with("fecha_"), lubridate::dmy),
    across(
      lenguaje:insuficiencia_respiratoria,
      ~ if_else(
        .x |> na_if("NS/NC") |> as.integer() |> between(0, 4),
        as.integer(.x), NA_integer_
      )
    ),
    across(
      kings,
      ~ if_else(
        .x |> na_if("NS/NC") |> as.integer() |> between(0, 4),
        as.integer(.x), NA_integer_
      )
    ),
    cortar = case_when(
      is.na(cortar_con_peg) & !is.na(cortar_sin_peg) ~ cortar_sin_peg,
      !is.na(cortar_con_peg) & is.na(cortar_sin_peg) ~ cortar_con_peg,
      cortar_con_peg == cortar_sin_peg ~ cortar_sin_peg,
      cortar_con_peg > 0 & cortar_sin_peg == 0 ~ cortar_con_peg,
      cortar_con_peg == 0 & cortar_sin_peg > 0 ~ cortar_sin_peg,
      TRUE ~ NA_integer_
    ),
    total_bulbar = lenguaje + salivacion + deglucion,
    total_motor_fino = escritura + cortar + vestido,
    total_motor_grosero = cama + caminar + subir_escaleras,
    total_respiratorio = disnea + ortopnea + insuficiencia_respiratoria,
    across(total, ~ coalesce(
      total_bulbar + total_motor_fino + total_motor_grosero + total_respiratorio,
      ifelse(.x |> na_if("NS/NC") |> as.integer() |> between(0, 48), as.integer(.x), NA_integer_)
    )),
    mitos = {
      movement <- (caminar <= 1) | (vestido <= 1)
      swallowing <- deglucion <= 1
      communication <- (lenguaje <= 1) & (escritura <= 1)
      breathing <- (disnea <= 1) | (insuficiencia_respiratoria <= 2)
      movement + swallowing + communication + breathing
    },
  ) |>
  arrange(pid, fecha_visita)

ufela_respi <- DBI::dbGetQuery(ufela_db, "SELECT * FROM fun_res") |>
  rename(fecha_visita = fecha_visita_fun_res) |>
  rows_update(tibble(id = "2c24ac7e-544b-11e8-b5ee-37c4bbc2f8c2", fecha_visita = "2018-05-11"), by = "id") |>
  rows_update(tibble(id = "5acbae18-e612-11eb-a8da-cb9b998413f1", fecha_visita = "2021-07-13"), by = "id") |>
  rows_update(tibble(id = "efd51e08-e57e-11ef-a972-73e96d46a9dc", fvc_sentado="155"), by = "id") |>
  rows_update(tibble(id = "f6f201c4-1284-11ea-badc-8b5db01dc70c", fvc_sentado_absoluto="5208"), by = "id") |>
  rows_update(tibble(id = "4308f4ec-d272-11eb-b08c-851f5b5d129e", fvc_sentado_absoluto="5420"), by = "id") |>
  mutate(
    across(starts_with("fecha_"), lubridate::dmy),
    across(starts_with("fvc_"), as.integer),
    across(c(pim, pem, pns), as.integer),
    across(c(ph_sangre_arterial, pao2, paco2, hco3, sao2_media, ct90), as.numeric),
    across(portador_vmni, ~ case_match(.x, "Sí" ~ TRUE, "No" ~ FALSE)),
    pao2 = if_else(pao2 |> between(60, 250), pao2, NA_real_),
    paco2 = if_else(paco2 |> between(10, 100), paco2, NA_real_),
    hco3 = if_else(hco3 |> between(10, 45), hco3, NA_real_),
    ph_sangre_arterial = case_when(
      ph_sangre_arterial |> between(700, 800) ~ ph_sangre_arterial / 100,
      ph_sangre_arterial |> between(7, 8) ~ ph_sangre_arterial,
      TRUE ~ NA_real_
    ),
  ) |>
  arrange(pid, fecha_visita)

ufela_nutri <- DBI::dbGetQuery(ufela_db, "SELECT * FROM datos_antro") |>
  rename(fecha_visita = fecha_visita_datos_antro) |>
  mutate(
    across(starts_with("fecha_"), lubridate::dmy),
    across(c(peso, peso_premorbido, estatura, imc_actual), as.numeric),
    across(c(starts_with("suplementacion_nutricional_"), portador_peg), ~case_match(.x, "Sí" ~ TRUE, "No" ~ FALSE)),
    estatura = if_else(estatura |> between(1, 2), estatura * 100, estatura),
    imc_actual = coalesce(round(peso / (estatura/100)^2, digits=1), imc_actual),
  ) |>
  arrange(pid, fecha_visita)

ufela_visitas <- bind_rows(
  ufela_alsfrs |> select(pid, fecha_visita),
  ufela_respi |> select(pid, fecha_visita),
  ufela_nutri |> select(pid, fecha_visita),
) |>
  drop_na() |>
  distinct() |>
  arrange(pid, fecha_visita)

ufela_seguimiento <- ufela_visitas |>
  summarize(fecha_ultima_visita = max(fecha_visita, na.rm = TRUE), .by = pid)

DBI::dbDisconnect(ufela_db)
