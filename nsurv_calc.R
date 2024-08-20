library(tidyr)
library(dplyr)

nsurv_calc <- function(df,
                       treatment,
                       replicate,
                       time,
                       survivor_columns,
                       dead_columns,
                       initial_population) {
  out_df <- df |>  # Rowwise summation across provided columns
    mutate(total_surv = rowSums(df |>
                                  select(all_of(survivor_columns)), na.rm = TRUE)) |>
    mutate(total_dead = rowSums(df |>
                                  select(all_of(dead_columns)), na.rm = TRUE)) |>
    arrange(treatment, replicate, time)  # order for cumulative calculations
  
  out_df <- out_df |>  # group and calculate cumulative values
    group_by(.data[[treatment]], .data[[replicate]]) |>
    mutate(known_survivors = rev(cummax(rev(total_surv)))) |>
    mutate(sum_of_dead = cumsum(total_dead))
  
  if (exists('initial_population')) {
    # use an initial population column
    out_df <- out_df |>
      mutate(max_possible_dead =
               pmax(0, .data[[initial_population]] - known_survivors)) |>
      mutate(expected_dead = pmin(max_possible_dead, sum_of_dead)) |>
      mutate(Nsurv = .data[[initial_population]] - expected_dead)
  } else {
    # Otherwise use max known_survivors
    out_df <- out_df |>
      mutate(max_possible_dead = max(known_survivors) - known_survivors) |>
      mutate(expected_dead = pmin(max_possible_dead, sum_of_dead)) |>
      mutate(Nsurv = max(known_survivors) - expected_dead)
  }
  
  return(
    out_df |>
      select(
        -expected_dead,
        -max_possible_dead,
        -sum_of_dead,
        -known_survivors,
        -total_surv,
        -total_dead
      )
  )
}

## Example

df <- read.csv('example.csv')

survivor_columns = c('living.larvae', 'living.pupae', 'emerged.midges')

dead_columns = c('dead.larvae', 'dead.pupae')

print(
  nsurv_calc(
    df |> select(
      c(
        treatment = conc,
        replicate,
        time,
        living.larvae,
        living.pupae,
        emerged.midges,
        dead.larvae,
        dead.pupae,
        initial.population
      )
    ),
    treatment = 'treatment',
    replicate = 'replicate',
    time = 'time',
    survivor_columns = c('living.larvae', 'living.pupae', 'emerged.midges'),
    dead_columns = c('dead.larvae', 'dead.pupae'),
    initial_population = 'initial.population'
  )
  ,
  n = 96
)
